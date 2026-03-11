"""
dynamics/propagator.py

Propagator: numerically integrates the equations of motion using SciPy's
RK45 solver with J2-perturbed two-body dynamics.

Maneuvers are handled by passing a control_law callable into propagate().
The control law returns a (3,) ECI thrust acceleration at any (t, state),
which is added directly to the equations of motion. This avoids trajectory
splitting and keeps the dynamics self-consistent throughout a burn.

ODE right-hand sides
--------------------
    dynamics_tb_j2          : state only, no thrust         (6-dim)
    dynamics_tb_j2_thrust   : state only, with thrust       (6-dim)
    dynamics_with_stm       : state + STM, no thrust        (42-dim)

The STM ODE does not include thrust because the EKF operates on the
uncontrolled (target) spacecraft or on the chaser between burn updates.
If you need an STM during a continuous burn, the Jacobian would need a
∂u/∂x term — that is left as a future extension.

Control law interface
---------------------
    def my_control_law(t: float, state: np.ndarray) -> np.ndarray:
        '''
        Args:
            t     : elapsed seconds since propagation start
            state : (6,) [rx, ry, rz, vx, vy, vz] ECI, metres / m·s⁻¹
        Returns:
            (3,) thrust acceleration in ECI frame, m/s²
        '''

Pass it to propagate():
    states = propagator.propagate(initial_state, duration,
                                  control_law=my_control_law)
"""

from __future__ import annotations

import numpy as np
from datetime import datetime, timedelta
from typing import Callable, List, Optional

from scipy.integrate import solve_ivp

from dynamics.orbital_state import OrbitalState
from utils.coordinates import MU_EARTH, RE_EARTH, J2


# ── Type alias ────────────────────────────────────────────────────────────────

ControlLaw = Callable[[float, np.ndarray], np.ndarray]


# ── J2 acceleration ───────────────────────────────────────────────────────────

def j2_acceleration(r_vec: np.ndarray, mu: float) -> np.ndarray:
    """
    Acceleration due to the J2 zonal harmonic in the inertial frame.

    Args:
        r_vec : (3,) position vector [x, y, z], metres
        mu    : gravitational parameter, m^3/s^2

    Returns:
        (3,) J2 acceleration, m/s^2
    """
    x, y, z = r_vec
    r      = np.linalg.norm(r_vec)
    factor = (1.5 * J2 * mu * RE_EARTH**2) / r**5
    zx2    = 5.0 * z**2 / r**2
    return np.array([
        factor * x * (zx2 - 1.0),
        factor * y * (zx2 - 1.0),
        factor * z * (zx2 - 3.0),
    ])


# ── Dynamics Jacobian (A matrix) ──────────────────────────────────────────────

def dynamics_jacobian(r: np.ndarray) -> np.ndarray:
    """
    Analytical 6×6 system Jacobian  A = ∂f/∂x  evaluated at position r.

    Used to propagate the STM:  Φ̇ = A · Φ

    Block structure:
        A = [ 0₃  I₃ ]
            [ G   0₃ ]
    where G = ∂a/∂r combines two-body and J2 position gradients.

    Args:
        r : (3,) ECI position vector, metres

    Returns:
        (6,6) system matrix A
    """
    rnorm = np.linalg.norm(r)
    z = r[2]

    I3  = np.eye(3)
    O3  = np.zeros((3, 3))
    rrT = np.outer(r, r)

    dadr_tb = MU_EARTH * (3.0 * rrT / rnorm**5 - I3 / rnorm**3)

    n3    = np.array([0.0, 0.0, 1.0])
    n3n3T = np.outer(n3, n3)
    rn3T  = np.outer(r, n3)
    n3rT  = np.outer(n3, r)

    dadr_j2 = (
        -3.0 * J2 * MU_EARTH * RE_EARTH**2 / 2.0
    ) * (
          I3          / rnorm**5
        + 2.0 * n3n3T / rnorm**5
        - 5.0 * (rrT + z**2 * I3 + 2.0 * z * (rn3T + n3rT)) / rnorm**7
        + 35.0 * z**2 * rrT / rnorm**9
    )

    return np.block([
        [O3,              I3],
        [dadr_tb + dadr_j2, O3],
    ])


# ── ODE right-hand sides ──────────────────────────────────────────────────────

def dynamics_tb_j2(t: float, y: np.ndarray) -> np.ndarray:
    """
    Two-body + J2 equations of motion, no thrust.

    y = [r (3), v (3)]  →  ẏ = [v, a_tb + a_J2]
    """
    r = y[0:3]
    v = y[3:6]
    a = -MU_EARTH * r / np.linalg.norm(r)**3 + j2_acceleration(r, MU_EARTH)
    dydt      = np.zeros(6)
    dydt[0:3] = v
    dydt[3:6] = a
    return dydt


def dynamics_tb_j2_thrust(
    t: float,
    y: np.ndarray,
    control_law: ControlLaw,
) -> np.ndarray:
    """
    Two-body + J2 + continuous thrust equations of motion.

    The thrust acceleration is supplied by control_law(t, y) and added
    directly to the natural acceleration. The control law is responsible
    for returning zero when the thruster is off.

    Args:
        t           : elapsed seconds since propagation start
        y           : (6,) state [r, v]
        control_law : callable (t, state) → (3,) ECI acceleration, m/s²

    Returns:
        (6,) state derivative [v, a_tb + a_J2 + u]
    """
    r = y[0:3]
    v = y[3:6]
    a = (-MU_EARTH * r / np.linalg.norm(r)**3
         + j2_acceleration(r, MU_EARTH)
         + control_law(t, y))
    dydt      = np.zeros(6)
    dydt[0:3] = v
    dydt[3:6] = a
    
    return dydt


def dynamics_with_stm(t: float, y: np.ndarray) -> np.ndarray:
    """
    Two-body + J2 equations of motion augmented with Φ̇ = A·Φ  (no thrust).

    Augmented state: y = [r (3), v (3), Φ.flatten() (36)]  →  42-dim ODE.
    Used by the EKF predict step.
    """
    r   = y[0:3]
    v   = y[3:6]
    Phi = y[6:].reshape((6, 6))

    a       = -MU_EARTH * r / np.linalg.norm(r)**3 + j2_acceleration(r, MU_EARTH)
    Phi_dot = dynamics_jacobian(r) @ Phi

    dydt      = np.zeros(42)
    dydt[0:3] = v
    dydt[3:6] = a
    dydt[6:]  = Phi_dot.flatten()
    return dydt


# ── Propagator ────────────────────────────────────────────────────────────────

class Propagator:
    """
    Numerical orbit propagator using RK45 with J2-perturbed two-body dynamics.

    Maneuvers are handled via an optional control_law passed to propagate().
    When control_law is None the standard unperturbed dynamics are used.
    When control_law is provided the thrusting dynamics are used instead,
    with the control law called at every integrator function evaluation.

    Args:
        rtol : relative tolerance for RK45 (default 1e-10)
        atol : absolute tolerance for RK45 (default 1e-10)

    Examples
    --------
    # Coast (no thrust)
    states = propagator.propagate(state0, 3600.0)

    # Burn using a control law
    def my_law(t, y):
        return np.array([0.0, 0.0, 1e-3])   # 1 mm/s² in ECI Z

    states = propagator.propagate(state0, 60.0, control_law=my_law)

    # EKF predict step
    new_state, stm = propagator.propagate_with_stm(state0, dt=10.0)
    P_new = stm @ P @ stm.T + Q
    """

    def __init__(self, rtol: float = 1e-10, atol: float = 1e-10) -> None:
        self.rtol = rtol
        self.atol = atol

    # ── Internal integration helper ───────────────────────────────────────────

    def _integrate(
        self,
        y0: np.ndarray,
        duration: float,
        rhs: Callable,
        eval_times: Optional[np.ndarray] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Thin wrapper around solve_ivp from t=0 to t=duration.

        Returns:
            (t_out, y_out) arrays from the solver.
        """
        sol = solve_ivp(
            rhs,
            (0.0, duration),
            y0,
            method="RK45",
            t_eval=eval_times,
            rtol=self.rtol,
            atol=self.atol,
            dense_output=False,
        )
        if not sol.success:
            raise RuntimeError(f"Integration failed: {sol.message}")
        return sol.t, sol.y

    # ── State-only propagation ────────────────────────────────────────────────

    def propagate(
        self,
        initial_state: OrbitalState,
        duration: float,
        eval_times: Optional[np.ndarray] = None,
        control_law: Optional[ControlLaw] = None,
    ) -> List[OrbitalState]:
        """
        Propagate state forward for `duration` seconds.

        When control_law is None the standard J2 dynamics are used.
        When control_law is provided the thrusting dynamics are used,
        with thrust added at every RK45 function evaluation — so the
        integrator naturally adapts its step size to the thrust profile.

        Args:
            initial_state : starting OrbitalState
            duration      : propagation time in seconds (must be > 0)
            eval_times    : (N,) output times in seconds from epoch.
                            Solver chooses its own steps if None.
            control_law   : optional callable (t, state) → (3,) ECI accel m/s²

        Returns:
            List[OrbitalState] at each output time.
        """
        if duration <= 0:
            raise ValueError(f"duration must be positive, got {duration}")

        if control_law is None:
            rhs = dynamics_tb_j2
        else:
            # Wrap to match scipy's f(t, y) signature
            def rhs(t: float, y: np.ndarray) -> np.ndarray:
                return dynamics_tb_j2_thrust(t, y, control_law)

        t_out, y_out = self._integrate(
            initial_state.state_vector,
            duration,
            rhs,
            eval_times,
        )

        return [
            OrbitalState.from_state_vector(
                y_out[:, k],
                initial_state.epoch + timedelta(seconds=float(t_out[k])),
            )
            for k in range(t_out.size)
        ]

    def propagate_to_epoch(
        self,
        initial_state: OrbitalState,
        target_epoch: datetime,
        control_law: Optional[ControlLaw] = None,
    ) -> OrbitalState:
        """Propagate to a specific UTC epoch and return a single OrbitalState."""
        dt = (target_epoch - initial_state.epoch).total_seconds()
        if dt < 0:
            raise ValueError(
                f"target_epoch {target_epoch} is before "
                f"initial epoch {initial_state.epoch}"
            )
        if dt == 0.0:
            return initial_state.copy()
        return self.propagate(
            initial_state, dt,
            eval_times=np.array([dt]),
            control_law=control_law,
        )[-1]

    def propagate_over_span(
        self,
        initial_state: OrbitalState,
        start_epoch: datetime,
        end_epoch: datetime,
        dt: float,
        control_law: Optional[ControlLaw] = None,
    ) -> List[OrbitalState]:
        """
        Propagate over a time span at a fixed output cadence.

        If initial_state.epoch < start_epoch the state is first brought
        forward to start_epoch (coasting, no thrust).

        Args:
            initial_state : reference OrbitalState (epoch ≤ start_epoch)
            start_epoch   : first output epoch
            end_epoch     : last output epoch
            dt            : output cadence in seconds
            control_law   : optional thrust callable

        Returns:
            List[OrbitalState] from start_epoch to end_epoch.
        """
        total      = (end_epoch - start_epoch).total_seconds()
        eval_times = np.arange(0.0, total + dt * 0.5, dt)

        if initial_state.epoch != start_epoch:
            state_at_start = self.propagate_to_epoch(initial_state, start_epoch)
        else:
            state_at_start = initial_state

        return self.propagate(
            state_at_start, total,
            eval_times=eval_times,
            control_law=control_law,
        )

    # ── State + STM propagation (EKF) ─────────────────────────────────────────

    def propagate_with_stm(
        self,
        initial_state: OrbitalState,
        duration: float,
    ) -> tuple[OrbitalState, np.ndarray]:
        """
        Propagate state and STM together over `duration` seconds (no thrust).

        Called by the EKF predict step. The 42-dimensional augmented ODE
        integrates [r, v, Φ] in a single solve_ivp call with Φ(0) = I₆.

        After propagation the covariance is updated externally:
            P⁺ = Φ · P · Φᵀ + Q

        Args:
            initial_state : OrbitalState at the current EKF epoch
            duration      : propagation interval in seconds

        Returns:
            new_state : OrbitalState at t + duration
            stm       : (6,6) State Transition Matrix Φ(t+dt, t)
        """
        if duration <= 0:
            raise ValueError(f"duration must be positive, got {duration}")

        y0      = np.zeros(42)
        y0[0:6] = initial_state.state_vector
        y0[6:]  = np.eye(6).flatten()

        _, y_out = self._integrate(
            y0,
            duration,
            dynamics_with_stm,
            eval_times=np.array([duration]),
        )

        final     = y_out[:, -1]
        new_state = OrbitalState.from_state_vector(
            final[0:6],
            initial_state.epoch + timedelta(seconds=duration),
        )
        stm = final[6:].reshape((6, 6))

        return new_state, stm