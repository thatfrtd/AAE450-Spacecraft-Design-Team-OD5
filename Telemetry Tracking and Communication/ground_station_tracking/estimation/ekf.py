""" 
hi
estimation/ekf.py

Extended Kalman Filter (EKF) for spacecraft orbit determination.

Formulation
-----------
This is a continuous-discrete EKF:
    - Dynamics are continuous  (J2 two-body ODE)
    - Measurements are discrete (radar observations at each timestep)

Predict step
    The state and covariance are propagated together through the ODE:

        ẋ   = f(x)                      (nonlinear dynamics)
        Ṗ   = F·P + P·Fᵀ + G·Q·Gᵀ      (continuous Riccati equation)
        Φ̇   = F·Φ                        (STM equation)

    where F = ∂f/∂x is the dynamics Jacobian (from dynamics_jacobian),
    Q is the continuous process noise spectral density (m²/s³),
    G maps process noise into the state (noise enters through acceleration).

    Integrating all three together in a single solve_ivp call (144-dim ODE)
    matches the approach in the reference MATLAB implementation exactly.

Update step (Joseph-stabilized form)
    Given measurement z and predicted measurement ẑ = h(x⁻):

        H   = ∂h/∂x   (numerical Jacobian via finite differences)
        W   = H·P⁻·Hᵀ + R      (innovations covariance)
        K   = P⁻·Hᵀ·W⁻¹        (Kalman gain)
        x⁺  = x⁻ + K·(z - ẑ)   (state update)
        P⁺  = (I - K·H)·P⁻·(I - K·H)ᵀ + K·R·Kᵀ   (Joseph form)

    The Joseph form is used (rather than the simpler P⁺ = (I-KH)P⁻) because
    it is numerically stable and guaranteed symmetric positive semi-definite
    even with finite-precision Kalman gain errors.

Measurement function
    h(x) maps ECI state → [range, range_rate, azimuth, elevation] for a
    given ground station. The Jacobian H is computed numerically.

Usage
-----
    ekf = EKF(
        initial_state      = orbital_state,
        initial_covariance = P0,
        process_noise_Q    = Q,
        measurement_noise_R= R,
        propagator         = propagator,
    )

    for state in truth_trajectory:
        obs = station.observe(state)
        ekf.predict(dt)
        if obs is not None:
            ekf.update(obs, station)

        # Access results
        ekf.state      → current OrbitalState estimate
        ekf.covariance → current 6×6 covariance matrix
"""

from __future__ import annotations

import numpy as np
from datetime import timedelta
from typing import List, Optional

from scipy.integrate import solve_ivp

from dynamics.orbital_state import OrbitalState
from dynamics.propagator import j2_acceleration, dynamics_jacobian
from sensors.observation_model import Observation, compute_observation
from utils.coordinates import MU_EARTH


# ── EKF ───────────────────────────────────────────────────────────────────────

class EKF:
    """
    Continuous-discrete Extended Kalman Filter for orbit determination.

    Tracks a single spacecraft in ECI J2000 using ground-based radar
    observations of range, range-rate, azimuth, and elevation.

    Attributes
    ----------
    state      : OrbitalState — current best estimate
    covariance : (6,6) ndarray — current state error covariance
    Q          : (2,2) continuous process noise spectral density
                 diagonal entries are [σ²_pos_accel, σ²_vel_accel] (m²/s³)
    R          : (4,4) measurement noise covariance
                 diagonal entries are [σ²_range (m²), σ²_range_rate (m²/s²),
                                       σ²_az (rad²), σ²_el (rad²)]
    rtol/atol  : integrator tolerances (inherited from Propagator or set here)

    Parameters
    ----------
    initial_state       : OrbitalState at filter start epoch
    initial_covariance  : (6,6) initial P₀
    process_noise_Q     : (2,2) or (6,6) process noise matrix.
                          If (2,2), interpreted as diag [q_pos, q_vel] and
                          expanded via G mapping. If (6,6), used directly.
    measurement_noise_R : (4,4) measurement noise covariance
    rtol                : ODE relative tolerance  (default 1e-10)
    atol                : ODE absolute tolerance  (default 1e-10)
    """

    def __init__(
        self,
        initial_state       : OrbitalState,
        initial_covariance  : np.ndarray,
        process_noise_Q     : np.ndarray,
        measurement_noise_R : np.ndarray,
        rtol                : float = 1e-10,
        atol                : float = 1e-10,
    ) -> None:

        self.state      = initial_state.copy()
        self.covariance = np.asarray(initial_covariance, dtype=float).copy()
        self.R          = np.asarray(measurement_noise_R, dtype=float).copy()
        self.rtol       = rtol
        self.atol       = atol

        # ── Process noise ─────────────────────────────────────────────────────
        # Q can be (2,2) [shorthand] or (6,6) [full].
        # The (2,2) shorthand is expanded via the G matrix:
        #   G = [0₃  |  I₃]ᵀ  — noise enters through the velocity/accel states
        # giving Q_full = G · Q_2x2 · Gᵀ
        Q = np.asarray(process_noise_Q, dtype=float)
        if Q.shape == (2, 2):
            self.Q = self._expand_process_noise(Q)
        elif Q.shape == (6, 6):
            self.Q = Q.copy()
        else:
            raise ValueError(
                f"process_noise_Q must be shape (2,2) or (6,6), got {Q.shape}"
            )

        # ── Logging ───────────────────────────────────────────────────────────
        # Store history for post-processing and plotting
        self.state_history      : List[OrbitalState] = [initial_state.copy()]
        self.covariance_history : List[np.ndarray]   = [self.covariance.copy()]
        self.innovation_history : List[Optional[np.ndarray]] = [None]
        self.innovations_cov_history: List[Optional[np.ndarray]] = [None]

    # ── Process noise expansion ───────────────────────────────────────────────

    @staticmethod
    def _expand_process_noise(Q_2x2: np.ndarray) -> np.ndarray:
        """
        Expand a (2,2) [q_pos, q_vel] noise matrix to full (6,6) via G mapping.

        G = [0₃  I₃]ᵀ  maps [position noise, velocity noise] into the
        6-element state.  Q_full = G · Q_2x2 · Gᵀ

        In practice for orbit determination, position process noise is often
        zero (position is determined by integration of velocity) and only
        velocity/acceleration noise is non-zero.
        """
        G = np.zeros((6, 2))
        G[0:3, 0] = 1.0   # position states receive pos-noise column
        G[3:6, 1] = 1.0   # velocity states receive vel-noise column
        return G @ Q_2x2 @ G.T

    # ── Predict step ──────────────────────────────────────────────────────────

    def predict(self, dt: float) -> None:
        """
        Propagate the state estimate and covariance forward by dt seconds.

        Integrates the augmented ODE:
            ẋ = f(x)
            Ṗ = F·P + P·Fᵀ + Q           (continuous Riccati)
            Φ̇ = F·Φ                        (STM)

        all in a single solve_ivp call (6 + 36 + 36 = 78 states).

        Args:
            dt : propagation interval in seconds (must be > 0)
        """
        if dt <= 0:
            raise ValueError(f"dt must be positive, got {dt}")

        # Pack augmented state: [x(6), P.flatten(36), Phi.flatten(36)]
        y0 = np.zeros(78)
        y0[0:6]    = self.state.state_vector
        y0[6:42]   = self.covariance.flatten()
        y0[42:78]  = np.eye(6).flatten()

        sol = solve_ivp(
            self._predict_ode,
            (0.0, dt),
            y0,
            method   = "RK45",
            t_eval   = np.array([dt]),
            rtol     = self.rtol,
            atol     = self.atol,
        )

        if not sol.success:
            raise RuntimeError(f"EKF predict ODE failed: {sol.message}")

        final = sol.y[:, -1]

        # Unpack
        new_state_vec = final[0:6]
        P_new         = final[6:42].reshape(6, 6)

        # Symmetrise to prevent numerical drift
        P_new = 0.5 * (P_new + P_new.T)

        # Update filter state
        self.state = OrbitalState.from_state_vector(
            new_state_vec,
            self.state.epoch + timedelta(seconds=dt),
        )
        self.covariance = P_new

    def _predict_ode(self, _t: float, y: np.ndarray) -> np.ndarray:
        """
        Augmented ODE RHS for the predict step.

        y = [x(6), P.flatten(36), Phi.flatten(36)]  →  ẏ (78,)
        """
        x   = y[0:6]
        P   = y[6:42].reshape(6, 6)
        Phi = y[42:78].reshape(6, 6)

        r = x[0:3]
        v = x[3:6]

        # State derivative
        a       = -MU_EARTH * r / np.linalg.norm(r)**3 + j2_acceleration(r, MU_EARTH)
        x_dot   = np.concatenate([v, a])

        # Dynamics Jacobian F = ∂f/∂x at current position
        F = dynamics_jacobian(r)

        # Continuous Riccati:  Ṗ = F·P + P·Fᵀ + Q
        P_dot   = F @ P + P @ F.T + self.Q

        # STM equation:  Φ̇ = F·Φ
        Phi_dot = F @ Phi

        dydt         = np.zeros(78)
        dydt[0:6]    = x_dot
        dydt[6:42]   = P_dot.flatten()
        dydt[42:78]  = Phi_dot.flatten()
        return dydt

    # ── Update step ───────────────────────────────────────────────────────────

    def update(
        self,
        observation  : Observation,
        station_ecef : np.ndarray,
        station_lat  : float,
        station_lon  : float,
    ) -> None:
        """
        Update the state estimate with a new radar observation.

        Applies the Joseph-stabilized EKF measurement update.

        Args:
            observation  : Observation dataclass from GroundStation.observe()
            station_ecef : (3,) station ECEF position, metres
            station_lat  : station geodetic latitude, radians
            station_lon  : station longitude, radians
        """
        # ── Predicted measurement ─────────────────────────────────────────────
        z_hat = self._predicted_measurement(
            self.state.state_vector,
            self.state.epoch,
            station_ecef, station_lat, station_lon,
        )

        # ── Measurement Jacobian ──────────────────────────────────────────────
        H = self._measurement_jacobian(
            self.state.state_vector,
            self.state.epoch,
            station_ecef, station_lat, station_lon,
        )

        # ── Innovation ────────────────────────────────────────────────────────
        z   = observation.as_vector()           # [range, range_rate, az, el]
        inn = z - z_hat

        # Wrap azimuth innovation to (-π, π]
        inn[2] = (inn[2] + np.pi) % (2 * np.pi) - np.pi

        # ── Innovations covariance ────────────────────────────────────────────
        W = H @ self.covariance @ H.T + self.R

        # ── Kalman gain ───────────────────────────────────────────────────────
        # K = P⁻ · Hᵀ · W⁻¹   (solve W·Kᵀ = H·P for numerical stability)
        K = np.linalg.solve(W.T, (self.covariance @ H.T).T).T

        # ── State update ──────────────────────────────────────────────────────
        x_new = self.state.state_vector + K @ inn
        self.state = OrbitalState.from_state_vector(x_new, self.state.epoch)

        # ── Joseph-stabilized covariance update ───────────────────────────────
        # P⁺ = (I - KH)·P⁻·(I - KH)ᵀ + K·R·Kᵀ
        I_KH = np.eye(6) - K @ H
        self.covariance = I_KH @ self.covariance @ I_KH.T + K @ self.R @ K.T
        self.covariance = 0.5 * (self.covariance + self.covariance.T)

        # ── Log ───────────────────────────────────────────────────────────────
        self.innovation_history.append(inn)
        self.innovations_cov_history.append(W)

    # ── Measurement function and Jacobian ─────────────────────────────────────

    def _predicted_measurement(
        self,
        state_vec    : np.ndarray,
        epoch,
        station_ecef : np.ndarray,
        station_lat  : float,
        station_lon  : float,
    ) -> np.ndarray:
        """
        Evaluate h(x): map state vector → [range, range_rate, az, el].
        """
        obs = compute_observation(
            r_eci          = state_vec[0:3],
            v_eci          = state_vec[3:6],
            epoch          = epoch,
            r_station_ecef = station_ecef,
            lat_rad        = station_lat,
            lon_rad        = station_lon,
            station_name   = "",
        )
        return obs.as_vector()

    def _measurement_jacobian(
        self,
        state_vec    : np.ndarray,
        epoch,
        station_ecef : np.ndarray,
        station_lat  : float,
        station_lon  : float,
        pert_pos     : float = 1.0,
        pert_vel     : float = 1e-3,
    ) -> np.ndarray:
        """
        Compute the 4×6 measurement Jacobian H = ∂h/∂x via central finite
        differences.

        Args:
            state_vec : (6,) current state estimate
            epoch     : current epoch (datetime)
            pert_pos  : position perturbation size, metres   (default 1 m)
            pert_vel  : velocity perturbation size, m/s      (default 1e-3)

        Returns:
            (4,6) Jacobian matrix H
        """
        pert  = np.array([pert_pos] * 3 + [pert_vel] * 3)
        H     = np.zeros((4, 6))

        for j in range(6):
            dx          = np.zeros(6)
            dx[j]       = pert[j]

            z_plus  = self._predicted_measurement(
                state_vec + dx, epoch,
                station_ecef, station_lat, station_lon,
            )
            z_minus = self._predicted_measurement(
                state_vec - dx, epoch,
                station_ecef, station_lat, station_lon,
            )

            diff      = z_plus - z_minus
            # Wrap azimuth difference
            diff[2]   = (diff[2] + np.pi) % (2 * np.pi) - np.pi

            H[:, j]   = diff / (2.0 * pert[j])

        return H

    # ── Logging ───────────────────────────────────────────────────────────────

    def log(self) -> None:
        """Append current state and covariance to history."""
        self.state_history.append(self.state.copy())
        self.covariance_history.append(self.covariance.copy())

    # ── Convenience ───────────────────────────────────────────────────────────

    @property
    def position(self) -> np.ndarray:
        """Current estimated ECI position, metres."""
        return self.state.position.copy()

    @property
    def velocity(self) -> np.ndarray:
        """Current estimated ECI velocity, m/s."""
        return self.state.velocity.copy()

    @property
    def sigma(self) -> np.ndarray:
        """Current 1-σ standard deviations (6,) from diagonal of covariance."""
        return np.sqrt(np.maximum(np.diag(self.covariance), 0.0))

    def __repr__(self) -> str:
        sig = self.sigma
        return (
            f"EKF(epoch={self.state.epoch.isoformat()}, "
            f"σ_pos=[{sig[0]:.2f}, {sig[1]:.2f}, {sig[2]:.2f}] m, "
            f"σ_vel=[{sig[3]:.4f}, {sig[4]:.4f}, {sig[5]:.4f}] m/s)"
        )