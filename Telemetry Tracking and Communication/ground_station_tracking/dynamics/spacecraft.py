"""
dynamics/spacecraft.py

Spacecraft: owns the truth state of a single vehicle.

Supports two initialization paths:
    - from_tle()          : target spacecraft (non-cooperative, TLE from space-track.org)
    - from_keplerian()    : maneuvering chaser spacecraft
    - from_state_vector() : direct Cartesian ECI initialisation

Maneuvers are handled by attaching a control_law to the chaser spacecraft.
The control law is a plain callable with signature:

    def control_law(t: float, state: np.ndarray) -> np.ndarray:
        '''
        Args:
            t     : elapsed seconds since the start of the current propagation
            state : (6,) [rx, ry, rz, vx, vy, vz] ECI, metres / m·s⁻¹
        Returns:
            (3,) thrust acceleration in ECI frame, m/s²
                 return np.zeros(3) when the thruster is off
        '''

The control law is passed directly into the propagator's ODE, so the
integrator sees the thrust at every function evaluation — no trajectory
splitting, no impulsive approximation.

Frame note:
    SGP4 returns states in TEME. The TEME→J2000 rotation is sub-arcsecond
    and is treated as negligible for this simulation (see coordinates.py).
"""

from __future__ import annotations

import numpy as np
from datetime import datetime
from typing import Callable, List, Optional

from sgp4.api import Satrec, jday

from dynamics.orbital_state import OrbitalState
from dynamics.propagator import Propagator
from utils.coordinates import MU_EARTH

_KM_TO_M = 1_000.0

# Type alias — matches the interface expected by Propagator.propagate()
ControlLaw = Callable[[float, np.ndarray], np.ndarray]


class Spacecraft:
    """
    Truth model for a single spacecraft.

    Attributes:
        name             : human-readable identifier
        initial_state    : OrbitalState at the reference epoch
        propagator       : Propagator instance used for truth integration
        control_law      : optional thrust callable (t, state) → (3,) accel m/s²
                           None for the non-maneuvering target spacecraft
        truth_trajectory : list of OrbitalState populated by
                           generate_truth_trajectory()

    Typical usage:
        # Target (non-cooperative, no thrust)
        target = Spacecraft.from_tle(line1, line2, epoch)
        target.generate_truth_trajectory(t_start, t_end, dt=10.0)

        # Chaser (maneuvering)
        chaser = Spacecraft.from_keplerian(kep0, epoch,
                                           control_law=my_control_law)
        chaser.generate_truth_trajectory(t_start, t_end, dt=10.0)
    """

    def __init__(
        self,
        initial_state: OrbitalState,
        name: str = "spacecraft",
        control_law: Optional[ControlLaw] = None,
        propagator: Optional[Propagator] = None,
    ) -> None:
        self.name            = name
        self.initial_state   = initial_state
        self.control_law     = control_law          # None → coast
        self.propagator      = propagator or Propagator()
        self.truth_trajectory: List[OrbitalState] = []

    # ── Constructors ───────────────────────────────────────────────────────────

    @classmethod
    def from_tle(
        cls,
        tle_line1: str,
        tle_line2: str,
        epoch: datetime,
        name: str = "target",
        control_law: Optional[ControlLaw] = None,
        propagator: Optional[Propagator] = None,
    ) -> Spacecraft:
        """
        Initialise from a TLE pair (space-track.org format).

        SGP4 evaluates the TLE at `epoch`, producing an initial Cartesian
        state in TEME (treated as ECI J2000 here — see module docstring).

        Args:
            tle_line1   : TLE line 1 string
            tle_line2   : TLE line 2 string
            epoch       : UTC datetime at which to evaluate the TLE
            name        : label for this spacecraft
            control_law : optional thrust callable (usually None for target)
            propagator  : optional custom Propagator

        Raises:
            RuntimeError: if SGP4 returns a non-zero error code
        """
        satellite = Satrec.twoline2rv(tle_line1, tle_line2)

        jd, fr = jday(
            epoch.year, epoch.month, epoch.day,
            epoch.hour, epoch.minute,
            epoch.second + epoch.microsecond * 1e-6,
        )
        error, r_teme_km, v_teme_km = satellite.sgp4(jd, fr)

        if error != 0:
            raise RuntimeError(
                f"SGP4 returned error code {error} for '{name}'. "
                "Check that the TLE is valid and the epoch is within its fit span."
            )

        position = np.array(r_teme_km) * _KM_TO_M
        velocity = np.array(v_teme_km) * _KM_TO_M

        initial_state = OrbitalState(position, velocity, epoch)
        return cls(initial_state, name=name,
                   control_law=control_law, propagator=propagator)

    @classmethod
    def from_keplerian(
        cls,
        kep: np.ndarray,
        epoch: datetime,
        mu: float = MU_EARTH,
        name: str = "chaser",
        control_law: Optional[ControlLaw] = None,
        propagator: Optional[Propagator] = None,
    ) -> Spacecraft:
        """
        Initialise from classical Keplerian elements.

        Args:
            kep         : (6,) [a (m), e, i (rad), RAAN (rad), argp (rad), ν (rad)]
            epoch       : UTC epoch for the Keplerian elements
            control_law : optional thrust callable
        """
        initial_state = OrbitalState.from_keplerian(kep, epoch, mu)
        return cls(initial_state, name=name,
                   control_law=control_law, propagator=propagator)

    @classmethod
    def from_state_vector(
        cls,
        state: np.ndarray,
        epoch: datetime,
        name: str = "chaser",
        control_law: Optional[ControlLaw] = None,
        propagator: Optional[Propagator] = None,
    ) -> Spacecraft:
        """
        Initialise directly from a (6,) Cartesian ECI state vector.

        Args:
            state       : (6,) [rx, ry, rz, vx, vy, vz] ECI, metres / m·s⁻¹
            epoch       : UTC epoch for the state
            control_law : optional thrust callable
        """
        initial_state = OrbitalState.from_state_vector(state, epoch)
        return cls(initial_state, name=name,
                   control_law=control_law, propagator=propagator)

    # ── Truth trajectory ───────────────────────────────────────────────────────

    def generate_truth_trajectory(
        self,
        start_epoch: datetime,
        end_epoch: datetime,
        dt: float,
    ) -> List[OrbitalState]:
        """
        Integrate the truth trajectory over [start_epoch, end_epoch].

        If a control_law is attached it is passed into the propagator so
        that thrust is applied continuously throughout the integration.
        The integrator adapts its internal step size to the thrust profile
        automatically — no manual splitting is required.

        Args:
            start_epoch : simulation start (UTC)
            end_epoch   : simulation end   (UTC)
            dt          : output cadence in seconds

        Returns:
            List[OrbitalState] at each output step.
            Also stored in self.truth_trajectory.
        """
        if end_epoch <= start_epoch:
            raise ValueError("end_epoch must be after start_epoch")

        # Bring state to start_epoch if needed
        if self.initial_state.epoch < start_epoch:
            current = self.propagator.propagate_to_epoch(
                self.initial_state, start_epoch
            )
        elif self.initial_state.epoch == start_epoch:
            current = self.initial_state.copy()
        else:
            raise ValueError(
                f"initial_state.epoch ({self.initial_state.epoch}) is after "
                f"start_epoch ({start_epoch})"
            )

        self.truth_trajectory = self.propagator.propagate_over_span(
            current, start_epoch, end_epoch, dt,
            control_law=self.control_law,
        )
        return self.truth_trajectory

    # ── Accessors ──────────────────────────────────────────────────────────────

    def get_state_at(self, epoch: datetime) -> Optional[OrbitalState]:
        """
        Return the truth state with the epoch closest to `epoch`.

        Returns None if the trajectory has not been generated yet.
        For exact-epoch queries use the propagator directly.
        """
        if not self.truth_trajectory:
            return None
        return min(
            self.truth_trajectory,
            key=lambda s: abs((s.epoch - epoch).total_seconds()),
        )

    def set_control_law(self, control_law: Optional[ControlLaw]) -> None:
        """
        Attach or replace the control law.

        Call before generate_truth_trajectory(). Set to None to coast.
        """
        self.control_law = control_law

    # ── Display ────────────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        return (
            f"Spacecraft(name='{self.name}', "
            f"controlled={self.control_law is not None}, "
            f"trajectory_points={len(self.truth_trajectory)}, "
            f"initial={self.initial_state})"
        )