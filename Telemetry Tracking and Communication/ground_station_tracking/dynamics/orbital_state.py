"""
dynamics/orbital_state.py

OrbitalState: the core data container for a spacecraft state at a single epoch.

Stores position and velocity in ECI J2000 (metres, m/s) plus an optional
6x6 covariance matrix. Provides conversions to Keplerian elements and ECEF.
"""

from __future__ import annotations

import numpy as np
from datetime import datetime          # import the CLASS not the module
from typing import Optional

from utils.coordinates import (
    cartesian_to_keplerian,
    keplerian_to_cartesian,
    eci_to_ecef,
    MU_EARTH,
)


class OrbitalState:
    """
    Spacecraft state at a single epoch in ECI J2000.

    Attributes:
        position   : (3,) ECI position, metres
        velocity   : (3,) ECI velocity, m/s
        epoch      : UTC datetime of the state
        covariance : (6,6) state covariance matrix (zeros if not set)

    All angles produced by conversion methods are in radians.
    """

    def __init__(
        self,
        position: np.ndarray,
        velocity: np.ndarray,
        epoch: datetime,
        covariance: Optional[np.ndarray] = None,
    ) -> None:
        self.position:   np.ndarray = np.asarray(position, dtype=float)
        self.velocity:   np.ndarray = np.asarray(velocity, dtype=float)
        self.epoch:      datetime   = epoch
        self.covariance: np.ndarray = (
            np.asarray(covariance, dtype=float)
            if covariance is not None
            else np.zeros((6, 6))
        )
        self._validate()              # FIX 1: was self.validate()

    # ── Validation ─────────────────────────────────────────────────────────────

    def _validate(self) -> None:
        if self.position.shape != (3,):
            raise ValueError(
                f"position must be shape (3,), got {self.position.shape}"
            )
        if self.velocity.shape != (3,):
            raise ValueError(
                f"velocity must be shape (3,), got {self.velocity.shape}"
            )
        if self.covariance.shape != (6, 6):
            raise ValueError(
                f"covariance must be shape (6,6), got {self.covariance.shape}"
            )

    # ── Properties ─────────────────────────────────────────────────────────────

    @property
    def state_vector(self) -> np.ndarray:
        """(6,) Cartesian state [rx, ry, rz, vx, vy, vz] in ECI."""
        return np.concatenate([self.position, self.velocity])  # FIX 2: added return

    @property
    def radius(self) -> float:
        """Distance from Earth centre in metres."""
        return float(np.linalg.norm(self.position))

    @property
    def speed(self) -> float:
        """Scalar speed in m/s."""
        return float(np.linalg.norm(self.velocity))

    # ── Constructors ───────────────────────────────────────────────────────────

    @classmethod
    def from_state_vector(
        cls,
        state: np.ndarray,
        epoch: datetime,
        covariance: Optional[np.ndarray] = None,
    ) -> OrbitalState:
        """Construct from a (6,) [r, v] array."""
        state = np.asarray(state, dtype=float)
        if state.shape != (6,):
            raise ValueError(f"state must be shape (6,), got {state.shape}")
        return cls(state[:3], state[3:], epoch, covariance)

    @classmethod
    def from_keplerian(
        cls,
        kep: np.ndarray,
        epoch: datetime,
        mu: float = MU_EARTH,
        covariance: Optional[np.ndarray] = None,
    ) -> OrbitalState:
        """
        Construct from classical Keplerian elements.

        Args:
            kep : (6,) [a (m), e, i (rad), RAAN (rad), argp (rad), nu (rad)]
        """
        state = keplerian_to_cartesian(np.asarray(kep, dtype=float), mu)
        return cls(state[:3], state[3:], epoch, covariance)

    # ── Conversions ────────────────────────────────────────────────────────────

    def to_keplerian(self, mu: float = MU_EARTH) -> np.ndarray:
        """
        Return classical Keplerian elements.

        Returns:
            (6,) [a (m), e, i (rad), RAAN (rad), argp (rad), nu (rad)]
        """
        return cartesian_to_keplerian(self.state_vector, mu)

    def to_ecef_position(self) -> np.ndarray:
        """Return (3,) position vector in ECEF, metres."""
        return eci_to_ecef(self.position, self.epoch)

    # ── Utilities ──────────────────────────────────────────────────────────────

    def copy(self) -> OrbitalState:
        """Return a deep copy of this state."""
        return OrbitalState(
            self.position.copy(),
            self.velocity.copy(),
            self.epoch,
            self.covariance.copy(),
        )

    def with_covariance(self, covariance: np.ndarray) -> OrbitalState:
        """Return a new OrbitalState with the given covariance."""
        return OrbitalState(
            self.position, self.velocity, self.epoch, covariance
        )

    def time_since(self, other: OrbitalState) -> float:
        """Seconds elapsed since another OrbitalState's epoch."""
        return (self.epoch - other.epoch).total_seconds()

    # ── Display ────────────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        kep = self.to_keplerian()
        return (
            f"OrbitalState("
            f"epoch={self.epoch.isoformat()}, "
            f"a={kep[0]/1e3:.3f} km, "
            f"e={kep[1]:.6f}, "
            f"i={np.rad2deg(kep[2]):.4f} deg, "
            f"|r|={self.radius/1e3:.3f} km, "
            f"|v|={self.speed:.4f} m/s)"
        )