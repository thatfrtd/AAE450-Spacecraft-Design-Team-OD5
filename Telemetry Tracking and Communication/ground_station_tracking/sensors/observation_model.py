"""
sensors/observation_model.py

Defines the Observation dataclass and the core geometry functions that map
a spacecraft ECI state + ground station position to radar observables.

Observables
-----------
    range       : scalar distance from station to spacecraft, metres
    range_rate  : radial velocity (positive = receding), m/s
    azimuth     : clockwise from North in the local NEU horizon frame, radians
    elevation   : angle above the local horizon, radians

Frame convention (NEU)
----------------------
    The local horizon frame at the station is defined as:
        N  — points geographic North along the ellipsoid surface
        E  — points geographic East  along the ellipsoid surface
        U  — points Up (normal to ellipsoid surface, away from Earth centre)

    Azimuth is measured clockwise from N toward E  (standard radar convention).
    Elevation is measured upward from the N-E plane.

    Visibility condition:  elevation >= station.el_mask_rad
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from datetime import datetime
from typing import Optional

from utils.coordinates import eci_to_ecef, eci_vel_to_ecef


# ── Observation dataclass ─────────────────────────────────────────────────────

@dataclass
class Observation:
    """
    A single radar observation of one spacecraft by one ground station.

    All angles in radians, distances in metres, rates in m/s.

    Attributes:
        station_name : name of the observing ground station
        epoch        : UTC time of the observation
        range        : slant range, metres
        range_rate   : radial velocity (positive = receding), m/s
        azimuth      : clockwise from North, radians  [0, 2π)
        elevation    : above local horizon, radians   (-π/2, π/2]
    """
    station_name : str
    epoch        : datetime
    range        : float
    range_rate   : float
    azimuth      : float
    elevation    : float

    def as_vector(self) -> np.ndarray:
        """Return (4,) observation vector [range, range_rate, az, el]."""
        return np.array([self.range, self.range_rate,
                         self.azimuth, self.elevation])

    def __repr__(self) -> str:
        return (
            f"Observation(station='{self.station_name}', "
            f"epoch={self.epoch.isoformat()}, "
            f"range={self.range/1e3:.3f} km, "
            f"range_rate={self.range_rate:.4f} m/s, "
            f"az={np.rad2deg(self.azimuth):.2f} deg, "
            f"el={np.rad2deg(self.elevation):.2f} deg)"
        )


# ── NEU rotation matrix ───────────────────────────────────────────────────────

def _neu_rotation(lat_rad: float, lon_rad: float) -> np.ndarray:
    """
    Build the 3×3 rotation matrix R_NEU that rotates an ECEF vector into
    the local North-East-Up (NEU) frame at (lat, lon).

    Row 0 → N component
    Row 1 → E component
    Row 2 → U component

    Args:
        lat_rad : geodetic latitude of the station, radians
        lon_rad : longitude of the station, radians

    Returns:
        (3,3) rotation matrix
    """
    sin_lat, cos_lat = np.sin(lat_rad), np.cos(lat_rad)
    sin_lon, cos_lon = np.sin(lon_rad), np.cos(lon_rad)

    # Unit vectors of the NEU frame expressed in ECEF
    N = np.array([-sin_lat * cos_lon,
                  -sin_lat * sin_lon,
                   cos_lat])
    E = np.array([-sin_lon,
                   cos_lon,
                   0.0])
    U = np.array([ cos_lat * cos_lon,
                   cos_lat * sin_lon,
                   sin_lat])

    return np.vstack([N, E, U])   # shape (3, 3)


# ── Core observation geometry ─────────────────────────────────────────────────

def compute_observation(
    r_eci    : np.ndarray,
    v_eci    : np.ndarray,
    epoch    : datetime,
    r_station_ecef : np.ndarray,
    lat_rad  : float,
    lon_rad  : float,
    station_name   : str,
) -> Observation:
    """
    Compute noiseless radar observables for a spacecraft seen from a station.

    This function assumes the spacecraft IS visible (elevation check is done
    in GroundStation.observe). It does NOT check the elevation mask.

    Args:
        r_eci          : (3,) spacecraft ECI position, metres
        v_eci          : (3,) spacecraft ECI velocity, m/s
        epoch          : UTC epoch of the observation
        r_station_ecef : (3,) station ECEF position, metres
        lat_rad        : station geodetic latitude, radians
        lon_rad        : station longitude, radians
        station_name   : label for the returned Observation

    Returns:
        Observation dataclass with range, range_rate, azimuth, elevation.
    """
    # ── Positions in ECEF ─────────────────────────────────────────────────────
    r_sc_ecef = eci_to_ecef(r_eci, epoch)
    v_sc_ecef = eci_vel_to_ecef(r_eci, v_eci, epoch)

    # Relative position and velocity (station is fixed in ECEF → v_station = 0)
    dr_ecef = r_sc_ecef - r_station_ecef
    dv_ecef = v_sc_ecef                       # station velocity is zero in ECEF

    # ── Range and range-rate ──────────────────────────────────────────────────
    rng      = float(np.linalg.norm(dr_ecef))
    rng_rate = float(np.dot(dr_ecef, dv_ecef) / rng)

    # ── NEU components ────────────────────────────────────────────────────────
    R_neu    = _neu_rotation(lat_rad, lon_rad)
    dr_neu   = R_neu @ dr_ecef               # (N, E, U) components

    n_comp, e_comp, u_comp = dr_neu

    # ── Elevation ─────────────────────────────────────────────────────────────
    # el = arcsin(U / |r|)
    elevation = float(np.arcsin(np.clip(u_comp / rng, -1.0, 1.0)))

    # ── Azimuth ───────────────────────────────────────────────────────────────
    # az = arctan2(E, N), wrapped to [0, 2π)
    azimuth = float(np.arctan2(e_comp, n_comp) % (2.0 * np.pi))

    return Observation(
        station_name = station_name,
        epoch        = epoch,
        range        = rng,
        range_rate   = rng_rate,
        azimuth      = azimuth,
        elevation    = elevation,
    )