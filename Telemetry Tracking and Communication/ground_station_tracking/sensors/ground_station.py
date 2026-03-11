"""
sensors/ground_station.py

GroundStation: models a single radar ground station.

Each station:
    - Has a fixed ECEF position derived from its geodetic coordinates
    - Has an elevation mask angle (per-station, from the JSON)
    - Can check whether a spacecraft is currently visible
    - Returns a noiseless Observation if the spacecraft is in view,
      or None if it is below the elevation mask

Loading
-------
    Use load_stations_from_json() to instantiate all stations from the
    data/ground_stations.json file. The JSON format is:

        [
          {"name": "Midland", "lat": 31.9, "lon": -102.2,
           "alt": 860, "el_mask": 20},
          ...
        ]

    Fields:
        name    : station identifier string
        lat     : geodetic latitude,  degrees
        lon     : longitude,          degrees
        alt     : altitude above WGS-84 ellipsoid, metres
        el_mask : elevation mask angle, degrees
                  spacecraft must be above this angle to be considered visible
"""

from __future__ import annotations

import json
import numpy as np
from datetime import datetime
from typing import List, Optional

from dynamics.orbital_state import OrbitalState
from sensors.observation_model import Observation, compute_observation
from utils.coordinates import geodetic_to_ecef


class GroundStation:
    """
    A single radar ground station.

    Attributes:
        name         : human-readable station identifier
        lat_rad      : geodetic latitude, radians
        lon_rad      : longitude, radians
        alt_m        : altitude above WGS-84 ellipsoid, metres
        el_mask_rad  : elevation mask angle, radians
                       spacecraft must be above this angle to be visible
        ecef_position: (3,) fixed ECEF position, metres
                       pre-computed at construction — does not change

    The station position is fixed in ECEF (Earth rotates with it).
    The ECI position of the station changes over time and is computed
    on-the-fly inside observe() via the ECI↔ECEF rotation.
    """

    def __init__(
        self,
        name      : str,
        lat_deg   : float,
        lon_deg   : float,
        alt_m     : float,
        el_mask_deg: float,
    ) -> None:
        self.name        = name
        self.lat_rad     = np.deg2rad(lat_deg)
        self.lon_rad     = np.deg2rad(lon_deg)
        self.alt_m       = float(alt_m)
        self.el_mask_rad = np.deg2rad(el_mask_deg)

        # Pre-compute fixed ECEF position (never changes)
        self.ecef_position: np.ndarray = geodetic_to_ecef(
            self.lat_rad, self.lon_rad, self.alt_m
        )

    # ── Visibility ────────────────────────────────────────────────────────────

    def is_visible(self, state: OrbitalState) -> bool:
        """
        Check whether a spacecraft is above the elevation mask.

        Computes the spacecraft's elevation angle as seen from this station
        at the state's epoch and compares it to el_mask_rad.

        Args:
            state : OrbitalState of the spacecraft to check

        Returns:
            True if elevation >= el_mask_rad, False otherwise.
        """
        obs = compute_observation(
            r_eci          = state.position,
            v_eci          = state.velocity,
            epoch          = state.epoch,
            r_station_ecef = self.ecef_position,
            lat_rad        = self.lat_rad,
            lon_rad        = self.lon_rad,
            station_name   = self.name,
        )
        return obs.elevation >= self.el_mask_rad

    # ── Observation ───────────────────────────────────────────────────────────

    def observe(self, state: OrbitalState) -> Optional[Observation]:
        """
        Attempt to observe a spacecraft at the given state.

        If the spacecraft is above the elevation mask, returns a noiseless
        Observation. Returns None if the spacecraft is not in view.

        Args:
            state : OrbitalState of the spacecraft to observe

        Returns:
            Observation if visible, None otherwise.
        """
        obs = compute_observation(
            r_eci          = state.position,
            v_eci          = state.velocity,
            epoch          = state.epoch,
            r_station_ecef = self.ecef_position,
            lat_rad        = self.lat_rad,
            lon_rad        = self.lon_rad,
            station_name   = self.name,
        )

        if obs.elevation >= self.el_mask_rad:
            return obs
        return None

    def observe_all(
        self,
        trajectory: List[OrbitalState],
    ) -> List[Optional[Observation]]:
        """
        Observe a spacecraft over an entire trajectory.

        Args:
            trajectory : list of OrbitalState (e.g. truth_trajectory)

        Returns:
            List of the same length as trajectory.
            Each entry is an Observation if visible at that epoch, else None.
        """
        return [self.observe(state) for state in trajectory]

    # ── Display ───────────────────────────────────────────────────────────────

    def __repr__(self) -> str:
        return (
            f"GroundStation(name='{self.name}', "
            f"lat={np.rad2deg(self.lat_rad):.2f} deg, "
            f"lon={np.rad2deg(self.lon_rad):.2f} deg, "
            f"alt={self.alt_m:.0f} m, "
            f"el_mask={np.rad2deg(self.el_mask_rad):.1f} deg)"
        )


# ── JSON loader ───────────────────────────────────────────────────────────────

def load_stations_from_json(filepath: str) -> List[GroundStation]:
    """
    Load a list of GroundStation objects from a JSON file.

    Expected JSON format:
        [
          {
            "name"    : "Midland",
            "lat"     : 31.9,       -- geodetic latitude,  degrees
            "lon"     : -102.2,     -- longitude,          degrees
            "alt"     : 860,        -- altitude,           metres
            "el_mask" : 20          -- elevation mask,     degrees
          },
          ...
        ]

    Args:
        filepath : path to the JSON file

    Returns:
        List[GroundStation], one per entry in the JSON.

    Raises:
        FileNotFoundError : if the file does not exist
        KeyError          : if a required field is missing from an entry
    """
    with open(filepath, "r") as f:
        entries = json.load(f)

    stations = []
    for entry in entries:
        stations.append(GroundStation(
            name       = entry["name"],
            lat_deg    = entry["lat"],
            lon_deg    = entry["lon"],
            alt_m      = entry["alt"],
            el_mask_deg= entry["el_mask"],
        ))

    print(f"  Loaded {len(stations)} ground stations from '{filepath}'")
    return stations