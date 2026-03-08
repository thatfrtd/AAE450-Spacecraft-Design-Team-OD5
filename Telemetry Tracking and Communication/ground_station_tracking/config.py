"""
config.py

Central configuration for the spacecraft tracking simulation.

All tuneable parameters live here. Units are SI (metres, seconds, radians)
unless explicitly noted otherwise. main.py reads this module and should
never have magic numbers of its own.
"""

from datetime import datetime, timezone
import numpy as np

# ── Target TLE ────────────────────────────────────────────────────────────────
# Source : space-track.org
# NORAD  : 37766  (object 11039B)
# Epoch  : 2026 day 62.27983304  →  2026-03-03 06:42:55.5 UTC

TARGET_TLE_LINE1 = (
    "1 37766U 11039B   26062.27983304  .00000976  00000-0  16844-3 0  9993"
)
TARGET_TLE_LINE2 = (
    "2 37766  98.1167 320.5520 0041526 288.1223  71.5463 14.70228878781238"
)
TARGET_NAME = "37766 (11039B)"

# ── Simulation window ─────────────────────────────────────────────────────────
# Start exactly at the TLE epoch for maximum SGP4 accuracy.
# TLE fractional day 62.27983304:
#   0.27983304 * 86400 = 24,175.5 s  →  06 h 42 m 55.5 s UTC

SIM_START = datetime(2026, 3, 3, 6, 42, 55, 500000, tzinfo=timezone.utc)
SIM_END   = datetime(2026, 3, 4, 6, 42, 55, 500000, tzinfo=timezone.utc)

# Truth propagation output cadence (seconds)
# 30 s → ~2880 states over 24 h, fine enough for smooth plots
TRUTH_DT = 30.0

# ── Chaser initial Keplerian elements ─────────────────────────────────────────
# Defined at SIM_START in ECI J2000.
#
# Only semi-major axis is specified here. All other elements (e, i, RAAN,
# argp, nu) are set to None and extracted from the target's propagated initial
# state at runtime in main.py, so the chaser exactly matches the target's
# osculating values rather than the TLE mean elements.
#
# Units: a in metres, angles in radians.

CHASER_NAME = "Chaser"

# Semi-major axis: 600 km altitude → a = RE + 600 km
# RE (WGS-84 equatorial) = 6 378 137.0 m
CHASER_A = (6_378_137.0 + 600_000.0)    # metres   →  6 978 137.0 m

CHASER_E = None                          # matched to target osculating state at runtime

CHASER_I = None                          # matched to target osculating state at runtime

# RAAN, argp, nu, e, i: all None — populated from target osculating state in main.py
CHASER_RAAN = None
CHASER_ARGP = None
CHASER_NU   = None

# ── Propagator tolerances ─────────────────────────────────────────────────────
PROP_RTOL = 1e-10
PROP_ATOL = 1e-10

# ── EKF parameters ────────────────────────────────────────────────────────────

# Initial covariance P0 — 1-sigma position and velocity uncertainty at filter start
EKF_INIT_POS_STD = 1_000.0      # metres   (1 km position uncertainty)
EKF_INIT_VEL_STD = 10.0         # m/s      (10 m/s velocity uncertainty)

EKF_P0 = np.diag([
    EKF_INIT_POS_STD**2,  EKF_INIT_POS_STD**2,  EKF_INIT_POS_STD**2,
    EKF_INIT_VEL_STD**2,  EKF_INIT_VEL_STD**2,  EKF_INIT_VEL_STD**2,
])

# Process noise Q — continuous spectral density
# Represents unmodelled accelerations (drag, SRP, maneuvers, model error)
# Shape (2,2): [position noise spectral density, velocity noise spectral density]
# Units: m^2/s^3
EKF_Q = np.diag([1e-8, 1e-6])

# Measurement noise R — diagonal, one entry per observable
# [range (m^2), range_rate (m^2/s^2), azimuth (rad^2), elevation (rad^2)]
_RANGE_STD       = 50.0               # metres
_RANGE_RATE_STD  = 0.1                # m/s
_AZ_STD          = np.deg2rad(0.1)    # radians
_EL_STD          = np.deg2rad(0.1)    # radians

EKF_R = np.diag([
    _RANGE_STD**2,
    _RANGE_RATE_STD**2,
    _AZ_STD**2,
    _EL_STD**2,
])

# ── Data paths ────────────────────────────────────────────────────────────────
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))

GROUND_STATIONS_JSON = _os.path.join(_HERE, "data", "ground_stations.json")

# ── Output ────────────────────────────────────────────────────────────────────
OUTPUT_DIR  = "outputs"
FIGURE_DPI  = 150