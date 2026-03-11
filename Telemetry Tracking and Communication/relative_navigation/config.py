from datetime import datetime, timezone
import numpy as np
from utils.orbital_frame_conversions import *

# Integration tolerance
tol = 1e-9

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

# Attitude information
TARGET_EP = np.array([0, 0, 0, 1]) # unit quaternion to start
TARGET_ANG_VEL = np.array([0, 0, 0]) # no ang vel to start

# Moment of Intertia
m_r = 3000   # kg  — rocket mass

R   = 1.85   # m   — rocket radius
L   = 30.0   # m   — rocket length

# Inertia of rocket about its own CM (solid cylinder)
TARGET_I = np.diag([
    (1/12) * m_r * (3*R**2 + L**2),   # Ixx
    (1/12) * m_r * (3*R**2 + L**2),   # Iyy
    (1/2)  * m_r * R**2,               # Izz
])



# -- Spacecraft Initial State --------------------------
# Translational State
# Done in main

# Attitude information
CHASER_EP = np.array([0, 0, 0, 1]) # unit quaternion to start
CHASER_ANG_VEL = np.array([0, 0, 0]) # no ang vel to start



# Moment of Intertia
m_s = 1000   # kg  — satellite mass
s   = 1.0    # m   — satellite side length
# Inertia of satellite about its own CM (solid cube)
CHASER_I = np.diag([
    (1/6) * m_s * s**2,   # Ixx
    (1/6) * m_s * s**2,   # Iyy
    (1/6) * m_s * s**2,   # Izz
])

# ── Simulation window ─────────────────────────────────────────────────────────
# Start exactly at the TLE epoch for maximum SGP4 accuracy.
# TLE fractional day 62.27983304:
#   0.27983304 * 86400 = 24,175.5 s  →  06 h 42 m 55.5 s UTC

SIM_START = 0 # datetime(2026, 3, 3, 6, 42, 55, 500000, tzinfo=timezone.utc)
SIM_END   = 24 * 60 * 60  # datetime(2026, 3, 4, 6, 42, 55, 500000, tzinfo=timezone.utc)

# Truth propagation output cadence (seconds)
# 30 s → ~2880 states over 24 h, fine enough for smooth plots
SIM_DT = 30.0