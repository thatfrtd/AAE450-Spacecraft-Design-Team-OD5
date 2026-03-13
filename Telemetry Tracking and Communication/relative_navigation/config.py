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

Target_kep = tle_to_keplerian(TARGET_TLE_LINE1, TARGET_TLE_LINE2)
tar_pos, tar_vel = keplerian_to_cartesian(Target_kep)

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
R_IN = inert_to_RTN_313(Target_kep.raan, Target_kep.i, Target_kep.argp)
chaser_pos_RTN = R_IN.T @ np.array([0, -0.08, 0])
chaser_pos = tar_pos + chaser_pos_RTN
chaser_vel = tar_vel  # + np.array([0, 0.001, 0])  

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

# -- Parameter Initial State ----------------------------
b_w_c = np.zeros(3) # gyro bias in the chaser frame
ep_S_S = np.zeros(3) # star-camra misalignment
ep_O_O = np.zeros(3) # optical camera misalignment

# tau values
tau_b = 0 # TODO: Update
tau_s = 0
tau_o = 0
sigma_b = 0
sigma_s = 0
sigma_o = 0

# ── Simulation window ─────────────────────────────────────────────────────────
# Start exactly at the TLE epoch for maximum SGP4 accuracy.
# TLE fractional day 62.27983304:
#   0.27983304 * 86400 = 24,175.5 s  →  06 h 42 m 55.5 s UTC

SIM_START = 0 # datetime(2026, 3, 3, 6, 42, 55, 500000, tzinfo=timezone.utc)
SIM_END   = 24 * 60 * 60  # datetime(2026, 3, 4, 6, 42, 55, 500000, tzinfo=timezone.utc)

# Truth propagation output cadence (seconds)
# 30 s → ~2880 states over 24 h, fine enough for smooth plots
SIM_DT = 10.0

N = int(SIM_END / SIM_DT)# number of iterations

# ── Model parameters ─────────────────────────────────────────────────────────
# Truth noise values
sigma_vel_T = (0.06 * 1e-6) # km/s, (0.06 mm/s)
sigma_omega_T = (0.001 * 1e-3) # rad/s
sigma_vel_C = (0.06 * 1e-6) # km/s, (0.06 mm/s)
sigma_omega_C = (0.001 * 1e-3) # rad/s

# Initial 1-sigma uncertainties (used to build P0 and scatter x_est)
sigma_r_T, sigma_v_T = 10.0, 0.1     # target pos/vel (m, m/s)
sigma_th_T, sigma_w_T = 0.01, 0.001  # target att (rad), ang vel (rad/s)
sigma_r_C, sigma_v_C = 5.0, 0.05
sigma_th_C = 0.001
sigma_bw, sigma_Ss, sigma_Oo = 1e-4, 1e-3, 1e-3

POS_TOL = 0.1    # m
ATT_TOL = 0.01   # rad


# -- Sensor parameters ---------------------------------------------------------
DT_GYRO         = 10.0    # gyro cadence (s)
f_w             = 3 * 1e-6  # scale factor bias

DT_STAR_TRACKER = 10.0

DT_OPT_CAM      = 60.0
sigma_alpha     = 1e-9 # mrad/axis
sigma_elev      = 1e-9 # mrad/axis

# -- Guidance parameters -------------------------------------------------------
r_dock_T = np.array([0, 0, -4]) * 1e-3 # Dock with the engine 4 m below centroid in Target frame
r_rel_des_D = np.zeros(3) # Desired relative position
r_attach_C = np.array([0, 2, 0]) * 1e-3 # Telescopic arm off of chaser top

v_rel_des_D = np.zeros(3) # no velocity at docking

q_D_T = np.array([0, 0, 0, 1]) # the docking port frame relative to the target frame

w_des_C = np.zeros(3)

