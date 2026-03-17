from datetime import datetime, timezone
import numpy as np
from utils.orbital_frame_conversions import *

# Integration tolerance
tol = 1e-6

# ── Target TLE ────────────────────────────────────────────────────────────────
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
TARGET_EP      = np.array([0, 0, 0, 1], dtype=float)
TARGET_ANG_VEL = np.array([0, 0, 0],    dtype=float)

# Moment of Inertia — solid cylinder
m_r = 3000.0   # kg
R   = 1.85     # m  radius
L   = 30.0     # m  length
TARGET_I = np.diag([
    (1/12) * m_r * (3*R**2 + L**2),
    (1/12) * m_r * (3*R**2 + L**2),
    (1/2)  * m_r * R**2,
])

# ── Chaser Initial State ──────────────────────────────────────────────────────
R_IN          = inert_to_RTN_313(Target_kep.raan, Target_kep.i, Target_kep.argp)
chaser_pos_RTN = R_IN.T @ np.array([0, -0.08, 0])
chaser_pos    = tar_pos + chaser_pos_RTN
chaser_vel    = tar_vel

CHASER_EP      = np.array([0, 0, 0, 1], dtype=float)
CHASER_ANG_VEL = np.array([0, 0, 0],    dtype=float)

# Moment of Inertia — solid cube
m_s = 1000.0   # kg
s   = 1.0      # m  side
CHASER_I = np.diag([
    (1/6) * m_s * s**2,
    (1/6) * m_s * s**2,
    (1/6) * m_s * s**2,
])

# ── Parameter Initial State ───────────────────────────────────────────────────
b_w_c  = np.zeros(3)   # gyro bias
ep_S_S = np.zeros(3)   # star-camera misalignment
ep_O_O = np.zeros(3)   # optical camera misalignment

# First-order Markov time constants and noise intensities
# Set tau=1e9, sigma=0 to freeze parameters (no drift) during initial testing
tau_b   = 1e9    # s  — gyro bias correlation time
tau_s   = 1e9    # s  — star tracker misalignment
tau_o   = 1e9    # s  — optical camera misalignment
sigma_b = 0.0    # rad/s — gyro bias random walk (frozen)
sigma_s = 0.0    # rad   — star tracker misalignment random walk
sigma_o = 0.0    # rad   — optical camera misalignment random walk

# ── Simulation Window ─────────────────────────────────────────────────────────
SIM_START = 0
SIM_END   = 0.5 * 60 * 60      

SIM_DT = 10.0              # s — timestep
N      = int(SIM_END / SIM_DT)

# ── Truth Process Noise ───────────────────────────────────────────────────────
# These feed into the noise vector passed to dynamics_truth.
# Named sigma_vel_* to avoid collision with EKF uncertainty sigmas below.
# These values match the paper. 
sigma_vel_T   = 0.06e-6 * np.ones(3)   # km/s  (0.06 mm/s per axis)
sigma_omega_T = 1e-6    * np.ones(3)   # rad/s
sigma_vel_C   = 0.06e-6 * np.ones(3)   # km/s
sigma_omega_C = 1e-6    * np.ones(3)   # rad/s

# ── EKF Initial 1-sigma Uncertainties (build P0 and scatter x_est) ───────────
# Scalars — each is applied to a 3-element block via * np.ones(3) in the sim
sigma_r_T  = 0.010    # km  (10 m)
sigma_v_T  = 1e-4     # km/s (0.1 m/s)
sigma_th_T = 0.01     # rad
sigma_w_T  = 0.001    # rad/s

sigma_r_C  = 0.005    # km  (5 m)
sigma_v_C  = 5e-5     # km/s (0.05 m/s)
sigma_th_C = 0.001    # rad

sigma_bw   = 1e-4
sigma_Ss   = 1e-3
sigma_Oo   = 1e-3

# ── Tolerances ────────────────────────────────────────────────────────────────
# These are to break the simulation when it has been effectively captured. 
POS_TOL = 0.0001    # km
ATT_TOL = 0.01   # rad

# ── Sensor Parameters ─────────────────────────────────────────────────────────
DT_GYRO         = 10.0    # s — gyro cadence
f_w             = 3e-6    # scale factor bias

DT_STAR_TRACKER = 10.0    # s
sigma_st        = np.deg2rad(0.01)          # rad  — 0.01° 1-sigma per axis
R_STAR_TRACKER  = (sigma_st**2) * np.eye(3) # (3×3) attitude noise covariance

DT_OPT_CAM      = 60.0    # s
sigma_alpha     = np.deg2rad(0.1 * 1e-3)    # rad  — 0.1 mrad azimuth noise
sigma_elev      = np.deg2rad(0.1 * 1e-3)    # rad  — 0.1 mrad elevation noise
R_OPT_CAM       = np.diag([sigma_alpha**2, sigma_elev**2])  # (2×2)

# ── Guidance Parameters ───────────────────────────────────────────────────────
# This provides information on where the docking port is and what the desired
# final state of the chaser is for sucessful rendezvous. 
r_dock_T    = np.array([0, 0, -4]) * 1e-3   # km — docking port in target frame
r_rel_des_D = np.zeros(3)
r_attach_C  = np.array([0, 2, 0])  * 1e-3   # km — attach point in chaser frame

v_rel_des_D = np.zeros(3)
q_D_T       = np.array([0, 0, 0, 1], dtype=float)
w_des_C     = np.zeros(3)