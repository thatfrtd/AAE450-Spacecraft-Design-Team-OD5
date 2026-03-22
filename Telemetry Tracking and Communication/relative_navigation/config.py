from datetime import datetime, timezone
from dataclasses import replace
import numpy as np
from utils.orbital_frame_conversions import *

# ── Integration tolerance ─────────────────────────────────────────────────────
tol = 1e-6

# ── Target TLE ────────────────────────────────────────────────────────────────
TARGET_TLE_LINE1 = (
    "1 37766U 11039B   26062.27983304  .00000976  00000-0  16844-3 0  9993"
)
TARGET_TLE_LINE2 = (
    "2 37766  98.1167 320.5520 0041526 288.1223  71.5463 14.70228878781238"
)
TARGET_NAME = "37766 (11039B)"

Target_kep       = tle_to_keplerian(TARGET_TLE_LINE1, TARGET_TLE_LINE2)
tar_pos, tar_vel = keplerian_to_cartesian(Target_kep)

TARGET_EP      = np.array([0, 0, 0, 1], dtype=float)
TARGET_ANG_VEL = np.array([0, 0, 0],    dtype=float)

# Moment of Inertia — solid cylinder [kg·km²]
m_r = 3000.0;  R = 1.85e-3;  L = 30.0e-3
TARGET_I = np.diag([
    (1/12) * m_r * (3*R**2 + L**2),
    (1/12) * m_r * (3*R**2 + L**2),
    (1/2)  * m_r * R**2,
])

# ── Chaser Initial State ──────────────────────────────────────────────────────
# R_IN           = inert_to_RTN_313(Target_kep.raan, Target_kep.i, Target_kep.argp)
# chaser_pos_RTN = R_IN.T @ np.array([0, 0, -0.3])   # 300 m behind in RTN [km]
# chaser_pos     = tar_pos + chaser_pos_RTN
# chaser_vel     = tar_vel
# How far behind the chaser starts (positive = behind target in orbit)
delta_nu = np.deg2rad(0.0003)   # ~300 m behind at 6500 km altitude (arc ≈ r*Δν)

# Create chaser elements with reduced true anomaly
# Also recompute mean anomaly M from the new nu using the same eccentricity
e   = Target_kep.e
nu  = Target_kep.nu - delta_nu

# Convert true anomaly → eccentric anomaly → mean anomaly
E   = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2),
                      np.sqrt(1 + e) * np.cos(nu / 2))
M   = E - e * np.sin(E)

chaser_kep = replace(Target_kep, nu=nu, M=M)

chaser_pos, chaser_vel = keplerian_to_cartesian(chaser_kep)
chaser_theta = np.deg2rad(180)
CHASER_EP      = np.array([0, 0, 0, 1], dtype=float)
CHASER_ANG_VEL = np.array([0, 0, 0],    dtype=float)

# Moment of Inertia — solid cube [kg·km²]
m_s = 1000.0;  s = 1.0e-3
CHASER_I = np.diag([
    (1/6) * m_s * s**2,
    (1/6) * m_s * s**2,
    (1/6) * m_s * s**2,
])

u_cmd_max = 10

# ── Parameter Initial State ───────────────────────────────────────────────────
b_w_c  = np.zeros(3)
ep_S_S = np.zeros(3)
ep_O_O = np.zeros(3)

tau_b   = 1e9;  tau_s   = 1e9;  tau_o   = 1e9
sigma_b = 1e-8; sigma_s = 1e-8; sigma_o = 1e-8

# ── Simulation Window ─────────────────────────────────────────────────────────
SIM_START = 0
SIM_END   = 0.5 * 60 * 60    # s
SIM_DT    = 1.0               # s
N         = int(SIM_END / SIM_DT)

# ── Truth Process Noise ───────────────────────────────────────────────────────
# Set to zero for no-noise testing — these only add disturbances to dynamics,
# they do not affect any matrix inversions so zero is safe here.
sigma_vel_T   = 0.0 * np.ones(3)
sigma_omega_T = 0.0 * np.ones(3)
sigma_vel_C   = 0.0 * np.ones(3)
sigma_omega_C = 0.0 * np.ones(3)

# ── EKF Initial 1-sigma Uncertainties ────────────────────────────────────────
# IMPORTANT: must be nonzero — P0 = diag(sigma^2) must be positive definite.
# Use 1e-10 for "no noise" testing (effectively zero but invertible).
# Restore realistic values (commented below) for full simulation.
sigma_r_T  = 1e-10   # km     # realistic: 0.010
sigma_v_T  = 1e-10   # km/s   # realistic: 1e-4
sigma_th_T = 1e-10   # rad    # realistic: 0.01
sigma_w_T  = 1e-10   # rad/s  # realistic: 0.001

sigma_r_C  = 1e-10   # km     # realistic: 0.005
sigma_v_C  = 1e-10   # km/s   # realistic: 5e-5
sigma_th_C = 1e-10   # rad    # realistic: 0.001

sigma_bw   = 1e-10           # realistic: 1e-4
sigma_Ss   = 1e-10           # realistic: 1e-3
sigma_Oo   = 1e-10           # realistic: 1e-3

# ── Rendezvous Tolerances ─────────────────────────────────────────────────────
POS_TOL = 0.0001    # km  (10 m) — relaxed for initial testing
ATT_TOL = 2     # rad

# ── Sensor Parameters ─────────────────────────────────────────────────────────
# IMPORTANT: R matrices must be nonzero — S = HPH^T + R must be invertible.
# Use 1e-10 for "no noise" testing. Restore realistic values for full sim.

DT_GYRO  = 1.0    # s
f_w      = 3e-6

DT_STAR_TRACKER = 10.0    # s
sigma_st        = 1e-10                      # realistic: np.deg2rad(0.01)
R_STAR_TRACKER  = (sigma_st**2) * np.eye(3)

DT_OPT_CAM  = 60.0    # s
sigma_alpha = 1e-10                          # realistic: np.deg2rad(0.1e-3)
sigma_elev  = 1e-10                          # realistic: np.deg2rad(0.1e-3)
R_OPT_CAM   = np.diag([sigma_alpha**2, sigma_elev**2])

DT_GPS    = 1.0       # s
sigma_gps = 1e-10     # km    # realistic: 10e-3
R_GPS     = (sigma_gps**2) * np.eye(3)

DT_LIDAR    = 1.0     # s
sigma_lidar = 1e-10   # km    # realistic: 1e-3
R_LIDAR     = np.array([[sigma_lidar**2]])

# Target attitude from feature tracking (proxy for pose estimation)
DT_TARGET_ATT       = 10.0                        # s — same cadence as star tracker
sigma_target_att    = np.deg2rad(0.1)             # rad — 0.1° (worse than star tracker)
R_TARGET_ATT        = (sigma_target_att**2) * np.eye(3)

# ── Guidance Parameters ───────────────────────────────────────────────────────
r_dock_T    = np.array([0,  0, -4]) * 1e-3   # km
r_rel_des_D = np.zeros(3)
r_attach_C  = np.array([0,  2,  0]) * 1e-3   # km

v_rel_des_D = np.zeros(3)
theta_dock = np.deg2rad(-90)
q_D_T       = np.array([0, 0, np.sin(theta_dock), np.cos(theta_dock)], dtype=float)
w_des_C     = np.zeros(3)