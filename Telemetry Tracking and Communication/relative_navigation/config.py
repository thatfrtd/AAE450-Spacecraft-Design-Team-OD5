from datetime import datetime, timezone
from dataclasses import replace
import numpy as np
from utils.orbital_frame_conversions import *

# ── Integration tolerance ─────────────────────────────────────────────────────
tol = 1e-12

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
# Cylinder axis = target body x-axis
# Nozzle end at x = -L/2 = -4 m = -0.004 km
m_r = 4000.0;  R = 1.85e-3;  L = 30.0e-3   # 8 m length
TARGET_I = np.diag([
    (1/12) * m_r * (3*R**2 + L**2),
    (1/12) * m_r * (3*R**2 + L**2),
    (1/2)  * m_r * R**2,
])

# ── Chaser Initial State ──────────────────────────────────────────────────────
delta_nu = -np.deg2rad(0.01)   # ~300 m behind

e  = Target_kep.e
nu = Target_kep.nu - delta_nu
E  = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2),
                     np.sqrt(1 + e) * np.cos(nu / 2))
M  = E - e * np.sin(E)
chaser_kep             = replace(Target_kep, nu=nu, M=M)
chaser_pos, chaser_vel = keplerian_to_cartesian(chaser_kep)

# Initial chaser attitude: x-axis aligned with target x-axis for insertion.
# The chaser will approach along the target cylinder axis from the nozzle end,
# so we align body x-axes from the start.
# At t=0 the target attitude is identity [0,0,0,1], so target x = inertial x.
# We initialise to identity as well — guidance will correct attitude rapidly.

CHASER_EP      = np.array([0., 0., 0., 1.])

CHASER_ANG_VEL = np.array([0., 0., 0.])

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

tau_b   = 1e9;    tau_s   = 1e9;    tau_o   = 1e9
sigma_b = 1e-8;   sigma_s = 1e-8;   sigma_o = 1e-8

# ── Simulation Window ─────────────────────────────────────────────────────────
SIM_START = 0
SIM_END   = 0.11 * 60 * 60    # s — 30 min
SIM_DT    = 1.0               # s
N         = int(SIM_END / SIM_DT)

# ── Q Matrix Process Noise (feeds EKF, not truth dynamics) ───────────────────
sigma_vel_T   = 1e-4 * np.ones(3)   # km/s
sigma_omega_T = 1e-4 * np.ones(3)   # rad/s
sigma_vel_C   = 1e-3 * np.ones(3)   # km/s
sigma_omega_C = 1e-4 * np.ones(3)   # rad/s

# ── EKF Initial 1-sigma Uncertainties ────────────────────────────────────────
sigma_r_T  = 0.01    # km  (10 m)
sigma_v_T  = 1e-3     # km/s
sigma_th_T = 0.01     # rad
sigma_w_T  = 0.001    # rad/s

sigma_r_C  = 0.01    # km  (5 m)
sigma_v_C  = 1e-3     # km/s
sigma_th_C = 0.001    # rad

sigma_bw   = 1e-4
sigma_Ss   = 1e-3
sigma_Oo   = 1e-3

# ── Rendezvous Tolerances ─────────────────────────────────────────────────────
POS_TOL = 0.002    # km  (1 m) — final insertion point
ATT_TOL = 1      # rad — tighter for aligned insertion

# ── Sensor Parameters ─────────────────────────────────────────────────────────

# Gyro
DT_GYRO = 1.0
f_w     = 3e-6

# Star tracker — chaser attitude
DT_STAR_TRACKER = 1.0
sigma_st        = np.deg2rad(0.01)
R_STAR_TRACKER  = (sigma_st**2) * np.eye(3)

# GPS — chaser absolute position
DT_GPS    = 1.0
sigma_gps = 1e-3      # km (10 m)
R_GPS     = (sigma_gps**2) * np.eye(3)

# Combined LIDAR + optical camera → relative position r^T - r^C
DT_REL_POS      = 5.0
sigma_lidar     = 1e-4                   # km (1 m)
sigma_opt_angle = np.deg2rad(0.1e-3)     # rad (0.1 mrad)
_nominal_range  = 36e-3
_sigma_rel_nom  = np.sqrt(sigma_lidar**2 + (_nominal_range * sigma_opt_angle)**2)
R_REL_POS       = (_sigma_rel_nom**2) * np.eye(3)

# Target attitude — feature-tracking proxy
DT_TARGET_ATT    = 10.0
sigma_target_att = np.deg2rad(0.1)
R_TARGET_ATT     = (sigma_target_att**2) * np.eye(3)

# ── Guidance Parameters — Nozzle Insertion ───────────────────────────────────
#
# Geometry (all in target body frame, km):
#
#   Target cylinder axis = body x-axis
#   Cylinder runs from x = -L/2 = -0.004 km  (nozzle end)
#                      to x = +L/2 = +0.004 km  (top end)
#
#   Approach strategy (two-phase):
#     Phase 1 — Approach:  hold at r_dock_T_standoff (9 m outside nozzle face)
#     Phase 2 — Insert:    drive to r_dock_T_insert  (at nozzle face)
#     Switch range: < INSERTION_RANGE_KM
#
#   Chaser docking probe tip at +x face (r_attach_C = [0.5, 0, 0] mm)
#   Insertion complete when chaser probe tip reaches nozzle plane.

# Phase 1 standoff: 5 m outside the nozzle face (-4m - 5m = -9m along x)
r_dock_T_standoff = np.array([-9.0, 0.0, 0.0]) * 1e-3   # km

# Phase 2 insertion: at the nozzle face (-4 m along target x)
r_dock_T_insert   = np.array([-4.0, 0.0, 0.0]) * 1e-3   # km

# Active docking point — updated dynamically in sim loop based on range
r_dock_T    = r_dock_T_standoff.copy()    # start in standoff mode

# Range threshold to switch from standoff to insertion [km]
INSERTION_RANGE_KM = 0.010    # 10 m

# Chaser probe tip in chaser body frame (front face along +x)
r_attach_C  = np.array([0.5, 0.0, 0.0]) * 1e-3   # km

r_rel_des_D = np.zeros(3)
v_rel_des_D = np.zeros(3)
q_D_T       = np.array([0, 0, 0, 1], dtype=float)   # docking frame = target frame
w_des_C     = np.zeros(3)