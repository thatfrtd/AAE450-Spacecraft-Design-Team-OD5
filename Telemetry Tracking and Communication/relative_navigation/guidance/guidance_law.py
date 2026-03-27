"""
Guidance Law — Nozzle Insertion, all quantities in km and km/s.

Insertion geometry
------------------
Target cylinder axis = target body x-axis.
Nozzle end is at x = -L/2 = -4 m in the target body frame.

The chaser must:
  1. Approach along the target x-axis from outside the nozzle end
  2. Align its own x-axis with the target x-axis (probe pointing at nozzle)
  3. Insert slowly along the axis

Two-phase guidance:
  Phase 1 (range > INSERTION_RANGE_KM):
    Drive to standoff point 5 m outside nozzle face.
    Attitude: align chaser x-axis with -target_x (pointing at nozzle).

  Phase 2 (range <= INSERTION_RANGE_KM):
    Drive slowly into the nozzle along the axis.
    Same attitude — x-axes aligned.
    Velocity gate tightened for slow insertion.

Velocity gate
-------------
Applied to relative closing speed along LOS only — the orbital velocity
component of v_des is never modified (would corrupt navigation).

    v_max = K_GATE * rho

  Phase 1: K_GATE = 0.01 /s  →  at 300 m: 3 m/s, at 10 m: 0.1 m/s
  Phase 2: K_GATE = 0.002 /s →  at 5 m: 0.01 m/s  (10 mm/s insertion speed)
"""
import numpy as np
from dynamics.attitude import *
import config

K_GATE_APPROACH  = 0.05    # 1/s — phase 1 approach
K_GATE_INSERTION = 0.005   # 1/s — phase 2 slow insertion


def guidance_law(x_est):
    """
    Parameters
    ----------
    x_est : (32,) estimated state vector [km, km/s, quaternions, rad/s]

    Returns
    -------
    r_des : (3,)  desired chaser position        [km]
    v_des : (3,)  desired chaser velocity        [km/s]
    q_des : (4,)  desired chaser attitude        [quaternion]
    w_des : (3,)  desired chaser angular velocity [rad/s]
    """
    r_T = x_est[0:3];   v_T = x_est[3:6]
    q_T = x_est[6:10];  w_T = x_est[10:13]
    r_C = x_est[13:16]; v_C = x_est[16:19]
    q_C = x_est[19:23]

    # Rotation matrices
    T_T_I = quat2dcm(quat_inv(q_T))   # inertial → target body
    T_I_T = quat2dcm(q_T)             # target body → inertial
    T_D_T = quat2dcm(config.q_D_T)    # target body → docking frame
    T_C_I = quat2dcm(quat_inv(q_C))   # inertial → chaser body

    # ── Desired position ──────────────────────────────────────────────────────
    # config.r_dock_T is updated by the sim loop (standoff vs insert)
    r_des = (r_T
             + T_T_I @ (config.r_dock_T + T_D_T @ config.r_rel_des_D)
             - T_C_I @ config.r_attach_C)

    # ── Nominal desired velocity ──────────────────────────────────────────────
    v_des_nominal = (v_T
                     + T_T_I @ (T_D_T @ config.v_rel_des_D
                                + np.cross(w_T, config.r_dock_T
                                           + T_D_T @ config.r_rel_des_D))
                     - T_C_I @ np.cross(config.w_des_C, config.r_attach_C))

    # ── Velocity gate — relative closing speed along LOS only ─────────────────
    rho_vec = r_des - r_C
    rho     = np.linalg.norm(rho_vec)

    # Select gate based on phase
    if rho <= config.INSERTION_RANGE_KM:
        k_gate = K_GATE_INSERTION
    else:
        k_gate = K_GATE_APPROACH

    if rho > 1e-10:
        rho_hat   = rho_vec / rho
        v_max     = k_gate * rho
        v_rel_los = np.dot(v_C - v_des_nominal, rho_hat)
        if v_rel_los > v_max:
            v_des = v_des_nominal - (v_rel_los - v_max) * rho_hat
        else:
            v_des = v_des_nominal
    else:
        v_des = v_des_nominal

    # ── Desired attitude: align chaser x-axis with target x-axis ─────────────
    #
    # For nozzle insertion the chaser x-axis must point INTO the nozzle,
    # i.e. align with the NEGATIVE target x-axis (approaching from outside).
    #
    # target x-axis in inertial frame:
    target_x_I   = T_I_T[:, 0]        # first column of body→inertial DCM

    # Chaser x-axis should point along -target_x (toward nozzle)
    approach_dir = -target_x_I

    x_body    = np.array([1., 0., 0.])
    axis      = np.cross(x_body, approach_dir)
    axis_norm = np.linalg.norm(axis)

    if axis_norm < 1e-10:
        # Already aligned or 180° flip
        if np.dot(x_body, approach_dir) > 0:
            q_des = np.array([0., 0., 0., 1.])
        else:
            q_des = np.array([0., 1., 0., 0.])
    else:
        axis  = axis / axis_norm
        angle = np.arccos(np.clip(np.dot(x_body, approach_dir), -1, 1))
        q_des = np.array([
            axis[0]*np.sin(angle/2), axis[1]*np.sin(angle/2),
            axis[2]*np.sin(angle/2), np.cos(angle/2)
        ])

    # ── Desired angular velocity: chaser matches target rotation ─────────────
    q_T_C = quat_multiply(q_T, q_C)
    w_des  = quat2dcm(q_T_C) @ w_T

    return r_des, v_des, q_des, w_des