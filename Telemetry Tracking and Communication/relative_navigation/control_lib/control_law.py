"""
Control Law — matches Woffinden & Geller JGCD 2007, Eqs. 51–54.
All quantities in km and km/s.

Eq 51:  Δv_rqd  = K_r(r_des − r_C) + K_v(v_des − v_C)   [inertial, km/s²]
Eq 52:  Δv_cmd  = T(q^{I→C}) @ Δv_rqd                    [body frame]
Eq 53:  τ_cmd   = K_θ · θ_C_des + K_ω · (ω_des − ω_C)   [kg·km²/s²]
Eq 54:  θ_C_des = q_des[0:3] − q^{C→I}[0:3]

Gain tuning notes
-----------------
Translational (km units):
  ω_n = 0.05 rad/s, ζ = 1  →  K_r = ω_n² = 2.5e-3 [1/s²]
                                K_v = 2ζω_n = 0.1    [1/s]

Attitude (must produce physically reasonable torques):
  Target angular acceleration α ~ 0.01 rad/s²
  I_c (km units) ~ 1/6 * 1000 * (1e-3)² = 1.67e-7 kg·km²
  Required τ ~ I_c * α ~ 1.67e-9 kg·km²/s²

  PD gains are therefore proportional to I_c:
    K_θ = I_c * ω_n_att²        [kg·km²/s²  per rad]
    K_ω = I_c * 2ζω_n_att       [kg·km²/s²  per rad/s]
  with ω_n_att = 0.05 rad/s, ζ = 1.
"""
import numpy as np
from dynamics.attitude import *
import config


def control_law(x_est, w_C, r_des, v_des, q_des, w_des, I_c):
    """
    Parameters
    ----------
    x_est   : (32,) estimated state (km, km/s, quaternions, rad/s)
    w_C     : (3,)  true angular velocity from gyro [rad/s]
    r_des   : (3,)  desired position  [km, inertial]
    v_des   : (3,)  desired velocity  [km/s, inertial]
    q_des   : (4,)  desired attitude  q^{I→C_des}
    w_des   : (3,)  desired angular velocity [rad/s]
    I_c     : (3,3) chaser MOI [kg·km²]

    Returns
    -------
    u_cmd   : (3,)  Δv command [km/s², chaser body frame]    Eq 52
    tau_cmd : (3,)  torque command [kg·km²/s²]               Eq 53
    """
    r_C = x_est[13:16]
    v_C = x_est[16:19]
    q_C = x_est[19:23]

    # ── Translational gains ───────────────────────────────────────────────────
    # ω_n = 0.002 rad/s → period ~3140 s, max accel ~4e-5 km/s² at 100m range
    # This keeps delta-v per step at ~4e-4 km/s = 0.4 mm/s — safe for proximity
    omega_n = 0.02
    zeta    = 3.0
    K_r = (omega_n**2)          * np.eye(3)   # 4e-6  [1/s²]
    K_v = (2 * zeta * omega_n)  * np.eye(3)   # 4e-3  [1/s]

    # ── Attitude gains — scaled by MOI so torques are physically sized ────────
    # Using mean diagonal of I_c as representative scalar MOI
    I_scalar    = np.mean(np.diag(I_c))        # [kg·km²]
    omega_n_att = 0.05                         # rad/s
    zeta_att    = 1.0
    K_theta = I_scalar * (omega_n_att**2)              * np.eye(3)
    K_w     = I_scalar * (2 * zeta_att * omega_n_att)  * np.eye(3)

    # ── Eq 51: required Δv in INERTIAL frame ─────────────────────────────────
    dv_rqd = K_r @ (r_des - r_C) + K_v @ (v_des - v_C)

    # ── Eq 52: rotate to CHASER BODY frame ───────────────────────────────────
    T_IC  = quat2dcm(q_C)
    u_cmd = T_IC @ dv_rqd

    # ── Eq 54: attitude error ─────────────────────────────────────────────────
    q_C_to_I    = quat_inv(q_C)
    theta_C_des = q_des[0:3] - q_C_to_I[0:3]

    # ── Eq 53: torque command ─────────────────────────────────────────────────
    tau_cmd = K_theta @ theta_C_des + K_w @ (w_des - w_C)

    return u_cmd, tau_cmd