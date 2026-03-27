"""
EKF Measurement Updates
=======================
Update functions and their measurement types:

  ekf_update_star_tracker      (4,)  → chaser attitude
  ekf_update_relative_position (3,)  → r^T - r^C [km, inertial]
  ekf_update_gps               (3,)  → r^C [km, inertial]
  ekf_update_target_attitude   (4,)  → target attitude

x_est layout (32 — nominal state):
  [0:3]  r^T  [3:6]  v^T  [6:10]  q^{I→T}  [10:13] ω^T
  [13:16] r^C [16:19] v^C [19:23] q^{I→C}
  [23:26] b_w [26:29] S_s [29:32] O_o

Error-state layout (30 — used in P, H, K):
  [0:3]  δr^T [3:6] δv^T [6:9]  δθ^T [9:12]  δω^T
  [12:15] δr^C [15:18] δv^C [18:21] δθ^C
  [21:24] δb_w [24:27] δS_s [27:30] δO_o
"""
import numpy as np

# ── Nominal-state slices ──────────────────────────────────────────────────────
IDX_rT = slice(0,   3);   IDX_vT = slice(3,   6)
IDX_qT = slice(6,  10);   IDX_wT = slice(10, 13)
IDX_rC = slice(13, 16);   IDX_vC = slice(16, 19)
IDX_qC = slice(19, 23)
IDX_bw = slice(23, 26);   IDX_Ss = slice(26, 29);   IDX_Oo = slice(29, 32)

# ── Error-state slices ────────────────────────────────────────────────────────
IDX_E_rT  = slice(0,  3);   IDX_E_vT  = slice(3,  6)
IDX_E_thT = slice(6,  9);   IDX_E_wT  = slice(9,  12)
IDX_E_rC  = slice(12, 15);  IDX_E_vC  = slice(15, 18)
IDX_E_thC = slice(18, 21)
IDX_E_bw  = slice(21, 24);  IDX_E_Ss  = slice(24, 27);  IDX_E_Oo  = slice(27, 30)

N_ERR   = 30
N_STATE = 32


# ── Math utilities ────────────────────────────────────────────────────────────

def quat2dcm(q):
    q1, q2, q3, q4 = q
    return np.array([
        [1-2*q2**2-2*q3**2,  2*(q1*q2+q3*q4),   2*(q1*q3-q2*q4)],
        [2*(q1*q2-q3*q4),    1-2*q1**2-2*q3**2,  2*(q2*q3+q1*q4)],
        [2*(q1*q3+q2*q4),    2*(q2*q3-q1*q4),    1-2*q1**2-2*q2**2],
    ])

def quat_mult(p, q):
    p1,p2,p3,p4 = p;  q1,q2,q3,q4 = q
    return np.array([
         p4*q1+p3*q2-p2*q3+p1*q4,
        -p3*q1+p4*q2+p1*q3+p2*q4,
         p2*q1-p1*q2+p4*q3+p3*q4,
        -p1*q1-p2*q2-p3*q3+p4*q4,
    ])

def quat_inv(q):
    return np.array([-q[0], -q[1], -q[2], q[3]])

def delta_theta_to_quat(dtheta):
    angle = np.linalg.norm(dtheta)
    if angle < 1e-10:
        return np.array([0., 0., 0., 1.])
    axis = dtheta / angle
    return np.concatenate([np.sin(angle/2)*axis, [np.cos(angle/2)]])

def normalise_quat(q):
    return q / np.linalg.norm(q)


# ── Shared helpers ────────────────────────────────────────────────────────────

def apply_correction(x, delta_x):
    """Apply 30-element error-state correction to 32-element nominal state."""
    xp = x.copy()
    xp[IDX_rT] += delta_x[IDX_E_rT];  xp[IDX_vT] += delta_x[IDX_E_vT]
    xp[IDX_wT] += delta_x[IDX_E_wT]
    xp[IDX_rC] += delta_x[IDX_E_rC];  xp[IDX_vC] += delta_x[IDX_E_vC]
    xp[IDX_bw] += delta_x[IDX_E_bw];  xp[IDX_Ss] += delta_x[IDX_E_Ss]
    xp[IDX_Oo] += delta_x[IDX_E_Oo]
    xp[IDX_qT] = normalise_quat(quat_mult(x[IDX_qT],
                                 delta_theta_to_quat(delta_x[IDX_E_thT])))
    xp[IDX_qC] = normalise_quat(quat_mult(x[IDX_qC],
                                 delta_theta_to_quat(delta_x[IDX_E_thC])))
    return xp

def _joseph_update(P, K, H, R):
    IKH = np.eye(N_ERR) - K @ H
    return IKH @ P @ IKH.T + K @ R @ K.T

def _kalman_gain(P, H, R):
    S = H @ P @ H.T + R
    return P @ H.T @ np.linalg.solve(S, np.eye(R.shape[0])).T

def _innovation_gate(innov, P, H, R, n_sigma=5.0):
    S     = H @ P @ H.T + R
    mahal = float(innov @ np.linalg.solve(S, innov))
    return mahal <= n_sigma**2


# ── 1. STAR TRACKER ───────────────────────────────────────────────────────────

def ekf_update_star_tracker(x_est, P, z_st, R_st):
    """
    Update chaser attitude from star tracker.
    z_st : (4,) measured q^{I→C}
    R_st : (3,3) rotation-vector noise covariance [rad²]
    """
    q_pred = x_est[IDX_qC]
    dq     = quat_mult(normalise_quat(z_st), quat_inv(q_pred))
    if dq[3] < 0:
        dq = -dq
    innov = 2.0 * dq[0:3]

    H = np.zeros((3, N_ERR))
    H[:, IDX_E_thC] = np.eye(3)

    if not _innovation_gate(innov, P, H, R_st):
        return x_est, P, innov

    K      = _kalman_gain(P, H, R_st)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_st)
    return x_plus, P_plus, innov


# ── 2. RELATIVE POSITION (LIDAR + camera combined) ───────────────────────────

def ekf_update_relative_position(x_est, P, z_rel, R_rel):
    """
    Update from combined LIDAR+camera relative position measurement.

    The LIDAR provides range and the optical camera provides bearing
    angles.  Together they produce a 3D relative position vector
    z = r^T - r^C + noise  expressed in the INERTIAL frame.

    This replaces the separate LIDAR and optical camera updates.
    Advantages:
      - No FOV singularity issues
      - No angle wrapping
      - Simple linear Jacobian
      - Directly constrains what the controller needs (r^T - r^C)

    Measurement model:  z = r^T - r^C + noise
    Jacobian (3×30):
      H[:, IDX_E_rT] = +I₃
      H[:, IDX_E_rC] = -I₃

    Parameters
    ----------
    z_rel : (3,)  measured r^T - r^C [km, inertial]
    R_rel : (3,3) measurement noise covariance [km²]
    """
    innov = z_rel - (x_est[IDX_rT] - x_est[IDX_rC])   # (3,)

    H = np.zeros((3, N_ERR))
    H[:, IDX_E_rT] =  np.eye(3)
    H[:, IDX_E_rC] = -np.eye(3)

    if not _innovation_gate(innov, P, H, R_rel):
        return x_est, P, innov

    K      = _kalman_gain(P, H, R_rel)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_rel)
    return x_plus, P_plus, innov


# ── 3. GPS — chaser absolute position ────────────────────────────────────────

def ekf_update_gps(x_est, P, z_gps, R_gps):
    """
    Update chaser position from GPS.
    z_gps : (3,) measured r^C [km]
    R_gps : (3,3) noise covariance [km²]
    """
    innov = z_gps - x_est[IDX_rC]

    H = np.zeros((3, N_ERR))
    H[:, IDX_E_rC] = np.eye(3)

    if not _innovation_gate(innov, P, H, R_gps):
        return x_est, P, innov

    K      = _kalman_gain(P, H, R_gps)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_gps)
    return x_plus, P_plus, innov


# ── 4. TARGET ATTITUDE ────────────────────────────────────────────────────────

def ekf_update_target_attitude(x_est, P, z_qT, R_target_att):
    """
    Update target attitude from feature-tracking proxy.
    z_qT        : (4,) measured q^{I→T}
    R_target_att: (3,3) rotation-vector noise covariance [rad²]
    """
    q_pred = x_est[IDX_qT]
    dq     = quat_mult(normalise_quat(z_qT), quat_inv(q_pred))
    if dq[3] < 0:
        dq = -dq
    innov = 2.0 * dq[0:3]

    H = np.zeros((3, N_ERR))
    H[:, IDX_E_thT] = np.eye(3)

    if not _innovation_gate(innov, P, H, R_target_att):
        return x_est, P, innov

    K      = _kalman_gain(P, H, R_target_att)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_target_att)
    return x_plus, P_plus, innov