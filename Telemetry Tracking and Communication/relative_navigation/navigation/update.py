"""
EKF Measurement Updates
=======================
  1. ekf_update_star_tracker  — chaser attitude (quaternion) from star tracker
  2. ekf_update_optical_cam   — angles-only (az/el) from optical camera
  3. ekf_update_gps           — chaser inertial position from GPS

State vector layout (32 nominal / 30 error-state):

  x_T  [0:13]   Target:  r^T(3), v^T(3), q^{I→T}(4), ω^T(3)
  x_C  [13:23]  Chaser:  r^C(3), v^C(3), q^{I→C}(4)
  x_P  [23:32]  Params:  b^C_ω(3), S_S(3), O_O(3)

Error-state (30-dim):
  δx_T [0:12]   δr^T(3), δv^T(3), δθ^T(3), δω^T(3)
  δx_C [12:21]  δr^C(3), δv^C(3), δθ^C(3)
  δx_P [21:30]  δb_ω(3), δS_S(3), δO_O(3)

Quaternion convention: q = [q1, q2, q3, q4]  (q4 = scalar)
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

def skew(v):
    return np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])

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
    xp[IDX_qT] = normalise_quat(quat_mult(x[IDX_qT], delta_theta_to_quat(delta_x[IDX_E_thT])))
    xp[IDX_qC] = normalise_quat(quat_mult(x[IDX_qC], delta_theta_to_quat(delta_x[IDX_E_thC])))
    return xp

def _joseph_update(P, K, H, R):
    IKH = np.eye(N_ERR) - K @ H
    return IKH @ P @ IKH.T + K @ R @ K.T

def _kalman_gain(P, H, R):
    S = H @ P @ H.T + R
    return P @ H.T @ np.linalg.solve(S, np.eye(R.shape[0])).T


# ── 1. STAR TRACKER ──────────────────────────────────────────────────────────

def ekf_update_star_tracker(x_est, P, z_st, R_st):
    """
    EKF update for star tracker quaternion measurement.

    z_st  : (4,)   measured q^{I→C}
    R_st  : (3,3)  attitude noise covariance (rotation-vector space)
    """
    q_pred = x_est[IDX_qC]
    dq     = quat_mult(normalise_quat(z_st), quat_inv(q_pred))
    if dq[3] < 0:
        dq = -dq
    innov = 2.0 * dq[0:3]              # (3,) rotation-vector innovation

    H = np.zeros((3, N_ERR))
    H[:, IDX_E_thC] = np.eye(3)        # only δθ^C observed

    K      = _kalman_gain(P, H, R_st)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_st)

    return x_plus, P_plus, innov


# ── 2. OPTICAL CAMERA (angles-only) ──────────────────────────────────────────

def _h_opt(x):
    """Predicted [azimuth, elevation] from state x."""
    rho_I = x[IDX_rT] - x[IDX_rC]
    rho_C = quat2dcm(x[IDX_qC]) @ rho_I
    az = np.arctan2(rho_C[1], rho_C[0])
    el = np.arctan2(rho_C[2], np.sqrt(rho_C[0]**2 + rho_C[1]**2))
    return np.array([az, el])

def _H_opt(x):
    """Analytic Jacobian of optical camera measurement (2×30)."""
    rho_I = x[IDX_rT] - x[IDX_rC]
    C_IC  = quat2dcm(x[IDX_qC])
    rho_C = C_IC @ rho_I
    rx, ry, rz = rho_C
    d2  = rx**2 + ry**2
    d   = np.sqrt(d2)
    rn2 = rx**2 + ry**2 + rz**2
    if d < 1e-10:
        raise ValueError(f"Azimuth singularity: LOS nearly vertical (d={d:.2e})")
    dz_drho = np.array([
        [-ry/d2,           rx/d2,           0.0   ],
        [-rx*rz/(rn2*d),  -ry*rz/(rn2*d),   d/rn2 ],
    ])
    H = np.zeros((2, N_ERR))
    H[:, IDX_E_rT]  =  dz_drho @ C_IC
    H[:, IDX_E_rC]  = -dz_drho @ C_IC
    H[:, IDX_E_thC] =  dz_drho @ (skew(rho_C) @ C_IC)
    return H


# ── 3. GPS (chaser inertial position) ────────────────────────────────────────

def ekf_update_gps(x_est, P, z_gps, R_gps):
    innov = z_gps - x_est[IDX_rC]

    H = np.zeros((3, N_ERR))
    H[:, IDX_E_rC] = np.eye(3)

    # Innovation gate — reject if > 5-sigma
    S = H @ P @ H.T + R_gps
    mahal = innov @ np.linalg.solve(S, innov)
    if mahal > 25.0:
        return x_est, P, innov

    K      = _kalman_gain(P, H, R_gps)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_gps)
    return x_plus, P_plus, innov

def ekf_update_gps_target(x_est, P, z_gps, R_gps):
    """GPS update for TARGET position — diagnostic only."""
    innov = z_gps - x_est[IDX_rT]
    H = np.zeros((3, N_ERR))
    H[:, IDX_E_rT] = np.eye(3)
    K      = _kalman_gain(P, H, R_gps)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_gps)
    return x_plus, P_plus, innov

def ekf_update_lidar(x_est, P, z_range, R_range):
    """
    LIDAR range measurement: z = ||r^T - r^C|| + noise

    H is (1×30):  ∂z/∂δr_T = +rho_hat,  ∂z/∂δr_C = -rho_hat
    """
    rho_I   = x_est[IDX_rT] - x_est[IDX_rC]
    rho     = np.linalg.norm(rho_I)
    rho_hat = rho_I / rho                       # unit LOS vector

    innov = z_range - rho                       # scalar

    H = np.zeros((1, N_ERR))
    H[0, IDX_E_rT] =  rho_hat
    H[0, IDX_E_rC] = -rho_hat

    if not _innovation_gate(innov.reshape(1), P, H, R_range.reshape(1,1)):
        return x_est, P, innov

    K      = _kalman_gain(P, H, R_range.reshape(1,1))
    x_plus = apply_correction(x_est, K @ innov.reshape(1))
    P_plus = _joseph_update(P, K, H, R_range.reshape(1,1))
    return x_plus, P_plus, innov

def _innovation_gate(innov, P, H, R, n_sigma=5.0):
    """Returns True if innovation passes the chi-squared gate."""
    S = H @ P @ H.T + R
    mahal = innov @ np.linalg.solve(S, innov)
    return mahal <= n_sigma**2

def ekf_update_optical_cam(x_est, P, z_opt, R_opt):
    innov    = z_opt - _h_opt(x_est)
    innov[0] = (innov[0] + np.pi) % (2*np.pi) - np.pi

    # Innovation gate — reject if innovation is implausibly large (5-sigma)
    S = _H_opt(x_est) @ P @ _H_opt(x_est).T + R_opt
    mahal = innov @ np.linalg.solve(S, innov)   # Mahalanobis distance²
    if mahal > 25.0:                             # 5-sigma gate
        return x_est, P, innov                  # skip update

    H      = _H_opt(x_est)
    K      = _kalman_gain(P, H, R_opt)
    x_plus = apply_correction(x_est, K @ innov)
    P_plus = _joseph_update(P, K, H, R_opt)
    return x_plus, P_plus, innov