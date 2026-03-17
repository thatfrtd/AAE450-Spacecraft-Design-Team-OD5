"""
EKF Measurement Updates
=======================
Implements two measurement update functions:

  1. ekf_update_star_tracker   — chaser attitude from star tracker
  2. ekf_measurement_update    — angles-only (azimuth/elevation) from optical camera

State vector layout (32 nominal / 30 error-state):

  x_T  [0:13]   Target:   r^T(3), v^T(3), q^{I→T}(4), ω^T(3)
  x_C  [13:23]  Chaser:   r^C(3), v^C(3), q^{I→C}(4)
  x_P  [23:32]  Params:   b^C_ω(3), S_S(3), O_O(3)

Error-state (30-dim, quaternions replaced by 3-param δθ rotation vector):

  δx_T [0:12]   δr^T(3), δv^T(3), δθ^T(3), δω^T(3)
  δx_C [12:21]  δr^C(3), δv^C(3), δθ^C(3)
  δx_P [21:30]  δb_ω(3), δS_S(3), δO_O(3)

Quaternion convention: q = [q1, q2, q3, q4]  (q4 = scalar part)
  v_body = C(q) @ v_inertial
"""

import numpy as np
# from dynamics.attitude import *

# ─────────────────────────────────────────────────────────────────────────────
# Nominal-state slices (into 32-element x vector)
# ─────────────────────────────────────────────────────────────────────────────
IDX_rT = slice(0,   3)
IDX_vT = slice(3,   6)
IDX_qT = slice(6,   10)
IDX_wT = slice(10,  13)
IDX_rC = slice(13,  16)
IDX_vC = slice(16,  19)
IDX_qC = slice(19,  23)
IDX_bw = slice(23,  26)
IDX_Ss = slice(26,  29)
IDX_Oo = slice(29,  32)

# ─────────────────────────────────────────────────────────────────────────────
# Error-state slices (into 30-element δx vector)
# ─────────────────────────────────────────────────────────────────────────────
IDX_E_rT  = slice(0,  3)
IDX_E_vT  = slice(3,  6)
IDX_E_thT = slice(6,  9)    # δθ^T  (rotation vector)
IDX_E_wT  = slice(9,  12)
IDX_E_rC  = slice(12, 15)
IDX_E_vC  = slice(15, 18)
IDX_E_thC = slice(18, 21)   # δθ^C  (rotation vector)
IDX_E_bw  = slice(21, 24)
IDX_E_Ss  = slice(24, 27)
IDX_E_Oo  = slice(27, 30)

N_ERR   = 30
N_STATE = 32



def quat2dcm(q: np.ndarray) -> np.ndarray:
    """q = [q1,q2,q3,q4] (q4 scalar) → 3×3 DCM  v_body = C @ v_inertial."""
    q1, q2, q3, q4 = q
    return np.array([
        [1-2*q2**2-2*q3**2,    2*(q1*q2+q3*q4),    2*(q1*q3-q2*q4)],
        [2*(q1*q2-q3*q4),    1-2*q1**2-2*q3**2,    2*(q2*q3+q1*q4)],
        [2*(q1*q3+q2*q4),    2*(q2*q3-q1*q4),    1-2*q1**2-2*q2**2],
    ])
 
 
def quat_multiply(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Hamilton product p ⊗ q, convention [vec; scalar]."""
    p1, p2, p3, p4 = p
    q1, q2, q3, q4 = q
    return np.array([
         p4*q1 + p3*q2 - p2*q3 + p1*q4,
        -p3*q1 + p4*q2 + p1*q3 + p2*q4,
         p2*q1 - p1*q2 + p4*q3 + p3*q4,
        -p1*q1 - p2*q2 - p3*q3 + p4*q4,
    ])
 
 
def quat_inv(q: np.ndarray) -> np.ndarray:
    """Unit-quaternion inverse = conjugate: [-vec; scalar]."""
    return np.array([-q[0], -q[1], -q[2], q[3]])
 
 
def skew(v: np.ndarray) -> np.ndarray:
    """3×3 skew-symmetric matrix: skew(v) @ w = v × w."""
    return np.array([
        [ 0,    -v[2],  v[1]],
        [ v[2],  0,    -v[0]],
        [-v[1],  v[0],  0   ],
    ])
 
 
def delta_theta_to_quat(dtheta: np.ndarray) -> np.ndarray:
    """Small rotation vector δθ → unit quaternion [vec; scalar]."""
    angle = np.linalg.norm(dtheta)
    if angle < 1e-10:
        return np.array([0.0, 0.0, 0.0, 1.0])
    axis = dtheta / angle
    s = np.sin(angle / 2.0)
    return np.concatenate([s * axis, [np.cos(angle / 2.0)]])
 
 
def normalize_quat(q: np.ndarray) -> np.ndarray:
    return q / np.linalg.norm(q)
 

# ─────────────────────────────────────────────────────────────────────────────
# Error-state correction  (shared by both update steps)
# ─────────────────────────────────────────────────────────────────────────────

def apply_correction(x: np.ndarray, delta_x: np.ndarray) -> np.ndarray:
    """
    Apply error-state correction δx to nominal state x.

    All states are corrected additively except quaternions, which use a
    right-multiplicative update:  q_new = q_nom ⊗ δq(δθ)
    """
    xp = x.copy()

    # Additive corrections
    xp[IDX_rT] += delta_x[IDX_E_rT]
    xp[IDX_vT] += delta_x[IDX_E_vT]
    xp[IDX_wT] += delta_x[IDX_E_wT]
    xp[IDX_rC] += delta_x[IDX_E_rC]
    xp[IDX_vC] += delta_x[IDX_E_vC]
    xp[IDX_bw] += delta_x[IDX_E_bw]
    xp[IDX_Ss] += delta_x[IDX_E_Ss]
    xp[IDX_Oo] += delta_x[IDX_E_Oo]

    # Multiplicative attitude corrections  q_new = q_nom ⊗ δq(δθ)
    xp[IDX_qT] = normalize_quat(
        quat_multiply(x[IDX_qT], delta_theta_to_quat(delta_x[IDX_E_thT]))
    )
    xp[IDX_qC] = normalize_quat(
        quat_multiply(x[IDX_qC], delta_theta_to_quat(delta_x[IDX_E_thC]))
    )

    return xp


# ─────────────────────────────────────────────────────────────────────────────
# ① STAR TRACKER — attitude update
# ─────────────────────────────────────────────────────────────────────────────

def ekf_update_star_tracker(
    x_est: np.ndarray,
    P:     np.ndarray,
    z_st:  np.ndarray,
    R_st:  np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    EKF measurement update for the star tracker.

    The star tracker directly measures the chaser inertial attitude as a
    quaternion.  Because the EKF error-state for attitude is a 3-element
    rotation vector δθ, the innovation is expressed as a rotation vector
    rather than a 4-element quaternion difference.

    Measurement model
    -----------------
      z_st  = q^{I→C}  (unit quaternion, 4-element)
      h(x)  = x[IDX_qC]

    Innovation (rotation-vector form, 3-element)
    --------------------------------------------
      δq    = q_meas ⊗ q_pred^{-1}          (attitude error quaternion)
      innov = 2 · δq[0:3]                   (small-angle rotation vector)

    The sign of δq is forced so the scalar part δq[3] > 0, ensuring we
    always take the short arc of the rotation.

    Jacobian
    --------
      H = ∂innov / ∂δx  →  (3 × 30)
      Only the δθ^C block is nonzero:  H[:, IDX_E_thC] = I₃

    Parameters
    ----------
    x_est : (32,)   nominal state estimate
    P     : (30,30) error-state covariance
    z_st  : (4,)    measured quaternion  q^{I→C}  [q1,q2,q3,q4]
    R_st  : (3,3)   attitude measurement noise covariance (rotation-vector space)

    Returns
    -------
    x_plus : (32,)   updated state
    P_plus : (30,30) updated covariance (Joseph form)
    innov  : (3,)    rotation-vector innovation [rad]
    """
    # ── Predicted and measured attitude
    q_pred = x_est[IDX_qC]
    q_meas = normalize_quat(z_st)

    # ── Quaternion error: δq = q_meas ⊗ q_pred^{-1}
    dq = quat_multiply(q_meas, quat_inv(q_pred))

    # Force short arc (scalar part positive)
    if dq[3] < 0.0:
        dq = -dq

    # Rotation-vector innovation:  innov ≈ 2 · δq_vec
    innov = 2.0 * dq[0:3]                         # (3,)

    # ── Measurement Jacobian  H : (3 × 30)
    #    h(x) = q^C  →  ∂(δθ innov)/∂δθ^C = I₃
    #    All other error-state blocks → zero
    H = np.zeros((3, N_ERR))
    H[:, IDX_E_thC] = np.eye(3)

    # ── Kalman update
    S = H @ P @ H.T + R_st                        # (3×3)
    K = P @ H.T @ np.linalg.solve(S, np.eye(3)).T # (30×3)

    delta_x = K @ innov                            # (30,)
    x_plus  = apply_correction(x_est, delta_x)

    # Joseph form for numerical stability
    IKH    = np.eye(N_ERR) - K @ H
    P_plus = IKH @ P @ IKH.T + K @ R_st @ K.T

    return x_plus, P_plus, innov


# ─────────────────────────────────────────────────────────────────────────────
# ② OPTICAL CAMERA — angles-only update
# ─────────────────────────────────────────────────────────────────────────────

def h_opt(x: np.ndarray) -> np.ndarray:
    """
    Nonlinear optical camera measurement  z = [azimuth, elevation].

    Algorithm
    ---------
    1. LOS in inertial frame:  ρ^I = r^T − r^C
    2. Rotate to chaser body:  ρ^C = C(q^{I→C}) @ ρ^I
    3. Compute angles from ρ^C components
    """
    rho_I = x[IDX_rT] - x[IDX_rC]
    C_IC  = quat2dcm(x[IDX_qC])
    rho_C = C_IC @ rho_I

    az = np.arctan2(rho_C[1], rho_C[0])
    el = np.arctan2(rho_C[2], np.sqrt(rho_C[0]**2 + rho_C[1]**2))

    return np.array([az, el])


def meas_Jacobian(x: np.ndarray) -> np.ndarray:
    """
    Analytic measurement Jacobian  H = ∂z/∂(δx)  (2 × 30).

    Nonzero blocks
    --------------
      ∂z/∂δr^T  =  (∂z/∂ρ^C) C_{IC}
      ∂z/∂δr^C  = −(∂z/∂ρ^C) C_{IC}
      ∂z/∂δθ^C  =  (∂z/∂ρ^C) [ρ^C]× C_{IC}

    Attitude perturbation (right-multiplicative):
      ρ^C_new = ρ^C + [ρ^C]× C_{IC} δθ^C
      ⟹  ∂ρ^C/∂δθ^C = [ρ^C]× C_{IC}

    All other blocks (δv^T, δθ^T, δω^T, δv^C, δb_ω, δS, δO) → zero.
    """
    rho_I = x[IDX_rT] - x[IDX_rC]       # ← r^T minus r^C (target relative to chaser)
    C_IC  = quat2dcm(x[IDX_qC])          # ← use IDX_qC slice, not raw index
    rho_C = C_IC @ rho_I

    rx, ry, rz = rho_C
    d2  = rx**2 + ry**2
    d   = np.sqrt(d2)
    rn2 = rx**2 + ry**2 + rz**2

    if d < 1e-10:
        raise ValueError(
            f"Azimuth singularity: LOS nearly vertical (d={d:.2e}). "
            "Consider switching angle parametrisation near the pole."
        )

    # ∂[az, el] / ∂ρ^C  →  (2×3)
    dz_drho = np.array([
        [-ry / d2,            rx / d2,            0.0     ],
        [-rx*rz / (rn2 * d), -ry*rz / (rn2 * d), d / rn2 ],
    ])

    dz_drT  =  dz_drho @ C_IC                  # (2×3)
    dz_drC  = -dz_drT                           # (2×3)
    dz_dthC =  dz_drho @ (skew(rho_C) @ C_IC)  # (2×3)

    H = np.zeros((2, N_ERR))
    H[:, IDX_E_rT]  = dz_drT
    H[:, IDX_E_rC]  = dz_drC
    H[:, IDX_E_thC] = dz_dthC

    return H


def ekf_measurement_update(
    x_est:  np.ndarray,
    P:      np.ndarray,
    z_opt:  np.ndarray,
    R_opt:  np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    EKF measurement update for the optical camera (angles-only).

    Parameters
    ----------
    x_est : (32,)   nominal state estimate
    P     : (30,30) error-state covariance
    z_opt : (2,)    [azimuth, elevation] measurement [rad]
    R_opt : (2,2)   optical camera noise covariance

    Returns
    -------
    x_plus : (32,)   updated state
    P_plus : (30,30) updated covariance (Joseph form)
    innov  : (2,)    innovation [az_err, el_err] in radians
    """
    z_pred = h_opt(x_est)
    innov  = z_opt - z_pred

    # Wrap azimuth innovation to (−π, π] to avoid 2π jumps
    innov[0] = (innov[0] + np.pi) % (2.0 * np.pi) - np.pi

    H = meas_Jacobian(x_est)                       # (2×30)
    S = H @ P @ H.T + R_opt                        # (2×2)
    K = P @ H.T @ np.linalg.solve(S, np.eye(2)).T  # (30×2)

    delta_x = K @ innov
    x_plus  = apply_correction(x_est, delta_x)

    # Joseph form
    IKH    = np.eye(N_ERR) - K @ H
    P_plus = IKH @ P @ IKH.T + K @ R_opt @ K.T

    return x_plus, P_plus, innov


# ─────────────────────────────────────────────────────────────────────────────
# Quick self-tests
# ─────────────────────────────────────────────────────────────────────────────

def _make_test_state() -> np.ndarray:
    x = np.zeros(N_STATE)
    x[IDX_rT] = np.array([100.0,  5.0, -3.0])
    x[IDX_vT] = np.array([ -0.1,  0.01, 0.0])
    ang = np.deg2rad(30)
    x[IDX_qT] = np.array([0.0, 0.0, np.sin(ang/2), np.cos(ang/2)])
    x[IDX_wT] = np.array([0.01, -0.005, 0.002])
    x[IDX_rC] = np.array([0.0, 0.0, 0.0])
    x[IDX_vC] = np.array([0.0, 0.0, 0.0])
    ang2 = np.deg2rad(15)
    x[IDX_qC] = np.array([0.0, np.sin(ang2/2), 0.0, np.cos(ang2/2)])
    x[IDX_bw] = np.array([ 1e-4, -2e-4, 1e-5])
    x[IDX_Ss] = np.array([ 1.0,   1.0,  1.0 ])
    x[IDX_Oo] = np.array([ 0.0,   0.0,  0.0 ])
    return x


def test_star_tracker_reduces_covariance():
    x = _make_test_state()
    P = np.eye(N_ERR)
    sigma = np.deg2rad(0.01)
    R_st  = (sigma**2) * np.eye(3)
    z_true  = x[IDX_qC].copy()
    z_noisy = normalize_quat(z_true + np.array([1e-4, -1e-4, 5e-5, 0.0]))
    x_p, P_p, innov = ekf_update_star_tracker(x, P, z_noisy, R_st)
    assert np.trace(P_p) < np.trace(P), "Covariance did not decrease (star tracker)"
    print(f"✓ star tracker:  trace(P) {np.trace(P):.4f} → {np.trace(P_p):.4f}")
    print(f"  innovation = {np.rad2deg(innov)} deg")


def test_optical_camera_reduces_covariance():
    x = _make_test_state()
    P = np.eye(N_ERR)
    sigma = np.deg2rad(0.1)
    R_opt = (sigma**2) * np.eye(2)
    z_true  = h_opt(x)
    z_noisy = z_true + np.array([1e-4, -5e-5])
    x_p, P_p, innov = ekf_measurement_update(x, P, z_noisy, R_opt)
    assert np.trace(P_p) < np.trace(P), "Covariance did not decrease (optical)"
    print(f"✓ optical camera: trace(P) {np.trace(P):.4f} → {np.trace(P_p):.4f}")
    print(f"  innovation = {np.rad2deg(innov)} deg")


def test_zero_innovation_star_tracker():
    x = _make_test_state()
    P = np.eye(N_ERR) * 0.01
    R_st = (np.deg2rad(0.01)**2) * np.eye(3)
    z_perfect = x[IDX_qC].copy()
    x_p, P_p, innov = ekf_update_star_tracker(x, P, z_perfect, R_st)
    assert np.allclose(innov, 0.0, atol=1e-12), "Non-zero innovation on perfect star tracker meas"
    print("✓ zero-innovation star tracker: state unchanged")


if __name__ == "__main__":
    np.random.seed(0)
    print("=" * 60)
    print("Measurement Update Tests")
    print("=" * 60)
    test_star_tracker_reduces_covariance()
    test_optical_camera_reduces_covariance()
    test_zero_innovation_star_tracker()
    print("\nAll tests passed ✓")