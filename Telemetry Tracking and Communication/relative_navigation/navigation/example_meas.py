"""
Angles-Only Navigation Measurement Model
Based on: Woffinden & Geller, "Relative Angles-Only Navigation and Pose
Estimation for Autonomous Orbital Rendezvous," JGCD 2007.

State vector layout (32 nominal / 30 error-state):

  x_T  [0:13]   Target:   r^T(3), v^T(3), q^{I→T}(4), ω^T(3)
  x_C  [13:23]  Chaser:   r^C(3), v^C(3), q^{I→C}(4)
  x_P  [23:32]  Params:   b^C_ω(3), S_S(3), O_O(3)

Error-state (30-dim, quaternions replaced by 3-param rotation-vector δθ):

  δx_T [0:12]   δr^T(3), δv^T(3), δθ^T(3), δω^T(3)
  δx_C [12:21]  δr^C(3), δv^C(3), δθ^C(3)
  δx_P [21:30]  δb_ω(3), δS_S(3), δO_O(3)

Quaternion convention: q = [q1, q2, q3, q4] where q4 is the scalar part.
  v_body = C(q) @ v_inertial   (q represents rotation I → body)
"""

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# Nominal-state slices  (into 32-element x vector)
# ─────────────────────────────────────────────────────────────────────────────
IDX_rT = slice(0,   3)    # target position  (inertial)
IDX_vT = slice(3,   6)    # target velocity  (inertial)
IDX_qT = slice(6,   10)   # target attitude  q^{I→T}
IDX_wT = slice(10,  13)   # target angular velocity
IDX_rC = slice(13,  16)   # chaser position  (inertial)
IDX_vC = slice(16,  19)   # chaser velocity  (inertial)
IDX_qC = slice(19,  23)   # chaser attitude  q^{I→C}
IDX_bw = slice(23,  26)   # gyro bias
IDX_Ss = slice(26,  29)   # scale factor
IDX_Oo = slice(29,  32)   # misalignment

N_STATE = 32

# ─────────────────────────────────────────────────────────────────────────────
# Error-state slices  (into 30-element δx vector)
# ─────────────────────────────────────────────────────────────────────────────
IDX_E_rT  = slice(0,  3)
IDX_E_vT  = slice(3,  6)
IDX_E_thT = slice(6,  9)   # δθ^T  (rotation vector)
IDX_E_wT  = slice(9,  12)
IDX_E_rC  = slice(12, 15)
IDX_E_vC  = slice(15, 18)
IDX_E_thC = slice(18, 21)  # δθ^C  (rotation vector)
IDX_E_bw  = slice(21, 24)
IDX_E_Ss  = slice(24, 27)
IDX_E_Oo  = slice(27, 30)

N_ERR = 30


# ─────────────────────────────────────────────────────────────────────────────
# Math utilities
# ─────────────────────────────────────────────────────────────────────────────

def quat2dcm(q: np.ndarray) -> np.ndarray:
    """
    Quaternion → Direction Cosine Matrix.

    q = [q1, q2, q3, q4]  (q4 = scalar)
    Returns C (3×3) such that v_body = C @ v_inertial.
    """
    q1, q2, q3, q4 = q
    return np.array([
        [1 - 2*q2**2 - 2*q3**2,  2*(q1*q2 + q3*q4),      2*(q1*q3 - q2*q4)],
        [2*(q1*q2 - q3*q4),      1 - 2*q1**2 - 2*q3**2,  2*(q2*q3 + q1*q4)],
        [2*(q1*q3 + q2*q4),      2*(q2*q3 - q1*q4),      1 - 2*q1**2 - 2*q2**2]
    ])


def quat_mult(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    """
    Hamilton product  p ⊗ q.
    Both vectors use [q1, q2, q3, q4] convention (q4 = scalar).
    """
    p1, p2, p3, p4 = p
    q1, q2, q3, q4 = q
    return np.array([
         p4*q1 + p3*q2 - p2*q3 + p1*q4,
        -p3*q1 + p4*q2 + p1*q3 + p2*q4,
         p2*q1 - p1*q2 + p4*q3 + p3*q4,
        -p1*q1 - p2*q2 - p3*q3 + p4*q4
    ])


def skew(v: np.ndarray) -> np.ndarray:
    """3×3 skew-symmetric matrix so that skew(v) @ w = v × w."""
    return np.array([
        [ 0,    -v[2],  v[1]],
        [ v[2],  0,    -v[0]],
        [-v[1],  v[0],  0   ]
    ])


def delta_theta_to_quat(dtheta: np.ndarray) -> np.ndarray:
    """
    Small rotation vector δθ  →  unit quaternion [vec; scalar].
    Uses exact (non-small-angle) formula for robustness.
    """
    angle = np.linalg.norm(dtheta)
    if angle < 1e-10:
        return np.array([0.0, 0.0, 0.0, 1.0])
    axis = dtheta / angle
    s = np.sin(angle / 2.0)
    return np.concatenate([s * axis, [np.cos(angle / 2.0)]])


def normalise_quat(q: np.ndarray) -> np.ndarray:
    return q / np.linalg.norm(q)


# ─────────────────────────────────────────────────────────────────────────────
# Measurement model
# ─────────────────────────────────────────────────────────────────────────────

def h(x: np.ndarray) -> np.ndarray:
    """
    Nonlinear measurement function  z = h(x).

    Algorithm
    ---------
    1. Form LOS vector in inertial frame:  ρ^I = r^T − r^C
    2. Rotate to chaser body frame:        ρ^C = C(q^{I→C}) ρ^I
    3. Extract azimuth and elevation from  ρ^C

    Returns
    -------
    z : (2,)  [azimuth, elevation]  in radians
    """
    rho_I = x[IDX_rT] - x[IDX_rC]              # LOS in inertial   (3,)
    C_IC  = quat2dcm(x[IDX_qC])                 # inertial → chaser body
    rho_C = C_IC @ rho_I                         # LOS in body frame (3,)

    az = np.arctan2(rho_C[1], rho_C[0])
    el = np.arctan2(rho_C[2], np.sqrt(rho_C[0]**2 + rho_C[1]**2))

    return np.array([az, el])


def H_analytic(x: np.ndarray) -> np.ndarray:
    """
    Analytic measurement Jacobian  H = ∂z/∂(δx)  (2 × 30).

    Derivation
    ----------
    z depends on ρ^C = C_{IC} (r^T − r^C).

    Partials w.r.t. each error-state block:

      ∂z/∂δr^T  =  (∂z/∂ρ^C) C_{IC}          [LOS position changes r^T]
      ∂z/∂δr^C  = −(∂z/∂ρ^C) C_{IC}          [opposite sign]
      ∂z/∂δθ^C  =  (∂z/∂ρ^C)(−[ρ^C]×)        [attitude perturbation rotates ρ^C]

    Attitude perturbation model (right-multiplicative):
      q_new = q_nom ⊗ δq(δθ)  →  C_new = C_nom (I − [δθ]×)
      ρ^C_new = C_nom (I − [δθ]×) ρ^I = ρ^C + C_nom [ρ^I]× δθ
             equivalently: ρ^C + [ρ^C]× C_nom δθ
      ⟹  ∂ρ^C/∂δθ^C = [ρ^C]× C_nom

    All other error-state blocks (δv, δω, δb, δS, δO) have zero partial.

    Returns
    -------
    H : (2, 30)
    """
    rho_I = x[IDX_rT] - x[IDX_rC]
    C_IC  = quat2dcm(x[IDX_qC])
    rho_C = C_IC @ rho_I

    rx, ry, rz = rho_C
    d2  = rx**2 + ry**2                    # horizontal range squared
    d   = np.sqrt(d2)                      # horizontal range
    rn2 = rx**2 + ry**2 + rz**2           # total range squared

    # Guard against singularities (target directly above/below)
    if d < 1e-10:
        raise ValueError(
            f"Azimuth singularity: LOS is nearly vertical (d = {d:.2e}). "
            "Consider switching to a different angle parametrisation near the pole."
        )

    # ∂[az, el] / ∂ρ^C   →   (2 × 3)
    #   az = atan2(ry, rx)   →   ∂az/∂ρ = [-ry/d², rx/d², 0]
    #   el = atan2(rz, d)    →   ∂el/∂ρ = [-rx*rz/(rn²*d), -ry*rz/(rn²*d), d/rn²]
    dz_drho = np.array([
        [-ry / d2,           rx / d2,           0.0      ],
        [-rx*rz / (rn2 * d), -ry*rz / (rn2 * d), d / rn2]
    ])                                           # (2×3)

    dz_drT  =  dz_drho @ C_IC              # (2×3)
    dz_drC  = -dz_drT                      # (2×3)
    # Perturbation model (right-multiply):  q_new = q_nom ⊗ δq(δθ)
    #   C_new  = C_nom · C(δq) ≈ C_nom (I − [δθ]×)
    #   ρ^C_new = ρ^C + C_nom [ρ^I]× δθ  =  ρ^C + [ρ^C]× C_nom δθ
    #   ∂ρ^C/∂δθ = [ρ^C]× C_nom
    dz_dthC = dz_drho @ (skew(rho_C) @ C_IC)   # (2×3)

    H = np.zeros((2, N_ERR))
    H[:, IDX_E_rT]  = dz_drT
    H[:, IDX_E_rC]  = dz_drC
    H[:, IDX_E_thC] = dz_dthC
    # δvT, δθT, δωT, δvC, δbw, δSs, δOo → all zero

    return H


def H_numerical(x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """
    Numerical (central-difference) Jacobian for validation against H_analytic.

    Perturbations are applied in error-state space:
      - position/velocity/angular-velocity/parameter blocks: additive
      - attitude blocks (δθ^T, δθ^C): multiplicative via delta_theta_to_quat
    """
    H = np.zeros((2, N_ERR))
    z0 = h(x)

    for i in range(N_ERR):
        # ── forward perturbation
        dxp = np.zeros(N_ERR); dxp[i] = +eps
        xp  = apply_correction(x, dxp)
        zp  = h(xp)

        # ── backward perturbation
        dxm = np.zeros(N_ERR); dxm[i] = -eps
        xm  = apply_correction(x, dxm)
        zm  = h(xm)

        H[:, i] = (zp - zm) / (2 * eps)

    return H


# ─────────────────────────────────────────────────────────────────────────────
# EKF measurement update
# ─────────────────────────────────────────────────────────────────────────────

def apply_correction(x: np.ndarray, delta_x: np.ndarray) -> np.ndarray:
    """
    Apply error-state correction δx to nominal state x.

    Additive for all states except quaternions, which use multiplicative
    (right-multiply by δq derived from rotation vector δθ).
    """
    xp = x.copy()

    # ── additive corrections
    xp[IDX_rT] += delta_x[IDX_E_rT]
    xp[IDX_vT] += delta_x[IDX_E_vT]
    xp[IDX_wT] += delta_x[IDX_E_wT]
    xp[IDX_rC] += delta_x[IDX_E_rC]
    xp[IDX_vC] += delta_x[IDX_E_vC]
    xp[IDX_bw] += delta_x[IDX_E_bw]
    xp[IDX_Ss] += delta_x[IDX_E_Ss]
    xp[IDX_Oo] += delta_x[IDX_E_Oo]

    # ── multiplicative attitude corrections  q_new = q_nom ⊗ δq(δθ)
    xp[IDX_qT] = normalise_quat(
        quat_mult(x[IDX_qT], delta_theta_to_quat(delta_x[IDX_E_thT]))
    )
    xp[IDX_qC] = normalise_quat(
        quat_mult(x[IDX_qC], delta_theta_to_quat(delta_x[IDX_E_thC]))
    )

    return xp


def ekf_measurement_update(
    x:      np.ndarray,
    P:      np.ndarray,
    z_meas: np.ndarray,
    R_meas: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    EKF measurement update (Joseph form for numerical stability).

    Parameters
    ----------
    x      : (32,)  nominal state
    P      : (30,30) error-state covariance
    z_meas : (2,)   [azimuth, elevation] measurement in radians
    R_meas : (2,2)  measurement noise covariance

    Returns
    -------
    x_plus : (32,)  updated state
    P_plus : (30,30) updated covariance  (Joseph form)
    innov  : (2,)   innovation  z_meas − h(x)
    """
    z_pred = h(x)
    innov  = z_meas - z_pred

    H = H_analytic(x)                         # (2×30)
    S = H @ P @ H.T + R_meas                  # (2×2)  innovation covariance
    K = P @ H.T @ np.linalg.solve(S.T, np.eye(2)).T  # (30×2)  Kalman gain

    delta_x = K @ innov                        # (30,)  error-state correction
    x_plus  = apply_correction(x, delta_x)

    # Joseph form:  P⁺ = (I − KH) P (I − KH)ᵀ + K R Kᵀ
    IKH    = np.eye(N_ERR) - K @ H
    P_plus = IKH @ P @ IKH.T + K @ R_meas @ K.T

    return x_plus, P_plus, innov


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

def _make_test_state() -> np.ndarray:
    """Construct a realistic, non-trivial 32-element state vector."""
    x = np.zeros(N_STATE)

    # Target: ~100 m ahead in x, slight cross-track offset
    x[IDX_rT] = np.array([100.0, 5.0, -3.0])
    x[IDX_vT] = np.array([-0.1,  0.01, 0.0])

    # Target attitude: 30° rotation about z-axis
    ang = np.deg2rad(30)
    x[IDX_qT] = np.array([0.0, 0.0, np.sin(ang/2), np.cos(ang/2)])
    x[IDX_wT] = np.array([0.01, -0.005, 0.002])

    # Chaser: at origin, slight velocity
    x[IDX_rC] = np.array([0.0, 0.0, 0.0])
    x[IDX_vC] = np.array([0.0, 0.0, 0.0])

    # Chaser attitude: 15° rotation about y-axis (camera not perfectly aligned)
    ang2 = np.deg2rad(15)
    x[IDX_qC] = np.array([0.0, np.sin(ang2/2), 0.0, np.cos(ang2/2)])

    # Parameter states near zero (small biases)
    x[IDX_bw] = np.array([1e-4, -2e-4,  1e-5])
    x[IDX_Ss] = np.array([1.0,   1.0,   1.0 ])   # unit scale factors
    x[IDX_Oo] = np.array([0.0,   0.0,   0.0 ])

    return x


def test_measurement_function():
    """Verify h(x) returns valid azimuth/elevation angles."""
    x = _make_test_state()
    z = h(x)

    assert z.shape == (2,), "z must be (2,)"
    assert -np.pi <= z[0] <= np.pi,  f"Azimuth out of range: {np.rad2deg(z[0]):.2f}°"
    assert -np.pi/2 <= z[1] <= np.pi/2, f"Elevation out of range: {np.rad2deg(z[1]):.2f}°"

    rho_I = x[IDX_rT] - x[IDX_rC]
    C_IC  = quat2dcm(x[IDX_qC])
    rho_C = C_IC @ rho_I

    az_expected = np.arctan2(rho_C[1], rho_C[0])
    el_expected = np.arctan2(rho_C[2], np.sqrt(rho_C[0]**2 + rho_C[1]**2))
    assert np.allclose(z, [az_expected, el_expected], atol=1e-12)

    print("✓ test_measurement_function passed")
    print(f"  Azimuth  = {np.rad2deg(z[0]):.4f}°")
    print(f"  Elevation= {np.rad2deg(z[1]):.4f}°")


def test_jacobian_vs_numerical():
    """Verify analytic H matches numerical finite-difference Jacobian."""
    x = _make_test_state()

    H_a = H_analytic(x)
    H_n = H_numerical(x, eps=1e-6)

    max_err = np.max(np.abs(H_a - H_n))
    assert max_err < 1e-7, f"Jacobian mismatch: max error = {max_err:.2e}"

    print("✓ test_jacobian_vs_numerical passed")
    print(f"  Max |H_analytic − H_numerical| = {max_err:.2e}")


def test_jacobian_sparsity():
    """Confirm that only position and chaser attitude columns are nonzero."""
    x   = _make_test_state()
    H   = H_analytic(x)
    tol = 1e-14

    nonzero_cols = set(np.where(np.any(np.abs(H) > tol, axis=0))[0])
    expected_nonzero = (
        set(range(*IDX_E_rT.indices(N_ERR)))    # δr^T
        | set(range(*IDX_E_rC.indices(N_ERR)))  # δr^C
        | set(range(*IDX_E_thC.indices(N_ERR))) # δθ^C
    )
    unexpected = nonzero_cols - expected_nonzero
    assert not unexpected, f"Unexpected nonzero columns: {unexpected}"

    print("✓ test_jacobian_sparsity passed")
    print(f"  Nonzero columns: {sorted(nonzero_cols)} "
          f"→ δr^T[0:3], δr^C[12:15], δθ^C[18:21]")


def test_measurement_update_reduces_uncertainty():
    """After an EKF update, covariance trace must decrease."""
    x = _make_test_state()
    P = np.eye(N_ERR) * 1.0        # large initial uncertainty

    sigma_angle_rad = np.deg2rad(0.1)          # 0.1° 1-sigma sensor noise
    R = (sigma_angle_rad**2) * np.eye(2)

    z_true  = h(x)
    z_noisy = z_true + np.random.default_rng(42).normal(0, sigma_angle_rad, 2)

    x_plus, P_plus, innov = ekf_measurement_update(x, P, z_noisy, R)

    assert np.trace(P_plus) < np.trace(P), "Covariance did not decrease after update"

    print("✓ test_measurement_update_reduces_uncertainty passed")
    print(f"  Trace(P)  before = {np.trace(P):.4f}")
    print(f"  Trace(P)  after  = {np.trace(P_plus):.4f}")
    print(f"  Innovation       = {np.rad2deg(innov)} deg")


def test_zero_innovation_no_state_change():
    """
    If measured angles exactly match predicted, the state should not change
    (within floating-point tolerance).
    """
    x = _make_test_state()
    P = np.eye(N_ERR) * 0.01
    R = (np.deg2rad(0.1)**2) * np.eye(2)

    z_perfect = h(x)
    x_plus, P_plus, innov = ekf_measurement_update(x, P, z_perfect, R)

    assert np.allclose(innov, 0.0, atol=1e-14), "Non-zero innovation on perfect measurement"
    assert np.allclose(x_plus[IDX_rT], x[IDX_rT], atol=1e-12), "Position changed on zero innovation"
    assert np.allclose(x_plus[IDX_rC], x[IDX_rC], atol=1e-12), "Chaser position changed"

    print("✓ test_zero_innovation_no_state_change passed")


if __name__ == "__main__":
    np.random.seed(0)
    print("=" * 60)
    print("Woffinden & Geller — Angles-Only EKF Measurement Model")
    print("=" * 60)
    test_measurement_function()
    test_jacobian_vs_numerical()
    test_jacobian_sparsity()
    test_measurement_update_reduces_uncertainty()
    test_zero_innovation_no_state_change()
    print("\nAll tests passed ✓")