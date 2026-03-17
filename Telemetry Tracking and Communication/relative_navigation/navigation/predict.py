import numpy as np
import config
from scipy.integrate import solve_ivp
from dynamics.dynamics import *
from dynamics.attitude import *


def ekf_predict(x_est, P, omega_gyro, u_applied, tau_applied, I_r, I_c, noise, param_noise, dt):
    """
    EKF predict step.

    x_est layout (32 elements):
      [0:13]  target  (r, v, q, w)
      [13:23] chaser  (r, v, q)       ← no omega in estimated state
      [23:32] params  (b_w, S_s, O_o)

    dynamics_predict expects (35 elements):
      [0:13]  target  (r, v, q, w)
      [13:26] chaser  (r, v, q, w)    ← omega inserted from gyro
      [26:35] params  (b_w, S_s, O_o)
    """

    # -- Unpack noise by slice (Bug fix: was tuple-unpacking 18-element arrays)
    tau_b   = param_noise[0:3]
    tau_s   = param_noise[3:6]
    tau_o   = param_noise[6:9]
    sigma_b = param_noise[9:12]
    sigma_s = param_noise[12:15]
    sigma_o = param_noise[15:18]

    sigma_v_T = noise[0:3]
    sigma_w_T = noise[3:6]
    sigma_v_C = noise[6:9]
    sigma_w_C = noise[9:12]

    # -- 1. Build 35-element propagation state --------------------------------
    # Bug fix: was np.concatenate([x_prop, omega_gyro]) which appended omega
    # to the end giving wrong layout. Must insert omega into chaser block.
    x_prop = np.concatenate([
        x_est[0:23],    # target (13) + chaser r,v,q (10)
        omega_gyro,     # chaser omega (3) — inserted here → chaser now 13 elem
        x_est[23:32],   # params (9)
    ])                  # total = 35

    # -- 2. Propagate nominal state -------------------------------------------
    sol = solve_ivp(
        dynamics_predict,
        [0, dt], x_prop, method='RK45',
        args=(u_applied, tau_applied, I_r, I_c, param_noise,),
        rtol=config.tol, atol=config.tol,
    )

    x_full = sol.y[:, -1]                  # (35,)

    # Bug fix: extract chaser WITHOUT omega (drop [23:26]) and rejoin params
    x_new = np.concatenate([
        x_full[0:23],   # target (13) + chaser r,v,q (10)
        x_full[26:35],  # params (9)  — skip chaser omega at [23:26]
    ])                  # back to 32

    x_new[6:10]  = normalize_quat(x_new[6:10])
    x_new[19:23] = normalize_quat(x_new[19:23])

    # -- 3. Build error-state vector for Jacobian (30-element) ----------------
    m_minus = np.concatenate([
        x_new[0:6],                      # δr_T, δv_T
        quat2rot_vet(x_new[6:10]),        # δθ_T  (3)
        x_new[10:19],                     # δω_T, δr_C, δv_C
        quat2rot_vet(x_new[19:23]),       # δθ_C  (3)
        x_new[23:32],                     # params (9)
    ])                                    # total = 30

    # -- 4. Compute continuous Jacobian and discretise ------------------------
    F_c = predict_jacobian(m_minus, omega_gyro, I_r, I_c, x_new[6:10],
                           tau_b, tau_s, tau_o)

    # State transition matrix: first-order hold  Φ ≈ I + F*dt
    # expm() is numerically unstable for this 30×30 orbital-mechanics matrix
    Phi = np.eye(30) + F_c * dt

    # -- 5. Compute process noise covariance ----------------------------------
    Q = compute_Q(dt, sigma_v_T, sigma_w_T, sigma_v_C,
                  tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o)

    # -- 6. Propagate covariance ----------------------------------------------
    P_new = Phi @ P @ Phi.T + Q

    # -- 7. NaN / Inf guard ---------------------------------------------------
    if np.any(np.isnan(P_new)) or np.any(np.isinf(P_new)):
        print("WARNING: P_new contains NaN/Inf — reverting to prior P")
        P_new = P.copy()
    if np.any(np.isnan(x_new)):
        print("WARNING: x_new contains NaN — check dynamics")

    return x_new, P_new


def compute_Q(dt, sigma_v_T, sigma_w_T, sigma_v_C,
              tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o):
    """
    Discrete process noise covariance Q (30×30).

    All sigma/tau inputs may be 3-element arrays (one per axis).
    """
    I = np.eye(3)

    # Helper: build 6×6 continuous white-noise position/velocity block
    def _pv_block(sig):
        s2 = np.mean(sig**2)          # scalar variance (isotropic approx)
        return s2 * np.block([
            [I * (dt**3 / 3), I * (dt**2 / 2)],
            [I * (dt**2 / 2), I *  dt         ],
        ])

    # Q_wT — target process noise (12×12)
    Q_wT = np.zeros((12, 12))
    Q_wT[0:6,  0:6]  = _pv_block(sigma_v_T)
    Q_wT[6:12, 6:12] = _pv_block(sigma_w_T)

    # Q_wC — chaser process noise (9×9)
    Q_wC = np.zeros((9, 9))
    Q_wC[0:6, 0:6] = _pv_block(sigma_v_C)
    # attitude block is zero — gyro drives attitude directly in predict step

    # Q_wP — parameter process noise (9×9), first-order Markov closed form
    def _markov(sig, tau):
        s2   = np.mean(sig**2)
        tau0 = np.mean(tau)
        if tau0 < 1e6:                 # non-frozen
            return s2 * (1 - np.exp(-2 * dt / tau0)) * I
        else:                          # frozen (tau → ∞) → zero noise
            return np.zeros((3, 3))

    Q_wP = np.block([
        [_markov(sigma_b, tau_b), np.zeros((3,3)),          np.zeros((3,3))         ],
        [np.zeros((3,3)),          _markov(sigma_s, tau_s),  np.zeros((3,3))         ],
        [np.zeros((3,3)),          np.zeros((3,3)),           _markov(sigma_o, tau_o)],
    ])

    # Assemble full Q (30×30)
    Q = np.zeros((30, 30))
    Q[0:12,  0:12]  = Q_wT
    Q[12:21, 12:21] = Q_wC
    Q[21:30, 21:30] = Q_wP

    return Q