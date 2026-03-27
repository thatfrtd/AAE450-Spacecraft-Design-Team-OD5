"""
EKF Predict Step
================
Noise enters the filter ONLY through the Q matrix — never through the
dynamics integrator.  The Q matrix models two distinct uncertainties:

  1. Unmodelled forces on target/chaser (drag, SRP, thruster imprecision)
     → position/velocity blocks of Q
     → tuned via sigma_vel_T, sigma_vel_C in config

  2. Parameter drift (gyro bias, scale factor, misalignment)
     → parameter block of Q via Gauss-Markov closed form
     → tuned via sigma_b, sigma_s, sigma_o in config
"""

import numpy as np
import config
from scipy.integrate import solve_ivp
from dynamics.dynamics import *
from dynamics.attitude import *


def ekf_predict(x_est, P, omega_gyro, u_applied, tau_applied, I_r, I_c,
                noise, param_noise, dt):
    """
    EKF predict step — deterministic state propagation + covariance update via Q.

    x_est layout (32 elements):
      [0:13]  target  (r, v, q, w)
      [13:23] chaser  (r, v, q)        ← no omega in estimated state
      [23:32] params  (b_w, S_s, O_o)

    dynamics_predict expects (35 elements):
      [0:13]  target  (r, v, q, w)
      [13:26] chaser  (r, v, q, w)    ← omega inserted from gyro
      [26:35] params  (b_w, S_s, O_o)
    """

    # ── Unpack Q tuning parameters ────────────────────────────────────────────
    sigma_v_T = noise[0:3]     # target velocity process noise [km/s]
    sigma_w_T = noise[3:6]     # target angular velocity process noise [rad/s]
    sigma_v_C = noise[6:9]     # chaser velocity process noise [km/s]
    # sigma_w_C = noise[9:12]  # (unused — attitude driven by gyro)

    tau_b   = param_noise[0:3];   sigma_b = param_noise[9:12]
    tau_s   = param_noise[3:6];   sigma_s = param_noise[12:15]
    tau_o   = param_noise[6:9];   sigma_o = param_noise[15:18]

    # ── 1. Insert gyro omega into propagation state ───────────────────────────
    if omega_gyro is None:
        omega_gyro = np.zeros(3)

    x_prop = np.concatenate([
        x_est[0:23],    # target(13) + chaser r,v,q(10)
        omega_gyro,     # chaser omega from gyro (3)
        x_est[23:32],   # params (9)
    ])                  # total = 35

    # ── 2. Propagate nominal state — fully deterministic ─────────────────────
    sol = solve_ivp(
        dynamics_predict,
        [0, dt], x_prop, method='RK45',
        args=(u_applied, tau_applied, I_r, I_c, param_noise,),
        rtol=config.tol, atol=config.tol,
    )

    x_full = sol.y[:, -1]   # (35,)

    # Drop chaser omega (gyro-driven, not in estimated state)
    x_new = np.concatenate([
        x_full[0:23],   # target(13) + chaser r,v,q(10)
        x_full[26:35],  # params(9)
    ])                  # back to 32

    x_new[6:10]  = normalize_quat(x_new[6:10])
    x_new[19:23] = normalize_quat(x_new[19:23])

    # ── 3. Build 30-element error-state vector for Jacobian ───────────────────
    m_minus = np.concatenate([
        x_new[0:6],
        np.zeros(3),    # δθ_T
        x_new[10:19],
        np.zeros(3),   # δθ_C
        x_new[23:32],
    ])

    # ── 4. Linearise and discretise ───────────────────────────────────────────
    F_c = predict_jacobian(m_minus, omega_gyro, I_r, I_c,
                           x_new[6:10], tau_b, tau_s, tau_o)

    # First-order hold: Φ ≈ I + F·dt  (adequate for dt = 1s orbital dynamics)
    Phi = np.eye(30) + F_c * dt

    # ── 5. Build process noise covariance Q ───────────────────────────────────
    Q = compute_Q(dt, sigma_v_T, sigma_w_T, sigma_v_C,
                  tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o)

    # ── 6. Propagate covariance ───────────────────────────────────────────────
    P_new = Phi @ P @ Phi.T + Q

    # Guard against numerical instability
    if np.any(np.isnan(P_new)) or np.any(np.isinf(P_new)):
        print("WARNING: P_new contains NaN/Inf — reverting to prior P")
        P_new = P.copy()

    return x_new, P_new


def compute_Q(dt, sigma_v_T, sigma_w_T, sigma_v_C,
              tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o):
    """
    Discrete-time process noise covariance Q (30×30).

    Derivation
    ----------
    For a continuous-time random walk on velocity with spectral density q:
      Q_rv = q * [[dt³/3,  dt²/2],
                  [dt²/2,  dt   ]]
    where q = sigma_v² / dt  (converting 1-sigma per step → spectral density).

    This gives the standard "continuous white noise acceleration" model where
    sigma_v is the 1-sigma velocity uncertainty injected PER TIMESTEP.

    Parameter noise uses the Gauss-Markov steady-state variance:
      Q_param = sigma² * (1 - exp(-2*dt/tau))
    """
    I  = np.eye(3)
    Z3 = np.zeros((3, 3))

    def _vel_block(sig):
        """
        6×6 position-velocity Q block for a velocity random walk.
        sig : 1-sigma velocity noise per step [km/s or rad/s]
        """
        q = np.mean(sig**2) / dt    # spectral density [km²/s³]
        return q * np.block([
            [I * (dt**3 / 3), I * (dt**2 / 2)],
            [I * (dt**2 / 2), I * dt          ],
        ])

    # ── Target (12×12): position/velocity + attitude/angular-velocity ─────────
    Q_T = np.zeros((12, 12))
    Q_T[0:6,  0:6]  = _vel_block(sigma_v_T)   # position + velocity
    Q_T[6:12, 6:12] = _vel_block(sigma_w_T)   # attitude + angular velocity

    # ── Chaser (9×9): position/velocity only (attitude driven by gyro) ────────
    Q_C = np.zeros((9, 9))
    Q_C[0:6, 0:6] = _vel_block(sigma_v_C)

    # ── Parameters (9×9): Gauss-Markov closed form ────────────────────────────
    def _markov(sig, tau):
        s2   = np.mean(sig**2)
        tau0 = np.mean(tau)
        if tau0 > 1e6:               # frozen (tau → ∞)
            return Z3
        return s2 * (1 - np.exp(-2 * dt / tau0)) * I

    Q_P = np.block([
        [_markov(sigma_b, tau_b), Z3,                    Z3                   ],
        [Z3,                       _markov(sigma_s, tau_s), Z3                ],
        [Z3,                       Z3,                    _markov(sigma_o, tau_o)],
    ])

    # ── Assemble (30×30) ──────────────────────────────────────────────────────
    Q = np.zeros((30, 30))
    Q[0:12,  0:12]  = Q_T
    Q[12:21, 12:21] = Q_C
    Q[21:30, 21:30] = Q_P

    return Q