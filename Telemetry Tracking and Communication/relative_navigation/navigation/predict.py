import numpy as np
import config
from scipy.integrate import solve_ivp
from dynamics.dynamics import *
from dynamics.attitude import *

def ekf_predict(x_est, P, omega_gyro, u_applied, tau_applied, I_r, I_c, noise, param_noise, dt):

    # -- 1. Propagate nominal state through dynamics --------
    x_prop = x_est.copy()
    x_prop = np.concatenate([x_prop, omega_gyro])  # replace chaser omega with gyro measurement

    sol = solve_ivp(
        dynamics_predict,
        [0, config.SIM_DT], x_prop, method='RK45', 
        args=(u_applied, tau_applied, I_r, I_c, param_noise, ),
        rtol=config.tol, atol=config.tol
    )
    x_new    = sol.y[0:23 , -1] # remove chaser omega from state
    x_params = sol.y[26:35 , -1]
    x_new = np.concatenate([x_new, x_params])
    x_new[6:10]  = normalize_quat(x_new[6:10])
    x_new[19:23] = normalize_quat(x_new[19:23])
    m_minus = np.concatenate([x_new[0:6], quat2rot_vet(x_new[6:10]), x_new[10:19], quat2rot_vet(x_new[19:23]), x_new[23:32]])

    # -- 2. Compute Jacobian ---------------------------------
    tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o = param_noise
    F = predict_jacobian(m_minus, omega_gyro, I_r, I_c, x_new[6:10], tau_b, tau_s, tau_o)

    # -- 3. Propagate covariance ------------------------------
    # Discrete: P = F P F^T + Q
    sigma_v_T, sigma_w_T, sigma_v_C, sigma_w_C = noise
    Q = compute_Q(dt, sigma_v_T, sigma_w_T, sigma_v_C, tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o)
    P_new = F @ P @ F.T + config.Q

    return m_minus, P_new


def compute_Q(dt, sigma_v_T, sigma_w_T, sigma_v_C, tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o):

    I = np.eye(3)

    # Q_wT — target process noise (12×12)
    Q_wT = np.zeros((12, 12))
    Q_wT[0:6, 0:6] = sigma_v_T**2 * np.block([
        [I * (dt**3/3), I * (dt**2/2)],
        [I * (dt**2/2), I *  dt      ],
    ])
    Q_wT[6:12, 6:12] = sigma_w_T**2 * np.block([
        [I * (dt**3/3), I * (dt**2/2)],
        [I * (dt**2/2), I *  dt      ],
    ])

    # Q_wC — chaser process noise (9×9)
    Q_wC = np.zeros((9, 9))
    Q_wC[0:6, 0:6] = sigma_v_C**2 * np.block([
        [I * (dt**3/3), I * (dt**2/2)],
        [I * (dt**2/2), I *  dt      ],
    ])
    # attitude block is zero — model replacement mode, gyro drives attitude directly

    # Q_wP — parameter process noise (9×9), first-order Markov closed form
    Q_wP = np.block([
        [sigma_b**2 * (1 - np.exp(-2*dt/tau_b)) * I, np.zeros((3,3)),                                np.zeros((3,3))                               ],
        [np.zeros((3,3)),                              sigma_s**2 * (1 - np.exp(-2*dt/tau_s)) * I,   np.zeros((3,3))                               ],
        [np.zeros((3,3)),                              np.zeros((3,3)),                               sigma_o**2 * (1 - np.exp(-2*dt/tau_o)) * I   ],
    ])

    # Assemble full Q (30×30) — eq A6
    Q = np.zeros((30, 30))
    Q[0:12, 0:12]  = Q_wT
    Q[12:21, 12:21] = Q_wC
    Q[21:30, 21:30] = Q_wP

    return Q