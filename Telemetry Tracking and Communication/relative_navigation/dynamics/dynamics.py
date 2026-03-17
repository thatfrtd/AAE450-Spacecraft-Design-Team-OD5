import numpy as np
from utils.constants import *
from dynamics.attitude import *


def target_dynamics(t, y, I, noise_g, noise_omega):
    """Two-body + J2 dynamics for the target. noise_g/omega are callables f(t)->(3,)."""
    r = y[0:3];  v = y[3:6];  q = y[6:10];  w = y[10:13]
    rnorm = np.linalg.norm(r)

    dist_g     = noise_g(t)
    dist_omega = noise_omega(t)

    a = -MU_EARTH * r / rnorm**3 + j2_acceleration(r, MU_EARTH) + dist_g

    q_dot = kde(q, w)
    w_dot = dde(I, w, q, r, dist_omega, np.zeros(3))

    dydt = np.zeros_like(y)
    dydt[0:3]  = v
    dydt[3:6]  = a
    dydt[6:10] = q_dot
    dydt[10:13] = w_dot
    return dydt


def chaser_dynamics(t, y, I, noise_g, noise_omega, u_dv, u_tau):
    """
    Two-body + J2 dynamics for the chaser.

    u_dv  : (3,) translational command in CHASER BODY frame [km/s²]  ← Eq 52
    u_tau : (3,) torque command [kg·km²/s²]

    u_dv is rotated to inertial before adding to acceleration.
    """
    r = y[0:3];  v = y[3:6];  q = y[6:10];  w = y[10:13]
    rnorm = np.linalg.norm(r)

    dist_g     = noise_g(t)
    dist_omega = noise_omega(t)

    # Rotate body-frame thrust to inertial frame
    T_CI = quat2dcm(quat_inv(q))      # chaser body → inertial  (C^T)
    u_dv_inertial = T_CI @ u_dv

    a = -MU_EARTH * r / rnorm**3 + j2_acceleration(r, MU_EARTH) + dist_g + u_dv_inertial

    q_dot = kde(q, w)
    w_dot = dde(I, w, q, r, dist_omega, u_tau)

    dydt = np.zeros_like(y)
    dydt[0:3]  = v
    dydt[3:6]  = a
    dydt[6:10] = q_dot
    dydt[10:13] = w_dot
    return dydt


def dynamics_truth(t, y, u_applied, tau_applied, I_r, I_c, noise, param_noise):
    """Full 35-element truth dynamics: target(13) + chaser(13) + params(9)."""
    dydt = np.zeros_like(y)

    sigma_v_T     = noise[0:3]
    sigma_omega_T = noise[3:6]
    sigma_v_C     = noise[6:9]
    sigma_omega_C = noise[9:12]

    tau_b   = param_noise[0:3];   tau_s   = param_noise[3:6];   tau_o   = param_noise[6:9]
    sigma_b = param_noise[9:12];  sigma_s = param_noise[12:15]; sigma_o = param_noise[15:18]

    noise_g_T     = lambda t: np.random.normal(0, sigma_v_T)
    noise_omega_T = lambda t: np.random.normal(0, sigma_omega_T)
    noise_g_C     = lambda t: np.zeros(3)
    noise_omega_C = lambda t: np.zeros(3)

    dydt[0:13]  = target_dynamics(t, y[0:13],  I_r, noise_g_T,     noise_omega_T)
    dydt[13:26] = chaser_dynamics(t, y[13:26], I_c, noise_g_C,     noise_omega_C,
                                  u_applied, tau_applied)
    dydt[26:35] = parameter_dynamics(t, y[26:35],
                                     tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o)
    return dydt


def dynamics_predict(t, y, u_applied, tau_applied, I_r, I_c, param_noise):
    """Full 35-element predict dynamics (no process noise on translational states)."""
    dydt = np.zeros_like(y)

    tau_b   = param_noise[0:3];   tau_s   = param_noise[3:6];   tau_o   = param_noise[6:9]
    sigma_b = param_noise[9:12];  sigma_s = param_noise[12:15]; sigma_o = param_noise[15:18]

    noise_zero = lambda t: np.zeros(3)

    dydt[0:13]  = target_dynamics(t, y[0:13],  I_r, noise_zero, noise_zero)
    dydt[13:26] = chaser_dynamics(t, y[13:26], I_c, noise_zero, noise_zero,
                                  u_applied, tau_applied)
    dydt[26:35] = parameter_dynamics(t, y[26:35],
                                     tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o)
    return dydt


def parameter_dynamics(t, y, tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o):
    """First-order Gauss-Markov bias/misalignment model."""
    b_w = y[0:3];  eps_s = y[3:6];  eps_o = y[6:9]

    def _drive(sigma):
        return np.zeros(3) if np.all(sigma == 0) else np.random.normal(0, sigma, 3)

    return np.concatenate([
        -b_w   / tau_b + _drive(sigma_b),
        -eps_s / tau_s + _drive(sigma_s),
        -eps_o / tau_o + _drive(sigma_o),
    ])


def predict_jacobian(m_minus, w_C, I_t, I_c, q_T, tau_b, tau_s, tau_o):
    """Linearised dynamics Jacobian F = ∂f/∂δx  (30×30 error-state)."""
    r_T = m_minus[0:3];  w_T = m_minus[9:12]
    r_C = m_minus[12:15];  b_w = m_minus[21:24]

    T_I_T  = quat2dcm(q_T)
    I_inv  = np.linalg.inv(I_t)
    r_norm = np.linalg.norm(r_T)
    r_hat_B = T_I_T @ r_T
    skew_rB = skew(r_hat_B)

    # Target block (12×12)
    F_T2 = g_t_Jacobian(r_T)
    F_T3 = -skew(w_T)
    term1 = (3*MU_EARTH/r_norm**5) * I_inv @ (
        skew_rB @ I_t @ T_I_T - skew(I_t @ r_hat_B) @ T_I_T)
    term2 = (-15*MU_EARTH/r_norm**7) * I_inv @ (
        skew_rB @ I_t @ r_hat_B[:,None] @ r_T[None,:] @ T_I_T.T)
    F_T5 = term1 + term2
    F_T6 = (3*MU_EARTH/r_norm**5) * I_inv @ (
        skew_rB @ I_t @ skew_rB - skew(I_t @ r_hat_B) @ skew_rB)
    F_T7 = np.linalg.solve(I_t, skew(I_t @ w_T) - skew(w_T) @ I_t)

    F_TT = np.block([
        [np.zeros((3,3)), np.eye(3),        np.zeros((3,3)), np.zeros((3,3))],
        [F_T2,            np.zeros((3,3)),  np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)),  F_T3,            np.eye(3)      ],
        [F_T5,            np.zeros((3,3)),  F_T6,            F_T7           ],
    ])

    # Chaser block (9×9)
    F_CC = np.block([
        [np.zeros((3,3)), np.eye(3),            np.zeros((3,3))],
        [g_t_Jacobian(r_C), np.zeros((3,3)),   np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)),       -skew(w_C - b_w)],
    ])

    # Chaser-parameter coupling (9×9)
    F_CP = np.block([
        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [-np.eye(3),      np.zeros((3,3)), np.zeros((3,3))],
    ])

    # Parameter block (9×9)
    F_PP = np.block([
        [(-1/tau_b[0])*np.eye(3), np.zeros((3,3)),       np.zeros((3,3))      ],
        [np.zeros((3,3)),          (-1/tau_s[0])*np.eye(3), np.zeros((3,3))   ],
        [np.zeros((3,3)),          np.zeros((3,3)),        (-1/tau_o[0])*np.eye(3)],
    ])

    F = np.zeros((30, 30))
    F[0:12,  0:12]  = F_TT
    F[12:21, 12:21] = F_CC
    F[12:21, 21:30] = F_CP
    F[21:30, 21:30] = F_PP
    return F


def g_t_Jacobian(r):
    """∂a/∂r for two-body + J2."""
    rnorm = np.linalg.norm(r);  z = r[2]
    I3 = np.eye(3);  rrT = np.outer(r, r)
    dadr_tb = MU_EARTH * (3*rrT/rnorm**5 - I3/rnorm**3)
    n3 = np.array([0., 0., 1.])
    n3n3T = np.outer(n3, n3);  rn3T = np.outer(r, n3);  n3rT = np.outer(n3, r)
    dadr_j2 = (-3*1.0826e-3*MU_EARTH*R_EARTH**2/2) * (
        I3/rnorm**5 + 2*n3n3T/rnorm**5
        - 5*(rrT + z**2*I3 + 2*z*(rn3T+n3rT))/rnorm**7
        + 35*z**2*rrT/rnorm**9)
    return dadr_tb + dadr_j2


def j2_acceleration(r_vec, mu):
    """J2 perturbation acceleration [km/s²]."""
    x, y, z = r_vec;  r = np.linalg.norm(r_vec)
    factor = (1.5 * J2 * mu * R_EARTH**2) / r**5
    zx2 = 5 * z**2 / r**2
    return np.array([factor*x*(zx2-1), factor*y*(zx2-1), factor*z*(zx2-3)])