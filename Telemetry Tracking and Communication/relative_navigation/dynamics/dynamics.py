import numpy as np
from utils.constants import *
from dynamics.attitude import *


def target_dynamics(t, y, I, u_tau=None):
    """
    Deterministic two-body + J2 dynamics for the target.
    No stochastic noise — noise enters only through the EKF Q matrix.
    """
    r = y[0:3];  v = y[3:6];  q = y[6:10];  w = y[10:13]
    rnorm = np.linalg.norm(r)

    a = -MU_EARTH * r / rnorm**3 + j2_acceleration(r, MU_EARTH)

    q_dot = kde(q, w)
    w_dot = dde(I, w, q, r, np.zeros(3), np.zeros(3))

    dydt = np.zeros_like(y)
    dydt[0:3]   = v
    dydt[3:6]   = a
    dydt[6:10]  = q_dot
    dydt[10:13] = w_dot
    return dydt


def chaser_dynamics(t, y, I, u_dv, u_tau):
    """
    Deterministic two-body + J2 dynamics for the chaser.
    u_dv : (3,) thrust in chaser BODY frame [km/s²] — rotated to inertial here.
    u_tau: (3,) torque [kg·km²/s²]
    """
    r = y[0:3];  v = y[3:6];  q = y[6:10];  w = y[10:13]
    rnorm = np.linalg.norm(r)

    T_CI = quat2dcm(quat_inv(q))          # body → inertial
    u_dv_inertial = T_CI @ u_dv

    a = -MU_EARTH * r / rnorm**3 + j2_acceleration(r, MU_EARTH) + u_dv_inertial

    q_dot = kde(q, w)
    w_dot = dde(I, w, q, r, np.zeros(3), u_tau)

    dydt = np.zeros_like(y)
    dydt[0:3]   = v
    dydt[3:6]   = a
    dydt[6:10]  = q_dot
    dydt[10:13] = w_dot
    return dydt


def parameter_dynamics_det(t, y, tau_b, tau_s, tau_o):
    """
    Deterministic parameter dynamics — only the mean-reversion term.
    Stochastic driving noise is accounted for in the Q matrix, NOT here.
    This makes the integrator deterministic and RK45-safe.
    """
    b_w   = y[0:3];  eps_s = y[3:6];  eps_o = y[6:9]
    return np.concatenate([
        -b_w   / tau_b,
        -eps_s / tau_s,
        -eps_o / tau_o,
    ])


def dynamics_truth(t, y, u_applied, tau_applied, I_r, I_c, param_noise):
    """
    Full 35-element DETERMINISTIC truth dynamics.
    target(13) + chaser(13) + params(9)
    Noise enters ONLY through sensor measurement functions (sensors.py).
    """
    dydt = np.zeros_like(y)

    tau_b = param_noise[0:3]
    tau_s = param_noise[3:6]
    tau_o = param_noise[6:9]

    dydt[0:13]  = target_dynamics(t, y[0:13],  I_r)
    dydt[13:26] = chaser_dynamics(t, y[13:26], I_c, u_applied, tau_applied)
    dydt[26:35] = parameter_dynamics_det(t, y[26:35], tau_b, tau_s, tau_o)
    return dydt


def dynamics_predict(t, y, u_applied, tau_applied, I_r, I_c, param_noise):
    """
    Full 35-element DETERMINISTIC predict dynamics for EKF propagation.
    Identical to dynamics_truth — both are now deterministic.
    """
    dydt = np.zeros_like(y)

    tau_b = param_noise[0:3]
    tau_s = param_noise[3:6]
    tau_o = param_noise[6:9]

    dydt[0:13]  = target_dynamics(t, y[0:13],  I_r)
    dydt[13:26] = chaser_dynamics(t, y[13:26], I_c, u_applied, tau_applied)
    dydt[26:35] = parameter_dynamics_det(t, y[26:35], tau_b, tau_s, tau_o)
    return dydt


# ── Keep old parameter_dynamics for backwards compat (not used in integrators)
def parameter_dynamics(t, y, tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o):
    """Legacy — do not call from inside solve_ivp."""
    return parameter_dynamics_det(t, y, tau_b, tau_s, tau_o)


def predict_jacobian(m_minus, w_C, I_t, I_c, q_T, tau_b, tau_s, tau_o):
    """Linearised dynamics Jacobian F = ∂f/∂δx  (30×30 error-state)."""
    r_T = m_minus[0:3];  w_T = m_minus[9:12]
    r_C = m_minus[12:15];  b_w = m_minus[21:24]

    T_I_T  = quat2dcm(q_T)
    I_inv  = np.linalg.inv(I_t)
    r_norm = np.linalg.norm(r_T)
    r_hat_B = T_I_T @ r_T
    skew_rB = skew(r_hat_B)

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
        [np.zeros((3,3)), np.eye(3),       np.zeros((3,3)), np.zeros((3,3))],
        [F_T2,            np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), F_T3,            np.eye(3)      ],
        [F_T5,            np.zeros((3,3)), F_T6,            F_T7           ],
    ])

    F_CC = np.block([
        [np.zeros((3,3)), np.eye(3),          np.zeros((3,3))],
        [g_t_Jacobian(r_C), np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)),     -skew(w_C - b_w)],
    ])

    F_CP = np.block([
        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [-np.eye(3),      np.zeros((3,3)), np.zeros((3,3))],
    ])

    F_PP = np.block([
        [(-1/tau_b[0])*np.eye(3), np.zeros((3,3)),         np.zeros((3,3))        ],
        [np.zeros((3,3)),          (-1/tau_s[0])*np.eye(3), np.zeros((3,3))        ],
        [np.zeros((3,3)),          np.zeros((3,3)),          (-1/tau_o[0])*np.eye(3)],
    ])

    F = np.zeros((30, 30))
    F[0:12,  0:12]  = F_TT
    F[12:21, 12:21] = F_CC
    F[12:21, 21:30] = F_CP
    F[21:30, 21:30] = F_PP
    return F


def g_t_Jacobian(r):
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
    x, y, z = r_vec;  r = np.linalg.norm(r_vec)
    factor = (1.5 * J2 * mu * R_EARTH**2) / r**5
    zx2 = 5 * z**2 / r**2
    return np.array([factor*x*(zx2-1), factor*y*(zx2-1), factor*z*(zx2-3)])