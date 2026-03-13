import numpy as np
from utils.constants import *
from dynamics.attitude import *

def target_dynamics(t, y, I, noise_g, noise_omega):
    """
    dynamics for two body problem and j2 pertubation. 
    """
    r = y[0:3]
    v = y[3:6]
    q = y[6:10]
    w = y[10:13]

    rnorm = np.linalg.norm(r)

    # Computing disturbances
    dist_g     = noise_g(t) # np.zeros(3) 
    dist_omega = noise_omega(t) # np.zeros(3) 

    # Accelerations
    a_tb = -MU_EARTH * r / rnorm**3
    a_j2 = j2_acceleration(r, MU_EARTH)
    a = a_tb + a_j2 + dist_g

    # Attitude Dynamics
    q_dot = kde(q, w)
    w_dot = dde(I, w, q, r, dist_omega, np.zeros(3))


    dydt = np.zeros_like(y)
    dydt[0:3] = v
    dydt[3:6] = a
    dydt[6:10] = q_dot
    dydt[10:13] = w_dot
    
    return dydt

def chaser_dynamics(t, y, I, noise_g, noise_omega, u_dv, u_tau):
    """
    dynamics for two body problem and j2 pertubation. 
    """
    r = y[0:3]
    v = y[3:6]
    q = y[6:10]
    w = y[10:13]

    rnorm = np.linalg.norm(r)

    # Computing disturbances
    dist_g     = np.zeros(3) #  noise_g(t)
    dist_omega = np.zeros(3) #  noise_omega(t)

    # Accelerations
    a_tb = -MU_EARTH * r / rnorm**3
    a_j2 = j2_acceleration(r, MU_EARTH)
    a = a_tb + a_j2 + dist_g + u_dv

    # Attitude Dynamics
    q_dot = kde(q, w)
    w_dot = dde(I, w, q, r, dist_omega, u_tau)


    dydt = np.zeros_like(y)
    dydt[0:3] = v
    dydt[3:6] = a
    dydt[6:10] = q_dot
    dydt[10:13] = w_dot
    
    return dydt

def dynamics_truth(t, y, u_applied, tau_applied, I_r, I_c, noise, param_noise):
    """
    Full 6-dof truth dynamics with guidance and control law for chaser rendezvous. 
    """
    dydt = np.zeros_like(y)
    # -- Unpack State ---------------------------------------------------
    # Target State
    r_T = y[0:3]
    v_T = y[3:6]
    q_T = y[6:10]
    w_T = y[10:13]
    # Chaser State
    r_C = y[13:16]
    v_C = y[16:19]
    q_C = y[19:23]
    w_C = y[23:26]
    # Parameter Values
    bias_C = y[26:29] # gyro bias
    ep_S   = y[29:32] # star-camera misalignment
    ep_O   = y[32:35] # optical camera misalignment

    # Unpack Parameters
    sigma_v_T, sigma_omega_T, sigma_v_C, sigma_omega_C = noise
    tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o = param_noise

    # Chaser Equations
    dydt[0:13] = target_dynamics(t, y, I_r, sigma_v_T, sigma_omega_T)

    # Target Equations
    dydt[13:26] = chaser_dynamics(t, y[13:26], I_c, sigma_v_C, sigma_omega_C, u_applied, tau_applied)

    # Parameter Equations
    dydt[26:35] = parameter_dynamics(t, y[26:35], tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o)

    return dydt


def dynamics_predict(t, y, u_applied, tau_applied, I_r, I_c, param_noise):
    """
    Full 6-dof prediction dynamics for rendezvous. 
    """
    dydt = np.zeros_like(y)
    # -- Unpack State ---------------------------------------------------
    # Target State
    r_T = y[0:3]
    v_T = y[3:6]
    q_T = y[6:10]
    w_T = y[10:13]
    # Chaser State
    r_C = y[13:16]
    v_C = y[16:19]
    q_C = y[19:23]
    w_C = y[23:26]

    # Unpack Parameters
    sigma_v_T     = 0
    sigma_omega_T = 0
    sigma_v_C     = 0
    sigma_omega_C = 0
    tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o = param_noise

    # Chaser Equations
    dydt[0:13] = target_dynamics(t, y, I_r, sigma_v_T, sigma_omega_T)

    # Target Equations
    dydt[13:26] = chaser_dynamics(t, y[13:26], I_c, sigma_v_C, sigma_omega_C, u_applied, tau_applied)
    
    # Parameter Equations
    dydt[26:35] = parameter_dynamics(t, y[26:35], tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o)

    return dydt


def parameter_dynamics(t, y, tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o):
    b_w   = y[0:3]
    eps_s = y[3:6]
    eps_o = y[6:9]

    b_w_dot   = -b_w   / tau_b + np.random.normal(0, sigma_b, size=(3,))
    eps_s_dot = -eps_s / tau_s + np.random.normal(0, sigma_s, size=(3,))
    eps_o_dot = -eps_o / tau_o + np.random.normal(0, sigma_o, size=(3,))

    dydt  = np.concatenate([b_w_dot, eps_s_dot, eps_o_dot])
    return dydt


def predict_jacobian(m_minus, w_C, I_t, I_c, q_T, tau_b, tau_s, tau_o):
    """
    Finds the jacobian for 30 parameter state. 
    """
    # Unpack state
    r_T      = m_minus[0:3]
    v_T      = m_minus[3:6]
    theta_T  = m_minus[6:9]
    w_T      = m_minus[9:12]
    r_C      = m_minus[12:15]
    v_C      = m_minus[15:18]
    theta_C  = m_minus[18:21]
    b_w      = m_minus[21:24]
    ep_S     = m_minus[24:27]
    ep_O     = m_minus[27:30]

    # DCM
    T_I_T = quat2dcm(q_T)
    T_T_I = quat2dcm(quat_inv(q_T))

    # Helpers
    r_norm  = np.linalg.norm(r_T)
    r_hat_B = T_I_T @ r_T                        # r_T expressed in target body frame
    skew_rB = skew(r_hat_B)
    I_inv   = np.linalg.inv(I_t)

    # Target State Jacobian
    F_T1 = np.eye(3)
    F_T2 = g_t_Jacobian(r_T)
    F_T3 = -skew(w_T)
    F_T4 = np.eye(3)
    term1 = (3 * MU_EARTH / r_norm**5) * I_inv @ (
        skew_rB @ I_t @ T_I_T - skew(I_t @ r_hat_B) @ T_I_T
    )
    term2 = (-15 * MU_EARTH / r_norm**7) * I_inv @ (
        skew_rB @ I_t @ r_hat_B[:, None] @ r_T[None, :] @ T_I_T.T
    )   # outer product gives (3×3)
    F_T5 = term1 + term2
    F_T6 = (3 * MU_EARTH / r_norm**5) * I_inv @ (
        skew_rB @ I_t @ skew_rB - skew(I_t @ r_hat_B) @ skew_rB
    )

    F_T7 = np.linalg.solve(I_t, skew(I_t @ w_T) - skew(w_T) @ I_t)
    F_TT = np.block([
        [np.zeros((3,3)), F_T1,             np.zeros((3,3)), np.zeros((3,3))],  # δṙ_T
        [F_T2,            np.zeros((3,3)),  np.zeros((3,3)), np.zeros((3,3))],  # δv̇_T
        [np.zeros((3,3)), np.zeros((3,3)),  F_T3,            F_T4           ],  # δθ̇_T
        [F_T5,            np.zeros((3,3)),  F_T6,            F_T7           ],  # δω̇_T
    ])

    # Chaser State Jacobian
    F_C1 = np.eye(3)
    F_C2 = g_t_Jacobian(r_C)
    F_C3 = -skew((w_C - b_w))
    F_CC = np.block([
        [np.zeros((3,3)), F_C1,            np.zeros((3,3))],
        [F_C2           , np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), F_C3],
    ])
    # Relationship between chaser and parameter Jacobian
    F_CP1 = -np.eye(3)
    F_CP = np.block([
        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), np.zeros((3,3))],
        [F_CP1          , np.zeros((3,3)), np.zeros((3,3))],
    ])
    # Parameter Jacobian
    F_P1 = (-1/tau_b) * np.eye(3)
    F_P2 = (-1/tau_s) * np.eye(3)
    F_P3 = (-1/tau_o) * np.eye(3)
    F_PP = np.block([
        [F_P1, np.zeros((3,3)), np.zeros((3,3))],
        [np.zeros((3,3)), F_P2, np.zeros((3,3))],
        [np.zeros((3,3)), np.zeros((3,3)), F_P3],
    ])

    # Pack Jacobian
    F = np.zeros((30, 30)) # Initialize Jacobian
    # F_TT (12×12) — top left
    F[0:12, 0:12] = F_TT

    # F_CC (9×9) — middle
    F[12:21, 12:21] = F_CC

    # F_CP (9×9) — middle right
    F[12:21, 21:30] = F_CP

    # F_PP (9×9) — bottom right
    F[21:30, 21:30] = F_PP

    return F

def g_t_Jacobian(r):

    rnorm = np.linalg.norm(r)
    z = r[2]

    I3 = np.eye(3)
    O3 = np.zeros((3, 3))
    rrT = np.outer(r, r)

    # Two-body Jacobian
    dadr_tb = MU_EARTH * (3 * rrT / rnorm**5 - I3 / rnorm**3)

    # J2 Jacobian
    n3 = np.array([0.0, 0.0, 1.0])
    n3n3T = np.outer(n3, n3)
    rn3T = np.outer(r, n3)
    n3rT = np.outer(n3, r)

    dadr_j2 = (
        -3 * 1.0826e-3 * MU_EARTH * R_EARTH**2 / 2
    ) * (
        I3 / rnorm**5
        + 2 * n3n3T / rnorm**5
        - 5 * (rrT + z**2 * I3 + 2*z*(rn3T + n3rT)) / rnorm**7
        + 35 * z**2 * rrT / rnorm**9
    )

    F = dadr_tb + dadr_j2

    return F

def j2_acceleration(r_vec, mu):
    """
    Compute the acceleration due to the J2 zonal harmonic.

    Parameters
    ----------
    r_vec : array_like, shape (3,)
        Position vector [x, y, z] in km
    mu : float
        Gravitational parameter [km^3/s^2]

    Returns
    -------
    a_J2 : ndarray, shape (3,)
        Acceleration due to J2 in inertial frame [km/s^2]
    """
    x, y, z = r_vec
    Re = R_EARTH
    r = np.linalg.norm(r_vec)
    
    factor = (1.5 * J2 * mu * Re**2) / r**5
    zx2 = 5 * z**2 / r**2
    
    ax = factor * x * (zx2 - 1)
    ay = factor * y * (zx2 - 1)
    az = factor * z * (zx2 - 3)
    
    return np.array([ax, ay, az])

