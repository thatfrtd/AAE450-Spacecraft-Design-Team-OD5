import numpy as np
from utils.constants import *
from utils.attitude import *

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

def truth_dynamics(t, y, I_t, noise_g, noise_omega, I_c, noise_g, noise_omega, u_dv, u_tau):

    return True

def parameter_dynamics(t, y, tau_b, tau_s, tau_o, sigma_b, sigma_s, sigma_o):
    b_w   = y[0:3]
    eps_s = y[3:6]
    eps_o = y[6:9]

    b_w_dot   = -b_w   / tau_b + np.random.normal(0, sigma_b, size=(3,))
    eps_s_dot = -eps_s / tau_s + np.random.normal(0, sigma_s, size=(3,))
    eps_o_dot = -eps_o / tau_o + np.random.normal(0, sigma_o, size=(3,))

    dydt  = np.concatenate([b_w_dot, eps_s_dot, eps_o_dot])
    return dydt

def dynamics_with_stm(t, y):
    """
    two body and J2 dynamics to get STM. 
    """
    # TODO: Clean up code and document
    r = y[0:3]
    v = y[3:6]
    Phi = y[6:].reshape((6, 6))

    rnorm = np.linalg.norm(r)

    # Accelerations
    a_tb = -MU_EARTH * r / rnorm**3
    a_j2 = j2_acceleration(r, MU_EARTH)
    a = a_tb + a_j2

    # Jacobian of two body and J2
    A = dynamics_Jacobian(r)

    Phi_dot = A @ Phi

    dydt = np.zeros_like(y)
    dydt[0:3] = v
    dydt[3:6] = a
    dydt[6:] = Phi_dot.flatten()

    return dydt

def dynamics_Jacobian(r):
    """
    Finds the jacobian for 32 parameter state
    """
    
    return True

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