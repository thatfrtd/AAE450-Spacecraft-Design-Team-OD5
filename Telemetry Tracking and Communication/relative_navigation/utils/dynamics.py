import numpy as np
from utils.constants import *
from utils.attitude import *

def dynamics(t, y, I, noise_g, noise_alpha, ):
    """
    dynamics for two body problem and j2 pertubation. 
    """
    r = y[0:3]
    v = y[3:6]
    q = y[6:10]
    w = y[10:13]

    rnorm = np.linalg.norm(r)

    # Accelerations
    a_tb = -MU_EARTH * r / rnorm**3
    a_j2 = j2_acceleration(r, MU_EARTH)
    a = a_tb + a_j2

    # Attitude Dynamics
    q_dot = kde(q, w)
    w_dot = dde(I, w, q, r, noise_alpha)


    dydt = np.zeros_like(y)
    dydt[0:3] = v
    dydt[3:6] = a
    dydt[6:10] = q_dot
    dydt[10:13] = w_dot
    
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

def grav_gradient(I, w, q, r):
    """
    Gravity gradient torque in the body frame.

    Args:
        I : (3,3) inertia tensor, body frame
        w : (3,)  angular velocity, body frame  (unused but kept for signature)
        q : (4,)  quaternion [qx, qy, qz, qw]
        r : (3,)  ECI position, km
    """
    # Rotate ECI position into body frame using quaternion
    r_body = eci_to_body(q, r)

    r_norm = np.linalg.norm(r_body)
    r_hat  = r_body / r_norm
    L_g    = (3.0 * MU_EARTH / r_norm**3) * np.cross(r_hat, I @ r_hat)
    return L_g


def eci_to_body(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """
    Rotate a vector from ECI to body frame using quaternion [qx, qy, qz, qw].
    """
    qx, qy, qz, qw = q
    # Rotation matrix from ECI to body  (transpose of body-to-ECI)
    R = np.array([
        [1 - 2*(qy**2 + qz**2),     2*(qx*qy + qw*qz),     2*(qx*qz - qw*qy)],
        [    2*(qx*qy - qw*qz), 1 - 2*(qx**2 + qz**2),     2*(qy*qz + qw*qx)],
        [    2*(qx*qz + qw*qy),     2*(qy*qz - qw*qx), 1 - 2*(qx**2 + qy**2)],
    ])
    return R @ v

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