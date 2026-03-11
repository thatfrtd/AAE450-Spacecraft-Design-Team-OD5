import numpy as np
from constants import *

def dynamics_tb_j2(t, y):
    """
    dynamics for two body problem and j2 pertubation. 
    """
    r = y[0:3]
    v = y[3:6]

    rnorm = np.linalg.norm(r)

    # Accelerations
    a_tb = -MU_EARTH * r / rnorm**3
    a_j2 = j2_acceleration(r, MU_EARTH)
    a = a_tb + a_j2

    dydt = np.zeros_like(y)
    dydt[0:3] = v
    dydt[3:6] = a
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
    Finds the jacobian of the two body and J2 at a position of the spacecraft.
    """
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

    G = dadr_tb + dadr_j2
    # System matrix
    A = np.block([
        [O3, I3],
        [G,  O3]
    ])
    return A

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