import numpy as np
import math

def solve_kepler(M, e, tol=1e-10, max_iter=50):
    """
    Solve Kepler's equation M = E - e sin(E) for eccentric anomaly E.

    Parameters
    ----------
    M : float
        Mean anomaly [rad]
    e : float
        Eccentricity (0 <= e < 1)
    tol : float
        Convergence tolerance
    max_iter : int
        Maximum Newton-Raphson iterations

    Returns
    -------
    E : float
        Eccentric anomaly [rad]
    """
    # Norm to [0, 2pi]
    M = np.mod(M, 2*np.pi)

    # Initial Guess
    E = M
    
    # Iterate
    for _ in range(max_iter):
        f = E - e*np.sin(E) - M
        f_prime = 1 - e*np.cos(E)

        dE = -f / f_prime
        E += dE
        if abs(dE) < tol:
            return E
        
    raise RuntimeError("Kepler solver did not converge. ")

def kepler_eqn(E, e):
    """
    Classic Kepler Equation
    
    Parameters
    ----------
    E : float
        Eccentric Anomaly [rad]
    e : float
        Eccentrity

    Returns
    -------
    M : float
        Mean Anomaly [rad]
    """
    M = E - e * np.sin(E)
    M = np.mod(M, 2 * np.pi)
    return M

def eccentric_to_true(E, e):
    """
    Convert Eccentric Anomaly to True Anomaly (nu). 

    Parameters
    ----------
    E : float
        Eccentric anomaly [rad]
    e : float
        Eccentricity 0 < e < 1. 

    Returns
    -------
    nu : float
        True anomaly [rad]
    """
    nu = 2 * math.atan2(math.sqrt(1+e)*math.sin(E/2), math.sqrt(1-e) * math.cos(E/2))
    nu = np.mod(nu, 2*np.pi)
    return nu

def true_to_eccentric(nu, e):
    """
    Convert True Anomaly (nu) to Eccentric Anomaly (E). 

    Parameters
    ----------
    nu : float
        True anomaly [rad]
    e : float
        Eccentricity 0 < e < 1. 

    Returns
    -------
    E : float
        Eccentric anomaly [rad]
    """
    E = 2 * np.arctan2(np.sqrt(1-e) * np.sin(nu/2), np.sqrt(1+e) * np.cos(nu/2))
    E = np.mod(E, 2*np.pi)
    return E