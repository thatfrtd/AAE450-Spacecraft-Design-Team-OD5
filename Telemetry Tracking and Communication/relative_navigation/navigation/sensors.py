import numpy as np
from dynamics.attitude import *
import config

def measure_gyro(w_T, w_C, b_w, ep_w):
    """Gyro measurement with gyro bias. """
    q_w = np.concatenate([ep_w, 1])
    return quat2dcm(q_w) @ ((np.eye(3) + np.diag(config.f_w)) @ w_C + b_w ) 

def measure_star_tracker(q_C):
    """Measurements from star tracker. """
    # TODO: Implement complex model
    return q_C

def measure_optical_camera(r_T, r_C, q_C):
    """Optical camera measuring azimuth and elevation angles."""
    # Find LOS between optical camera and target
    rho_I = r_T - r_C

    T_IC = quat2dcm(q_C)
    rho_C = T_IC @ rho_I

    # Extract azimuth and elevation
    az = np.arctan2(rho_C[1], rho_C[0])
    el = np.arctan2(rho_C[2], np.sqrt(rho_C[0]**2 + rho_C[1]**2))

    # Add measurement noise
    noise = np.random.multivariate_normal(np.zeros(2), config.sigma_alpha)

    return np.array([az, el]) + noise
