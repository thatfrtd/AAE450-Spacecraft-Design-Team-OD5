import numpy as np
from dynamics.attitude import *
import config


def measure_gyro(w_C, b_w, S_s, ep_w):
    """
    Gyro measurement model.

    Parameters
    ----------
    w_C  : (3,) true angular velocity of chaser (from truth state)
    b_w  : (3,) estimated gyro bias
    S_s  : (3,) estimated scale factor (diagonal elements)
    ep_w : (3,) estimated misalignment (small-angle rotation vector)

    Returns
    -------
    omega_gyro : (3,) corrupted angular velocity measurement
    """
    # Bug fix: np.concatenate needs a list of arrays, scalar 1 → [1.0]
    q_w = np.concatenate([ep_w, [1.0]])
    q_w = q_w / np.linalg.norm(q_w)           # normalise just in case

    C_misalign = quat2dcm(q_w)                 # small misalignment DCM
    scale_mat  = np.eye(3) + np.diag(S_s)      # scale factor matrix

    omega_gyro = C_misalign @ (scale_mat @ w_C + b_w)

    return omega_gyro


def measure_star_tracker(q_C):
    """
    Star tracker measurement — returns noisy quaternion.

    Parameters
    ----------
    q_C : (4,) true chaser attitude quaternion [q1,q2,q3,q4]

    Returns
    -------
    z_st : (4,) measured quaternion with additive attitude noise
    """
    # Add small attitude noise (rotation-vector perturbation)
    sigma = config.sigma_st if hasattr(config, 'sigma_st') else np.deg2rad(0.01)
    d_theta = np.random.normal(0, sigma, 3)

    angle = np.linalg.norm(d_theta)
    if angle < 1e-10:
        dq = np.array([0.0, 0.0, 0.0, 1.0])
    else:
        axis = d_theta / angle
        dq   = np.concatenate([np.sin(angle/2) * axis, [np.cos(angle/2)]])

    # q_meas = q_true ⊗ δq  (right-multiply noise)
    q_meas = quat_multiply(q_C, dq)
    return q_meas / np.linalg.norm(q_meas)


def measure_optical_camera(r_T, r_C, q_C):
    """
    Optical camera measuring azimuth and elevation angles to the target.

    Parameters
    ----------
    r_T : (3,) target position in inertial frame
    r_C : (3,) chaser position in inertial frame
    q_C : (4,) chaser attitude quaternion [q1,q2,q3,q4]

    Returns
    -------
    z : (2,) [azimuth, elevation] in radians
    """
    rho_I = r_T - r_C

    T_IC  = quat2dcm(q_C)
    rho_C = T_IC @ rho_I

    az = np.arctan2(rho_C[1], rho_C[0])
    el = np.arctan2(rho_C[2], np.sqrt(rho_C[0]**2 + rho_C[1]**2))

    # Bug fix: R_OPT_CAM is the (2×2) covariance; use it directly
    noise = np.random.multivariate_normal(np.zeros(2), config.R_OPT_CAM)

    return np.array([az, el]) + noise

def measure_gps(r_C_true):
    noise = np.random.multivariate_normal(np.zeros(3), config.R_GPS)
    return r_C_true + noise

def measure_lidar(r_T, r_C):
    """LIDAR range measurement [km]."""
    true_range = np.linalg.norm(r_T - r_C)
    noise = np.random.normal(0, config.sigma_lidar)
    return true_range + noise