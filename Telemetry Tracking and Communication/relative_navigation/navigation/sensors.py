"""
Sensor Measurement Models
=========================
All positions in km, velocities in km/s, angles in radians.

Key design: LIDAR and optical camera are combined into a single
relative position measurement.  The camera provides az/el angles
and the LIDAR provides range.  Together they give a 3D relative
position vector r^T - r^C in the inertial frame.

  measure_gyro              → (3,)  angular velocity [rad/s]
  measure_star_tracker      → (4,)  quaternion q^{I→C}
  measure_relative_position → (3,)  r^T - r^C [km, inertial]
  measure_gps               → (3,)  r^C [km, inertial]
  measure_target_attitude   → (4,)  quaternion q^{I→T}
"""
import numpy as np
from dynamics.attitude import *
import config


def measure_gyro(w_C, b_w, S_s, ep_w):
    """
    Gyro measurement model.
    w_C  : (3,) true chaser angular velocity [rad/s]
    b_w  : (3,) estimated gyro bias [rad/s]
    S_s  : (3,) estimated scale factors
    ep_w : (3,) estimated misalignment (small rotation vector)
    Returns corrupted angular velocity (3,).
    """
    q_w = np.concatenate([ep_w, [1.0]])
    q_w /= np.linalg.norm(q_w)
    C_misalign = quat2dcm(q_w)
    scale_mat  = np.eye(3) + np.diag(S_s)
    return C_misalign @ (scale_mat @ w_C + b_w)


def measure_star_tracker(q_C):
    """
    Star tracker — noisy quaternion measurement of chaser attitude.
    q_C : (4,) true q^{I→C}
    Returns (4,) measured quaternion.
    """
    sigma   = config.sigma_st
    d_theta = np.random.normal(0, sigma, 3)
    angle   = np.linalg.norm(d_theta)
    if angle < 1e-10:
        dq = np.array([0., 0., 0., 1.])
    else:
        axis = d_theta / angle
        dq   = np.concatenate([np.sin(angle/2)*axis, [np.cos(angle/2)]])
    q_meas = quat_multiply(q_C, dq)
    return q_meas / np.linalg.norm(q_meas)


def measure_relative_position(r_T, r_C, q_C):
    """
    Combined LIDAR + optical camera measurement.

    The camera provides bearing angles (az, el) and the LIDAR provides
    range.  Together they reconstruct the 3D relative position vector
    r^T - r^C expressed in the INERTIAL frame.

    Noise model
    -----------
    At range ρ with angle noise σ_θ and range noise σ_r:
      - Along-LOS error:    σ_r  [km]
      - Cross-LOS error:    ρ * σ_θ  [km]

    We model this as isotropic with σ_rel = sqrt(σ_r² + (ρ*σ_θ)²)
    for simplicity, which is conservative (slightly overestimates
    cross-LOS error at short range, underestimates at long range).

    Parameters
    ----------
    r_T : (3,) true target position [km, inertial]
    r_C : (3,) true chaser position [km, inertial]
    q_C : (4,) true chaser attitude q^{I→C}

    Returns
    -------
    z : (3,) noisy measurement of r^T - r^C [km, inertial]
    """
    rho_I_true = r_T - r_C
    rho        = np.linalg.norm(rho_I_true)

    # Combined noise: range noise + angle-induced cross-LOS noise
    sigma_along = config.sigma_lidar                       # km
    sigma_cross = rho * config.sigma_opt_angle             # km = range * angle_noise
    sigma_rel   = np.sqrt(sigma_along**2 + sigma_cross**2) # isotropic approximation

    noise = np.random.multivariate_normal(np.zeros(3),
                                          sigma_rel**2 * np.eye(3))
    return rho_I_true + noise                              # (3,) inertial


def measure_gps(r_C):
    """
    GPS — chaser absolute inertial position.
    r_C : (3,) true chaser position [km]
    Returns (3,) noisy position [km].
    """
    noise = np.random.multivariate_normal(np.zeros(3), config.R_GPS)
    return r_C + noise


def measure_target_attitude(q_T):
    """
    Target attitude proxy — simulates feature-tracking pose estimation.
    q_T : (4,) true target attitude q^{I→T}
    Returns (4,) noisy quaternion measurement.
    """
    sigma   = config.sigma_target_att
    d_theta = np.random.normal(0, sigma, 3)
    angle   = np.linalg.norm(d_theta)
    if angle < 1e-10:
        dq = np.array([0., 0., 0., 1.])
    else:
        axis = d_theta / angle
        dq   = np.concatenate([np.sin(angle/2)*axis, [np.cos(angle/2)]])
    q_meas = quat_multiply(q_T, dq)
    return q_meas / np.linalg.norm(q_meas)