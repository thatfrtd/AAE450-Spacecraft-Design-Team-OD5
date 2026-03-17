"""
Guidance Law — all quantities in km and km/s.
"""
import numpy as np
from dynamics.attitude import *
import config


def guidance_law(x_est):
    """
    Compute desired translational and attitude states for the chaser.

    Parameters
    ----------
    x_est : (32,) estimated state vector (km, km/s, quaternions, rad/s)

    Returns
    -------
    r_des : (3,)  desired chaser position        [km]
    v_des : (3,)  desired chaser velocity        [km/s]
    q_des : (4,)  desired chaser attitude        [quaternion]
    w_des : (3,)  desired chaser angular velocity [rad/s]
    """
    # Unpack states
    r_T = x_est[0:3]
    v_T = x_est[3:6]
    q_T = x_est[6:10]
    w_T = x_est[10:13]
    r_C = x_est[13:16]
    v_C = x_est[16:19]
    q_C = x_est[19:23]

    # Rotation matrices
    T_T_I = quat2dcm(quat_inv(q_T))   # inertial → target body
    T_D_T = quat2dcm(config.q_D_T)    # target body → docking frame
    T_C_I = quat2dcm(quat_inv(q_C))   # inertial → chaser body

    # Desired attitude: chaser aligns with target
    q_des = quat_inv(q_T)

    # Desired angular velocity: chaser matches target rate expressed in chaser frame
    q_T_C = quat_multiply(q_T, q_C)
    w_des = quat2dcm(q_T_C) @ w_T

    # Desired position [km]:
    #   r_des = r_T + R_TI @ (r_dock_T + R_DT @ r_rel_des_D) - R_CI @ r_attach_C
    # All config distances are in km
    r_des = (r_T
             + T_T_I @ (config.r_dock_T + T_D_T @ config.r_rel_des_D)
             - T_C_I @ config.r_attach_C)

    # Desired velocity [km/s]:
    #   v_des = v_T + R_TI @ (v_rel + w_T × r_dock_T) - R_CI @ (w_C × r_attach_C)
    v_des = (v_T
             + T_T_I @ (T_D_T @ config.v_rel_des_D
                        + np.cross(w_T, config.r_dock_T + T_D_T @ config.r_rel_des_D))
             - T_C_I @ np.cross(config.w_des_C, config.r_attach_C))

    return r_des, v_des, q_des, w_des