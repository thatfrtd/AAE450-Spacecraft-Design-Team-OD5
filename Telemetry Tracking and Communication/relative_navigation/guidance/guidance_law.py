import numpy as np
from dynamics.attitude import *
import config


def guidance_law(x_est):
    # Unpack Target State
    r_T = x_est[0:3]
    v_T = x_est[3:6]
    q_T = x_est[6:10]
    w_T = x_est[10:13]
    # Unpack Chaser State
    r_C = x_est[13:16]
    v_C = x_est[16:19]
    q_C = x_est[19:23]

    # Desired quaternion
    q_des = quat_inv(q_T)

    # Desired angular velocity
    q_T_C = quat_multiply(q_T, q_C) # TODO: Come back and check the order! 
    w_des = quat2dcm(q_T_C) @ w_T

    # Desired positon
    T_T_I = quat2dcm(quat_inv(q_T))
    T_D_T = quat2dcm(config.q_D_T)
    T_C_I = quat2dcm(quat_inv(q_C))
    r_des = r_T + T_T_I @ (config.r_dock_T + T_D_T @ config.r_rel_des_D) - T_C_I @ config.r_attach_C

    # Desired Velocity
    v_des = v_T + T_T_I @ (T_D_T @ config.v_rel_des_D + np.cross(w_T, (config.r_dock_T + T_D_T @ config.r_rel_des_D))) - T_C_I @ (np.cross(config.w_des_C, config.r_attach_C))

    return r_des, v_des, q_des, w_des