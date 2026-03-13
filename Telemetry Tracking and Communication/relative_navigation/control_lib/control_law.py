import numpy as np
from dynamics.attitude import *
import config

def control_law(x_est, w_C, r_des, v_des, q_des, w_des, I_c):
    # Chaser states
    r_C = x_est[13:16]
    v_C = x_est[16:19]
    q_C = x_est[19:23]

    # Gain matricies # TODO: Update at some point
    K_r     = np.eye(3)
    K_v     = np.eye(3)
    K_theta = np.eye(3)
    K_w     = np.eye(3)

    # Translational Control Input
    dv_rqd = K_r @ (r_des - r_C) + K_v @ (v_des - v_C)
    u_cmd = quat2dcm(q_C) @ dv_rqd

    # Torque Control Input
    theta_full = quat_multiply(q_des, quat_inv(q_C))
    theta_des = theta_full[0:3]
    tau_cmd = K_theta @ (theta_des) + K_w @ (w_des - w_C)


    return u_cmd, tau_cmd
