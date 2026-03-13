import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import config
from scipy.integrate import solve_ivp
from utils.constants import *
from dynamics.dynamics import *
from analysis.plotting import *
import matplotlib.pyplot as plt
from utils.helper import *
from control_lib.actuator_model import *
from guidance.guidance_law import *
from control_lib.control_law import *
from navigation.sensors import *
from navigation.predict import *


np.random.seed(0) # Set random seed

def main():
    # Run monte carlo using gnc_sim call. 
    gnc_sim()

def gnc_sim():
    # -- Initialization -----------------------
    # Target Init
    target_translational_state = np.concatenate([config.tar_pos, config.tar_vel]) # Translational State
    target_attitude = np.concatenate([config.TARGET_EP, config.TARGET_ANG_VEL]) # Attitude State
    target_state = np.concatenate([target_translational_state, target_attitude]) # Full State
    I_r = config.TARGET_I # MOI Tensor

    # Chaser Init
    chaser_translational_state = np.concatenate([config.chaser_pos, config.chaser_vel]) # Translational State
    chaser_attitude = np.concatenate([config.CHASER_EP, config.CHASER_ANG_VEL]) # Attitude State
    chaser_state = np.concatenate([chaser_translational_state, chaser_attitude]) # Full State
    I_c = config.CHASER_I # MOI Tensor

    # Bias Parameters
    param_init = np.concatenate([config.b_w_c, config.ep_S_S, config.ep_O_O])

    # Model parameters


    # -- Trajectory History Storage ----------- 
    last_opt_cam_meas = 0
    last_star_track_meas = 0
    last_gyro_meas = 0

    hist_target_truth   = np.zeros((config.N, 13))
    hist_chaser_truth   = np.zeros((config.N, 10))
    hist_u_cmd          = np.zeros((config.N, 3))
    hist_tau_cmd        = np.zeros((config.N, 3))
    hist_u_input        = np.zeros((config.N, 3))
    hist_tau_input      = np.zeros((config.N, 3))
    hist_x_est          = np.zeros((config.N, 32))
    hist_P_diag         = np.zeros((config.N, 30))
    hist_z_opt          = np.zeros((config.N, 2))
    hist_z_st           = np.zeros((config.N, 4))


    # -- 6-dof sim -----------------------------

    # Initial True State Vector
    x_true = np.concatenate([target_state, chaser_state])
    
    # Initial Estimation State Vector
    x_est = np.concatenate([
        target_state  + np.concatenate([
            config.sigma_r_T  * np.random.randn(3),
            config.sigma_v_T  * np.random.randn(3),
            delta_theta_to_quat(config.sigma_th_T * np.random.randn(3)),
            config.sigma_w_T  * np.random.randn(3)
        ]),
        chaser_state + np.concatenate([
            config.sigma_r_C  * np.random.randn(3),
            config.sigma_v_C  * np.random.randn(3),
            delta_theta_to_quat(config.sigma_th_C * np.random.randn(3)),
        ]),
        param_init + np.concatenate([
            config.sigma_bw   * np.random.randn(3),
            config.sigma_Ss   * np.random.randn(3),
            config.sigma_Oo   * np.random.randn(3),
        ])
    ])

    # Initial Estimation Covariance
    P0 = np.diag(np.concatenate([
        (config.sigma_r_T  * np.ones(3))**2,
        (config.sigma_v_T  * np.ones(3))**2,
        (config.sigma_th_T * np.ones(3))**2,
        (config.sigma_w_T  * np.ones(3))**2,
        (config.sigma_r_C  * np.ones(3))**2,
        (config.sigma_v_C  * np.ones(3))**2,
        (config.sigma_th_C * np.ones(3))**2,
        (config.sigma_bw   * np.ones(3))**2,
        (config.sigma_Ss   * np.ones(3))**2,
        (config.sigma_Oo   * np.ones(3))**2,
    ]))
    P = P0.copy()

    # Noise initialization
    noise  = np.concatenate([config.sigma_v_T, config.sigma_omega_T, config.sigma_v_C, config.sigma_omega_C])
    param_noise = np.concatenate([config.tau_b, config.tau_s, config.tau_o, config.sigma_b, config.sigma_s, config.sigma_o])

    for k in range(config.N):
        # Sim step time
        t_start, t_end = k * config.SIM_DT, (k+1) * config.SIM_DT
        tspan = [t_start, t_end]

        # -- Guidance --------------------------------------------
        # Compute Desired States using Guidance Law
        r_des, v_des, q_des, w_des = guidance_law(x_est)

        # -- Control ----------------------------------------------
        # Compute Control Input Command 
        w_C = x_true[23:26]
        u_cmd, tau_cmd = control_law(x_est, w_C, r_des, v_des, q_des, w_des, I_c)
        
        # -- Actuator models --------------------------------------
        u_applied = translational_control(u_cmd)
        tau_applied = momentum_wheel_model(tau_cmd)

        # -- Simulate Truth Chaser and Target State ---------------
        sol = solve_ivp(
            dynamics_truth, 
            tspan, x_true, method='RK45', 
            args=(u_applied, tau_applied, I_r, I_c, noise, param_noise, ), 
            rtol=config.tol, atol=config.tol,
        )
        x_true = sol.y[:, -1]
        x_true[6:10]  = normalize_quat(x_true[6:10]) # Seems like the wrong thing to do 
        x_true[19:23] = normalize_quat(x_true[19:23])

        # -- Sensor Measurements (z) ------------------------------
        t_now = t_end
        omega_gyro = None
        z_st       = None
        z_opt      = None
        # Chaser Gyro
        # Feeds into predict step as a input param not part of estimated state
        if (t_now - last_gyro_meas) >= config.DT_GYRO:
            omega_gyro = measure_gyro(x_est[23:26], x_est[26:29], x_est[29:32])
            last_gyro_meas = t_now

        # Chaser Star Tracker
        # if 10 seconds have passed since last measurement get measurement
        # two seperate updates one for star tracker one for opt cam
        if (t_now - last_star_track_meas) >= config.DT_STAR_TRACKER:
            z_st = measure_star_tracker(x_true[19:23])
            last_star_track_meas = t_now


        # Chaser Optical Camera
        # if 60 seconds have passed since last measurement get measurement
        if (t_now - last_opt_cam_meas) >= config.DT_OPT_CAM:
            z_opt = measure_optical_camera(x_true[0:3], x_true[13:16],
                                           x_true[19:23])
            last_opt_cam_meas = t_now
        
        # -- Navigation -------------------------------------------
        # Predict Step
        x_est, P = ekf_predict(x_est, P, omega_gyro, u_applied, tau_applied, I_r, I_c, noise, param_noise, config.SIM_DT)

        # Update Step
        if z_st is not None:
            x_est, P, _ = ekf_update_star_tracker(x_est, P, z_st, config.R_STAR_TRACKER)
 
        if z_opt is not None:
            x_est, P, _ = ekf_measurement_update(x_est, P, z_opt, config.R_OPT_CAM)
 
        x_est[6:10]  = normalize_quat(x_est[6:10])
        x_est[19:23] = normalize_quat(x_est[19:23])
        
        # -- Compute Delta e --------------------------------------
        delta_r = x_true[13:16] - x_true[0:3] - (x_est[13:16] - x_est[0:3])
        delta_q = quat_multiply(x_true[19:23], quat_inv(x_est[19:23]))
 

        # -- Simulation Variable Storage --------------------------
        # Save target Truth trajectory history
        hist_target_truth[k, :]  = x_true[0:13]

        # Save chaser Truth trajectory history
        hist_chaser_truth[k, :]  = x_true[13:23]

        # Save Control Command
        hist_u_cmd[k, :]         = u_cmd
        hist_tau_cmd[k, :]       = tau_cmd


        # Save Control input (delayed due to actuator model)
        hist_u_input[k, :]       = u_applied
        hist_tau_input[k, :]     = tau_applied

        # Save Estimated State (m)
        hist_x_est[k, :]         = x_est

        # Save Covariance (P)
        hist_P_diag[k, :]        = np.diag(P)

        # Save measurements
        if z_opt is not None:
            hist_z_opt[k, :] = z_opt
        if z_st is not None:
            hist_z_st[k, :] = z_st

        # -- Check for Rendezvous Sucess ---------------------------
        pos_err = np.linalg.norm(x_true[13:16] - x_true[0:3])
        att_err = 2 * np.arccos(np.clip(abs(quat_multiply(
                    x_true[19:23], quat_inv(x_true[6:10]))[3]), 0, 1))
 
        if pos_err < config.POS_TOL and att_err < config.ATT_TOL:
            print(f"Rendezvous achieved at t = {t_end:.1f} s (k={k})")
            break

    print("Simulation Complete after {k} time steps. ")
    return True
    """
    return {
        'target':    hist_target_truth[:k+1],
        'chaser':    hist_chaser_truth[:k+1],
        'u_cmd':     hist_control_cmd[:k+1],
        'u_applied': hist_control_input[:k+1],
        'x_est':     hist_x_est[:k+1],
        'P_diag':    hist_P_diag[:k+1],
        'z_opt':     hist_z_opt[:k+1],
        'z_st':      hist_z_st[:k+1],
    }
    """



if __name__ == "__main__":
    main()