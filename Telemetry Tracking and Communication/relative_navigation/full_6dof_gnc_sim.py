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
from navigation.update import *
from tqdm import tqdm
import time


np.random.seed(0) # Set random seed

def main():
    # Run monte carlo using gnc_sim call. 
    results = gnc_sim()

    # Extract histories
    target_hist  = results['target']   # (N, 13)
    chaser_hist  = results['chaser']   # (N, 10) ← NOTE: only 10 states stored
    N_sim        = target_hist.shape[0]

    # Extract estimated states
    x_est_hist = results['x_est']   # (N, 32)

    # Split estimated target & chaser
    target_est  = x_est_hist[:, 0:13]
    chaser_est  = x_est_hist[:, 13:23]

    # Transpose for plotting
    target_est_states = target_est.T
    chaser_est_states = chaser_est.T

    # Build time vector
    t = np.arange(N_sim) * config.SIM_DT

    # ---- Fix shape for plotting ----
    target_states  = target_hist.T  # (13, N)
    
    # Your chaser only stores 10 states → pad to 13 for plotting compatibility
    chaser_states = np.zeros((13, N_sim))
    chaser_states[0:10, :] = chaser_hist.T

    # ---- Plot ----
    plot_relative_state(t, target_states, chaser_states)
    plot_relative_state_overlay(
        t,
        target_states,
        chaser_states,
        target_est_states,
        chaser_est_states
    )

    plot_orbits_3d(t, target_states, chaser_states)

    plt.show()

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
    x_true = np.concatenate([target_state, chaser_state, param_init])  # 13 + 13 + 9 = 35
    
    # Initial Estimation State Vector
    x_est = np.concatenate([
        target_state  + np.concatenate([
            config.sigma_r_T  * np.random.randn(3),
            config.sigma_v_T  * np.random.randn(3),
            delta_theta_to_quat(config.sigma_th_T * np.random.randn(3)),
            config.sigma_w_T  * np.random.randn(3)
        ]),
        chaser_state[:10] + np.concatenate([
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

    # Noise initialization  (use truth noise values sigma_vel_*, not EKF sigma_v_*)
    noise = np.concatenate([
        np.atleast_1d(config.sigma_vel_T)   * np.ones(3),
        np.atleast_1d(config.sigma_omega_T) * np.ones(3),
        np.atleast_1d(config.sigma_vel_C)   * np.ones(3),
        np.atleast_1d(config.sigma_omega_C) * np.ones(3),
    ])
    param_noise = np.concatenate([
        np.atleast_1d(config.tau_b)   * np.ones(3),
        np.atleast_1d(config.tau_s)   * np.ones(3),
        np.atleast_1d(config.tau_o)   * np.ones(3),
        np.atleast_1d(config.sigma_b) * np.ones(3),
        np.atleast_1d(config.sigma_s) * np.ones(3),
        np.atleast_1d(config.sigma_o) * np.ones(3),
    ])

    # -- Main Loop -----------------------------
    delta_r = np.zeros(3)  # initialise so status print works at k=0

    t0_wall = time.time()

    for k in tqdm(range(config.N), desc="GNC Sim", unit="step", ncols=80):
        # Sim step time
        t_start, t_end = k * config.SIM_DT, (k+1) * config.SIM_DT
        tspan = [t_start, t_end]

        # -- Guidance --------------------------------------------
        r_des, v_des, q_des, w_des = guidance_law(x_est)

        # -- Control ----------------------------------------------
        w_C = x_true[23:26]
        u_cmd, tau_cmd = control_law(x_est, w_C, r_des, v_des, q_des, w_des, I_c)
        
        # -- Actuator models --------------------------------------
        u_applied = translational_control(u_cmd)
        tau_applied = momentum_wheel_model(tau_cmd)

        if k < 5:
            print(f"\n[k={k}] DIAG")
            print(f"  r_T   = {x_est[0:3]}")
            print(f"  r_C   = {x_est[13:16]}")
            print(f"  r_des = {r_des}")
            print(f"  pos_err_km = {np.linalg.norm(r_des - x_est[13:16]):.6f} km")
            print(f"  u_cmd      = {u_cmd}  km/s²")
            print(f"  u_applied  = {u_applied}  km/s²")   # ← add this
            print(f"  tau_cmd    = {tau_cmd}")

        # -- Simulate Truth Chaser and Target State ---------------
        sol = solve_ivp(
            dynamics_truth, 
            tspan, x_true, method='RK45', 
            args=(u_applied, tau_applied, I_r, I_c, noise, param_noise,), 
            rtol=config.tol, atol=config.tol,
        )
        x_true = sol.y[:, -1]
        x_true[6:10]  = normalize_quat(x_true[6:10])
        x_true[19:23] = normalize_quat(x_true[19:23])

        # -- Sensor Measurements (z) ------------------------------
        t_now = t_end
        omega_gyro = None
        z_st       = None
        z_opt      = None

        if (t_now - last_gyro_meas) >= config.DT_GYRO:
            omega_gyro = measure_gyro(
                x_true[23:26],   # w_C  — true angular velocity
                x_est[23:26],    # b_w  — estimated bias
                x_est[26:29],    # S_s  — estimated scale factor
                x_est[29:32],    # O_o  — estimated misalignment
            )
            last_gyro_meas = t_now

        if (t_now - last_star_track_meas) >= config.DT_STAR_TRACKER:
            z_st = measure_star_tracker(x_true[19:23])
            last_star_track_meas = t_now

        if (t_now - last_opt_cam_meas) >= config.DT_OPT_CAM:
            z_opt = measure_optical_camera(x_true[0:3], x_true[13:16],
                                           x_true[19:23])
            last_opt_cam_meas = t_now
        
        # Add GPS for the chaser

        # -- Navigation -------------------------------------------
        x_est, P = ekf_predict(x_est, P, omega_gyro, u_applied, tau_applied, I_r, I_c, noise, param_noise, config.SIM_DT)

        if z_st is not None:
            x_est, P, _ = ekf_update_star_tracker(x_est, P, z_st, config.R_STAR_TRACKER)
 
        if z_opt is not None:
            x_est, P, _ = ekf_measurement_update(x_est, P, z_opt, config.R_OPT_CAM)

        # GPS chaser update this will improve the estimate substantially
 
        x_est[6:10]  = normalize_quat(x_est[6:10])
        x_est[19:23] = normalize_quat(x_est[19:23])
        
        # -- Compute Errors --------------------------------------
        delta_r = x_true[13:16] - x_true[0:3] - (x_est[13:16] - x_est[0:3])
        delta_q = quat_multiply(x_true[19:23], quat_inv(x_est[19:23]))

        # -- Status Print (every 100 steps) ----------------------
        if k % 100 == 0:
            elapsed = time.time() - t0_wall
            eta     = (elapsed / (k + 1)) * (config.N - k - 1)

            pos_err_m  = np.linalg.norm(x_true[13:16] - x_true[0:3]) * 1e3   # km → m
            vel_err_ms = np.linalg.norm(x_true[16:19] - x_true[3:6]) * 1e3   # km/s → m/s
            nav_err_m  = np.linalg.norm(delta_r) * 1e3                        # km → m

            tqdm.write(
                f"  t={t_end:7.1f}s | "
                f"rel_pos={pos_err_m:8.3f}m | "
                f"rel_vel={vel_err_ms:7.4f}m/s | "
                f"nav_err={nav_err_m:8.3f}m | "
                f"elapsed={elapsed:5.1f}s  ETA={eta:6.1f}s"
            )

        # -- Simulation Variable Storage --------------------------
        hist_target_truth[k, :]  = x_true[0:13]
        hist_chaser_truth[k, :]  = x_true[13:23]
        hist_u_cmd[k, :]         = u_cmd
        hist_tau_cmd[k, :]       = tau_cmd
        hist_u_input[k, :]       = u_applied
        hist_tau_input[k, :]     = tau_applied
        hist_x_est[k, :]         = x_est
        hist_P_diag[k, :]        = np.diag(P)

        if z_opt is not None:
            hist_z_opt[k, :] = z_opt
        if z_st is not None:
            hist_z_st[k, :] = z_st

        # -- Check for Rendezvous Success -------------------------
        pos_err = np.linalg.norm(x_true[13:16] - x_true[0:3])
        att_err = 2 * np.arccos(np.clip(abs(quat_multiply(
                    x_true[19:23], quat_inv(x_true[6:10]))[3]), 0, 1))
 
        if pos_err < config.POS_TOL and att_err < config.ATT_TOL:
            tqdm.write(f"  ✓ Rendezvous achieved at t = {t_end:.1f} s (k={k})")
            break

    elapsed_total = time.time() - t0_wall
    print(f"Simulation complete after {k+1} steps ({(k+1) * config.SIM_DT:.1f} s simulated, {elapsed_total:.1f} s wall time).")

    return {
        'target':    hist_target_truth[:k+1],
        'chaser':    hist_chaser_truth[:k+1],
        'u_cmd':     hist_u_cmd[:k+1],
        'u_applied': hist_u_input[:k+1],
        'x_est':     hist_x_est[:k+1],
        'P_diag':    hist_P_diag[:k+1],
        'z_opt':     hist_z_opt[:k+1],
        'z_st':      hist_z_st[:k+1],
    }


if __name__ == "__main__":
    main()