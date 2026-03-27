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

np.random.seed(0)


def perturb_quat(q, sigma_th):
    """Perturb quaternion multiplicatively — correct MEKF initialisation."""
    dtheta = sigma_th * np.random.randn(3)
    angle  = np.linalg.norm(dtheta)
    if angle < 1e-10:
        dq = np.array([0., 0., 0., 1.])
    else:
        axis = dtheta / angle
        dq   = np.concatenate([np.sin(angle/2)*axis, [np.cos(angle/2)]])
    return quat_multiply(q, dq)


def main():
    results = gnc_sim()

    target_hist = results['target']      # (N, 13)
    chaser_hist = results['chaser']      # (N, 10)
    x_est_hist  = results['x_est']       # (N, 32)
    u_applied   = results['u_applied']   # (N, 3)
    N_sim       = target_hist.shape[0]

    target_states     = target_hist.T
    chaser_states     = np.zeros((13, N_sim))
    chaser_states[0:10, :] = chaser_hist.T
    target_est_states = x_est_hist[:, 0:13].T
    chaser_est_states = x_est_hist[:, 13:23].T
    t = np.arange(N_sim) * config.SIM_DT

    plot_relative_state(t, target_states, chaser_states)
    plot_relative_state_overlay(t, target_states, chaser_states,
                                target_est_states, chaser_est_states)
    plot_orbits_3d(t, target_states, chaser_states)
    plot_chaser_quaternion_truth_vs_est(t, chaser_states, chaser_est_states)
    plot_chaser_omega_truth_vs_est(t, results['omega_true'], results['omega_gyro'])
    plot_target_vs_chaser_quaternion(t, target_states, chaser_states)
    plot_target_quaternion_truth_vs_est(t, target_states, x_est_hist)
    plot_relative_with_arrows(t, target_states, chaser_states,
                              x_est_hist, u_applied,
                              arrow_stride=50, att_scale=0.015)
    plot_estimated_states_with_covariance(t, x_est_hist, results['P_diag'],
                                          target_hist=target_hist,
                                          chaser_hist=chaser_hist)
    """
    animate_rendezvous(
        target_hist, chaser_hist, x_est_hist,
        stride=5,
        target_radius=1.85e-3,
        target_length=8e-3,
        chaser_size=(2e-3, 2e-3, 2e-3),
        axis_scale=8e-3,
        save_path='rendezvous.mp4',
    )
    """

    # Save all figures to plotting/ directory
    save_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plotting')
    os.makedirs(save_dir, exist_ok=True)

    fig_names = [
        'relative_trajectory',
        'relative_truth_vs_est',
        'eci_orbit',
        'chaser_quaternion_truth_vs_est',
        'chaser_gyro_bias',
        'target_vs_chaser_quaternion',
        'target_quaternion_truth_vs_est',
        'relative_with_arrows',
        'target_translational_cov',
        'target_attitude_angvel_cov',
        'chaser_translational_cov',
        'chaser_attitude_gyrobias_cov',
        'nav_parameters_cov',
    ]

    figs = [plt.figure(i) for i in plt.get_fignums()]
    for fig, name in zip(figs, fig_names):
        # fig.savefig(os.path.join(save_dir, f'{name}.png'), dpi=150, bbox_inches='tight')
        print(f"  Saved: {name}.png")

    print(f"\nAll figures saved to: {save_dir}")
    plt.show()


def gnc_sim():
    # -- Initialization ---------------------------------------------------
    target_state = np.concatenate([config.tar_pos, config.tar_vel,
                                   config.TARGET_EP, config.TARGET_ANG_VEL])
    chaser_state = np.concatenate([config.chaser_pos, config.chaser_vel,
                                   config.CHASER_EP, config.CHASER_ANG_VEL])
    I_r        = config.TARGET_I
    I_c        = config.CHASER_I
    param_init = np.concatenate([config.b_w_c, config.ep_S_S, config.ep_O_O])

    # -- History storage --------------------------------------------------
    last_star_track_meas = 0
    last_gyro_meas       = 0
    last_gps_meas        = 0
    last_rel_pos_meas    = 0
    last_target_att_meas = 0

    N = config.N
    hist_target_truth = np.zeros((N, 13))
    hist_chaser_truth = np.zeros((N, 10))
    hist_u_cmd        = np.zeros((N, 3))
    hist_tau_cmd      = np.zeros((N, 3))
    hist_u_input      = np.zeros((N, 3))
    hist_tau_input    = np.zeros((N, 3))
    hist_x_est        = np.zeros((N, 32))
    hist_P_diag       = np.zeros((N, 30))
    hist_z_st         = np.zeros((N, 4))
    hist_z_gps        = np.zeros((N, 3))
    hist_z_rel_pos    = np.zeros((N, 3))
    hist_omega_gyro   = np.zeros((N, 3))
    hist_omega_true   = np.zeros((N, 3))


    # -- Initial states ---------------------------------------------------
    x_true = np.concatenate([target_state, chaser_state, param_init])  # 35

    x_est = np.concatenate([
        target_state + np.concatenate([
            config.sigma_r_T * np.random.randn(3),
            config.sigma_v_T * np.random.randn(3),
            np.zeros(4),               # placeholder — replaced below
            config.sigma_w_T * np.random.randn(3),
        ]),
        chaser_state[:10] + np.concatenate([
            config.sigma_r_C * np.random.randn(3),
            config.sigma_v_C * np.random.randn(3),
            np.zeros(4),               # placeholder — replaced below
        ]),
        param_init + np.concatenate([
            config.sigma_bw * np.random.randn(3),
            config.sigma_Ss * np.random.randn(3),
            config.sigma_Oo * np.random.randn(3),
        ]),
    ])

    # Proper multiplicative quaternion perturbation
    x_est[6:10]  = perturb_quat(target_state[6:10], config.sigma_th_T)
    x_est[19:23] = perturb_quat(chaser_state[6:10], config.sigma_th_C)

    P = np.diag(np.concatenate([
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

    # noise feeds Q matrix only (truth dynamics are deterministic)
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

    delta_r = np.zeros(3)
    t0_wall = time.time()

    for k in tqdm(range(N), desc="GNC Sim", unit="step", ncols=80):
        t_start, t_end = k * config.SIM_DT, (k+1) * config.SIM_DT

        # -- Phase switching: standoff → insertion ────────────────────────
        rho_true = np.linalg.norm(x_true[0:3] - x_true[13:16])
        if rho_true <= config.INSERTION_RANGE_KM:
            config.r_dock_T = config.r_dock_T_insert.copy()
        else:
            config.r_dock_T = config.r_dock_T_standoff.copy()

        # -- Guidance & Control -------------------------------------------
        r_des, v_des, q_des, w_des = guidance_law(x_est)
        w_C        = x_true[23:26]
        u_cmd, tau_cmd = control_law(x_est, w_C, r_des, v_des, q_des, w_des, I_c)
        u_applied      = translational_control(u_cmd)
        tau_applied    = momentum_wheel_model(tau_cmd)

        if k < 5:
            print(f"\n[k={k}] DIAG")
            print(f"  r_T        = {x_est[0:3]}")
            print(f"  r_C        = {x_est[13:16]}")
            print(f"  r_des      = {r_des}")
            print(f"  pos_err_km = {np.linalg.norm(r_des - x_est[13:16]):.6f} km")
            print(f"  u_cmd      = {u_cmd}  km/s²")
            print(f"  u_applied  = {u_applied}  km/s²")
            print(f"  tau_cmd    = {tau_cmd}")

        # -- Truth propagation (fully deterministic) ----------------------
        sol = solve_ivp(
            dynamics_truth, [t_start, t_end], x_true, method='RK45',
            args=(u_applied, tau_applied, I_r, I_c, param_noise,),
            rtol=config.tol, atol=config.tol,
        )
        x_true = sol.y[:, -1]
        x_true[6:10]  = normalize_quat(x_true[6:10])
        x_true[19:23] = normalize_quat(x_true[19:23])

        # -- Sensor measurements ------------------------------------------
        t_now        = t_end
        omega_gyro   = None
        z_st         = None
        z_gps        = None
        z_rel_pos    = None
        z_target_att = None

        if (t_now - last_gyro_meas) >= config.DT_GYRO:
            omega_gyro = measure_gyro(x_true[23:26], x_est[23:26],
                                      x_est[26:29], x_est[29:32])
            last_gyro_meas = t_now

        if (t_now - last_star_track_meas) >= config.DT_STAR_TRACKER:
            z_st = measure_star_tracker(x_true[19:23])
            last_star_track_meas = t_now

        if (t_now - last_gps_meas) >= config.DT_GPS:
            z_gps = measure_gps(x_true[13:16])
            last_gps_meas = t_now

        if (t_now - last_rel_pos_meas) >= config.DT_REL_POS:
            z_rel_pos = measure_relative_position(
                x_true[0:3], x_true[13:16], x_true[19:23])
            last_rel_pos_meas = t_now

        if (t_now - last_target_att_meas) >= config.DT_TARGET_ATT:
            z_target_att = measure_target_attitude(x_true[6:10])
            last_target_att_meas = t_now

        # -- Navigation ---------------------------------------------------
        x_est, P = ekf_predict(x_est, P, omega_gyro, u_applied, tau_applied,
                                I_r, I_c, noise, param_noise, config.SIM_DT)

        if z_st is not None:
            x_est, P, _ = ekf_update_star_tracker(
                x_est, P, z_st, config.R_STAR_TRACKER)

        if z_gps is not None:
            x_est, P, innov_gps = ekf_update_gps(
                x_est, P, z_gps, config.R_GPS)
            if k < 10:
                tqdm.write(
                    f"  [k={k}] GPS innov={np.linalg.norm(innov_gps)*1e3:.3f} m"
                    f"  rC_err={np.linalg.norm(x_true[13:16]-x_est[13:16])*1e3:.3f} m")

        if z_rel_pos is not None:
            rho = np.linalg.norm(x_true[0:3] - x_true[13:16])
            sigma_rel = np.sqrt(config.sigma_lidar**2 +
                                (rho * config.sigma_opt_angle)**2)
            R_rel = (sigma_rel**2) * np.eye(3)
            x_est, P, innov_rel = ekf_update_relative_position(
                x_est, P, z_rel_pos, R_rel)
            if k < 20:
                tqdm.write(
                    f"  [k={k}] REL_POS innov={np.linalg.norm(innov_rel)*1e3:.3f} m"
                    f"  rT_err={np.linalg.norm(x_true[0:3]-x_est[0:3])*1e3:.3f} m")

        if z_target_att is not None:
            x_est, P, _ = ekf_update_target_attitude(
                x_est, P, z_target_att, config.R_TARGET_ATT)

        x_est[6:10]  = normalize_quat(x_est[6:10])
        x_est[19:23] = normalize_quat(x_est[19:23])

        # -- Errors & logging ---------------------------------------------
        delta_r = x_true[13:16] - x_true[0:3] - (x_est[13:16] - x_est[0:3])

        if k % 100 == 0:
            elapsed    = time.time() - t0_wall
            eta        = (elapsed / (k+1)) * (N - k - 1)
            pos_err_m  = np.linalg.norm(x_true[13:16] - x_true[0:3]) * 1e3
            vel_err_ms = np.linalg.norm(x_true[16:19] - x_true[3:6]) * 1e3
            nav_err_m  = np.linalg.norm(delta_r) * 1e3
            phase      = "INSERT" if rho_true <= config.INSERTION_RANGE_KM else "APPROACH"
            tqdm.write(
                f"  t={t_end:7.1f}s | "
                f"rel_pos={pos_err_m:8.3f}m | "
                f"rel_vel={vel_err_ms:7.4f}m/s | "
                f"nav_err={nav_err_m:8.3f}m | "
                f"{phase} | "
                f"elapsed={elapsed:5.1f}s  ETA={eta:6.1f}s"
            )

        hist_target_truth[k, :] = x_true[0:13]
        hist_chaser_truth[k, :] = x_true[13:23]
        hist_u_cmd[k, :]        = u_cmd
        hist_tau_cmd[k, :]      = tau_cmd
        hist_u_input[k, :]      = u_applied
        hist_tau_input[k, :]    = tau_applied
        hist_x_est[k, :]        = x_est
        hist_P_diag[k, :]       = np.diag(P)
        if z_st      is not None: hist_z_st[k, :]      = z_st
        hist_omega_true[k, :]  = x_true[23:26]
        if omega_gyro is not None: hist_omega_gyro[k, :] = omega_gyro
        if z_gps     is not None: hist_z_gps[k, :]     = z_gps
        if z_rel_pos is not None: hist_z_rel_pos[k, :] = z_rel_pos

        pos_err = np.linalg.norm(x_true[13:16] - x_true[0:3])
        att_err = 2 * np.arccos(np.clip(abs(quat_multiply(
                      x_true[19:23], quat_inv(x_est[19:23]))[3]), 0, 1))
        if pos_err < config.POS_TOL and att_err < config.ATT_TOL:
            tqdm.write(f"  ✓ Rendezvous achieved at t = {t_end:.1f} s (k={k})")
            break

    elapsed_total = time.time() - t0_wall
    print(f"Simulation complete after {k+1} steps "
          f"({(k+1)*config.SIM_DT:.1f} s simulated, {elapsed_total:.1f} s wall time).")

    return {
        'target':    hist_target_truth[:k+1],
        'chaser':    hist_chaser_truth[:k+1],
        'u_cmd':     hist_u_cmd[:k+1],
        'u_applied': hist_u_input[:k+1],
        'x_est':     hist_x_est[:k+1],
        'P_diag':    hist_P_diag[:k+1],
        'z_st':      hist_z_st[:k+1],
        'z_gps':     hist_z_gps[:k+1],
        'z_rel_pos':    hist_z_rel_pos[:k+1],
        'omega_gyro':   hist_omega_gyro[:k+1],
        'omega_true':   hist_omega_true[:k+1],
    }


if __name__ == "__main__":
    main()