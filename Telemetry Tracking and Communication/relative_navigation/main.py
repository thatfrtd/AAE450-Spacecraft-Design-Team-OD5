import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import config
from scipy.integrate import solve_ivp
from utils.constants import *
from utils.dynamics import *
from utils.plotting import *
import matplotlib.pyplot as plt
from utils.helper import *

def main():
    np.random.seed(0)
    # -- Initialization -------------------------------------------
    # Target Initialization
    # Translational State
    target_translational_state = np.concatenate([config.tar_pos, config.tar_vel])
    # Attitude State
    target_attitude = np.concatenate([config.TARGET_EP, config.TARGET_ANG_VEL])
    # Full State
    target_state = np.concatenate([target_translational_state, target_attitude])
    I_r = config.TARGET_I

    # Chaser Initialization
    # Translational State
    chaser_translational_state = np.concatenate([config.chaser_pos, config.chaser_vel])
    # Attitude State
    chaser_attitude = np.concatenate([config.CHASER_EP, config.CHASER_ANG_VEL])
    # Full State
    chaser_state = np.concatenate([chaser_translational_state, chaser_attitude])
    I_c = config.CHASER_I

    # Parameter dynamics initialization
    param_init = np.concatenate([config.b_w_c, config.ep_S_S, config.ep_O_O])

    # Noise initialization
    t_span = [config.SIM_START, config.SIM_END]
    t_eval = np.arange(config.SIM_START, config.SIM_END, config.SIM_DT)

    noise_g_target     = make_noise_interpolator(t_span, config.SIM_DT, config.sigma_vel_T)
    noise_omega_target = make_noise_interpolator(t_span, config.SIM_DT, config.sigma_omega_T)
    noise_g_chaser     = make_noise_interpolator(t_span, config.SIM_DT, config.sigma_vel_C)
    noise_omega_chaser = make_noise_interpolator(t_span, config.SIM_DT, config.sigma_omega_C)

    # -- Truth Propogation ------------------------------------------

    t_eval = np.arange(config.SIM_START, config.SIM_END, config.SIM_DT)

    target_hist = solve_ivp(
        target_dynamics,
        t_span,
        target_state,
        method="RK45",
        t_eval=t_eval,
        args=(I_r, noise_g_target, noise_omega_target, ),
        rtol=config.tol,
        atol=config.tol,
    )

    chaser_u_dv = np.zeros(3)
    chaser_u_tau = np.zeros(3)
    
    chaser_hist = solve_ivp(
        chaser_dynamics,
        t_span,
        chaser_state,
        method="RK45",
        t_eval=t_eval,
        args=(I_c, noise_g_chaser, noise_omega_chaser, chaser_u_dv, chaser_u_tau, ),
        rtol=config.tol,
        atol=config.tol,
    )

    # paramter dynamics
    """
    param_hist = solve_ivp(
        parameter_dynamics,
        t_span,
        param_init,
        method="RK45",
        args=(config.tau_b, config.tau_s, config.tau_o, config.sigma_b, config.sigma_s, config.sigma_o,  )
    )
    """

    # sol.y is (13, N) — transpose to (N, 13) for plotting
    plot_all(target_hist.t, target_hist.y.T)
    plot_relative_state(target_hist.t, target_hist.y, chaser_hist.y)
    plot_orbits_3d(target_hist.t, target_hist.y, chaser_hist.y)
    plt.show()


    return None

if __name__ == "__main__":
    main()
