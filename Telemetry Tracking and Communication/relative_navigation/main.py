import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import config
from scipy.integrate import solve_ivp
from utils.constants import *
from utils.dynamics import dynamics
from utils.orbital_frame_conversions import *
from utils.plotting import *
import matplotlib.pyplot as plt

def main():
    # -- Initialization -------------------------------------------
    ## Initialize the Target from config
    # Translational State
    Target_kep = tle_to_keplerian(config.TARGET_TLE_LINE1, config.TARGET_TLE_LINE2)
    tar_pos, tar_vel = keplerian_to_cartesian(Target_kep)
    target_translational_state = np.concatenate([tar_pos, tar_vel])
    # Attitude State
    target_attitude = np.concatenate([config.TARGET_EP, config.TARGET_ANG_VEL])
    # Full State
    target_state = np.concatenate([target_translational_state, target_attitude])
    I_r = config.TARGET_I

    # Chaser Initialization
    # Translational State
    R_IN = inert_to_RTN_313(Target_kep.raan, Target_kep.i, Target_kep.argp)
    chaser_pos_RTN = R_IN.T @ np.array([0, -0.08, 0])
    print(chaser_pos_RTN)
    chaser_pos = tar_pos + chaser_pos_RTN
    chaser_vel = tar_vel + np.array([0, 0.1, 0.1])  
    chaser_translational_state = np.concatenate([chaser_pos, chaser_vel])
    # Attitude State
    chaser_attitude = np.concatenate([config.CHASER_EP, config.CHASER_ANG_VEL])
    # Full State
    chaser_state = np.concatenate([chaser_translational_state, chaser_attitude])
    I_c = config.CHASER_I

    # -- Truth Propogation ------------------------------------------

    t_eval = np.arange(config.SIM_START, config.SIM_END, config.SIM_DT)

    target_hist = solve_ivp(
        dynamics,
        [config.SIM_START, config.SIM_END],
        target_state,
        method="RK45",
        t_eval=t_eval,
        args=(I_r,),
        rtol=config.tol,
        atol=config.tol,
    )

    chaser_hist = solve_ivp(
        dynamics,
        [config.SIM_START, config.SIM_END],
        chaser_state,
        method="RK45",
        t_eval=t_eval,
        args=(I_c,),
        rtol=config.tol,
        atol=config.tol,
    )

    # sol.y is (13, N) — transpose to (N, 13) for plotting
    plot_all(target_hist.t, target_hist.y.T)
    plot_relative_state(target_hist.t, target_hist.y.T, chaser_hist.y.T)

    plt.show()


    return None

if __name__ == "__main__":
    main()
