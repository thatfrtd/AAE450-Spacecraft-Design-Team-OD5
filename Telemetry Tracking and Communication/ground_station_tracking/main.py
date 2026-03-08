"""
main.py

Entry point for the spacecraft tracking simulation.

Simulation stages
-----------------
    1.  Initialise target (TLE) and chaser (Keplerian)
    2.  Load ground stations from JSON
    3.  Propagate truth trajectories for both spacecraft
    4.  Initialise one EKF per spacecraft with a perturbed initial state
    5.  Run the filter loop:
            for each truth timestep
                predict both EKFs forward by dt
                for each station
                    observe target  → update target EKF  (if visible)
                    observe chaser  → update chaser EKF  (if visible)
                log both EKFs
    6.  Print summary statistics
    7.  Plot results

Chaser control law
------------------
    Stubbed as zero thrust. Replace the body of chaser_control_law() with
    your own law when ready. Signature must be:
        def control_law(t: float, state: np.ndarray) -> np.ndarray
            # t     : elapsed seconds since propagation start
            # state : (6,) [rx, ry, rz, vx, vy, vz] ECI, m / m·s⁻¹
            # return: (3,) thrust acceleration in ECI frame, m/s²

Run from project root:
    python main.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt

import config
from dynamics.spacecraft    import Spacecraft
from dynamics.propagator    import Propagator
from dynamics.orbital_state import OrbitalState
from sensors.ground_station import load_stations_from_json
from estimation.ekf         import EKF
from analysis.plotting import (
    print_keplerian,
    print_trajectory_stats,
    plot_joint_ground_track,
    plot_joint_orbit_3d,
    plot_joint_keplerian,
)

os.makedirs(config.OUTPUT_DIR, exist_ok=True)


# ── Chaser control law stub ───────────────────────────────────────────────────

def chaser_control_law(t: float, state: np.ndarray) -> np.ndarray:
    """
    Zero-thrust placeholder. Replace with your real control law when ready.
    """
    return np.zeros(3)


# ── EKF initialisation helper ─────────────────────────────────────────────────

def _init_ekf(true_state: OrbitalState, label: str) -> EKF:
    """
    Create an EKF for one spacecraft.

    The initial state estimate is perturbed from truth by drawing a sample
    from the initial covariance P0. This simulates realistic initial
    uncertainty — the filter should converge from this noisy starting point.

    Args:
        true_state : true OrbitalState at the filter start epoch
        label      : name printed to console

    Returns:
        Initialised EKF instance.
    """
    rng = np.random.default_rng(seed=42)

    # Draw perturbation from N(0, P0)
    L      = np.linalg.cholesky(config.EKF_P0)
    perturb = L @ rng.standard_normal(6)

    init_state = OrbitalState.from_state_vector(
        true_state.state_vector + perturb,
        true_state.epoch,
    )

    ekf = EKF(
        initial_state        = init_state,
        initial_covariance   = config.EKF_P0.copy(),
        process_noise_Q      = config.EKF_Q,
        measurement_noise_R  = config.EKF_R,
        rtol                 = config.PROP_RTOL,
        atol                 = config.PROP_ATOL,
    )

    pos_err = np.linalg.norm(perturb[:3])
    vel_err = np.linalg.norm(perturb[3:])
    print(f"  {label} EKF init error — pos: {pos_err:.1f} m, "
          f"vel: {vel_err:.4f} m/s")
    return ekf


# ── Main ──────────────────────────────────────────────────────────────────────

def run_simulation():

    print("=" * 62)
    print("  Spacecraft Tracking Simulation")
    print("=" * 62)
    print(f"  Start : {config.SIM_START.isoformat()}")
    print(f"  End   : {config.SIM_END.isoformat()}")
    print(f"  dt    : {config.TRUTH_DT} s")
    print("=" * 62)

    propagator = Propagator(rtol=config.PROP_RTOL, atol=config.PROP_ATOL)

    # ── 1. Spacecraft initialisation ──────────────────────────────────────────
    print("\n[1/6] Initialising spacecraft ...")

    target = Spacecraft.from_tle(
        tle_line1  = config.TARGET_TLE_LINE1,
        tle_line2  = config.TARGET_TLE_LINE2,
        epoch      = config.SIM_START,
        name       = config.TARGET_NAME,
        propagator = propagator,
    )
    kep_target = target.initial_state.to_keplerian()
    print_keplerian("Target initial elements (osculating, ECI J2000)", kep_target)

    chaser_kep = np.array([
        config.CHASER_A,
        kep_target[1],
        kep_target[2],
        kep_target[3],
        kep_target[4],
        kep_target[5] + np.deg2rad(30),
    ])
    chaser = Spacecraft.from_keplerian(
        kep         = chaser_kep,
        epoch       = config.SIM_START,
        name        = config.CHASER_NAME,
        control_law = chaser_control_law,
        propagator  = propagator,
    )
    print_keplerian("Chaser initial elements", chaser_kep)
    print(f"\n  Only semi-major axis differs from target: "
          f"{(config.CHASER_A - kep_target[0])/1e3:+.4f} km")

    # ── 2. Ground stations ────────────────────────────────────────────────────
    print("\n[2/6] Loading ground stations ...")
    stations = load_stations_from_json(config.GROUND_STATIONS_JSON)
    for s in stations:
        print(f"  {s}")

    # ── 3. Truth propagation ──────────────────────────────────────────────────
    print("\n[3/6] Propagating truth trajectories ...")

    target_traj = target.generate_truth_trajectory(
        config.SIM_START, config.SIM_END, config.TRUTH_DT
    )
    chaser_traj = chaser.generate_truth_trajectory(
        config.SIM_START, config.SIM_END, config.TRUTH_DT
    )
    print(f"  Target : {len(target_traj)} states")
    print(f"  Chaser : {len(chaser_traj)} states")
    print_trajectory_stats(config.TARGET_NAME, target_traj)
    print_trajectory_stats(config.CHASER_NAME, chaser_traj)

    # ── 4. EKF initialisation ─────────────────────────────────────────────────
    print("\n[4/6] Initialising EKFs ...")

    ekf_target = _init_ekf(target_traj[0], config.TARGET_NAME)
    ekf_chaser = _init_ekf(chaser_traj[0], config.CHASER_NAME)

    # ── 5. Filter loop ────────────────────────────────────────────────────────
    print("\n[5/6] Running EKF loop ...")

    n_steps           = len(target_traj)
    target_obs_counts = np.zeros(n_steps, dtype=int)
    chaser_obs_counts = np.zeros(n_steps, dtype=int)

    for k in range(1, n_steps):

        # Time step between consecutive truth states
        dt = (target_traj[k].epoch - target_traj[k-1].epoch).total_seconds()

        # ── Predict ───────────────────────────────────────────────────────────
        ekf_target.predict(dt)
        ekf_chaser.predict(dt)

        # ── Update — one pass per station ─────────────────────────────────────
        for station in stations:

            # Target observation
            obs_t = station.observe(target_traj[k])
            if obs_t is not None:
                ekf_target.update(
                    observation  = obs_t,
                    station_ecef = station.ecef_position,
                    station_lat  = station.lat_rad,
                    station_lon  = station.lon_rad,
                )
                target_obs_counts[k] += 1

            # Chaser observation
            obs_c = station.observe(chaser_traj[k])
            if obs_c is not None:
                ekf_chaser.update(
                    observation  = obs_c,
                    station_ecef = station.ecef_position,
                    station_lat  = station.lat_rad,
                    station_lon  = station.lon_rad,
                )
                chaser_obs_counts[k] += 1

        # ── Log state after all updates at this timestep ──────────────────────
        ekf_target.log()
        ekf_chaser.log()

    # ── Filter summary ────────────────────────────────────────────────────────
    print(f"\n  Target EKF — total updates : {target_obs_counts.sum()}")
    print(f"  Chaser EKF — total updates : {chaser_obs_counts.sum()}")
    print(f"  Steps with no target obs   : "
          f"{(target_obs_counts[1:] == 0).sum()} / {n_steps - 1}")
    print(f"  Steps with no chaser obs   : "
          f"{(chaser_obs_counts[1:] == 0).sum()} / {n_steps - 1}")

    # Final filter state
    print(f"\n  Final target EKF : {ekf_target}")
    print(f"  Final chaser EKF : {ekf_chaser}")

    # ── 6. Plots ──────────────────────────────────────────────────────────────
    print("\n[6/6] Generating plots ...")

    # Truth orbit and ground track
    plot_joint_orbit_3d(
        target_traj, chaser_traj,
        target_name = config.TARGET_NAME,
        chaser_name = config.CHASER_NAME,
        save_path   = os.path.join(config.OUTPUT_DIR, "joint_orbit_3d.png"),
    )
    plot_joint_ground_track(
        target_traj, chaser_traj,
        target_name = config.TARGET_NAME,
        chaser_name = config.CHASER_NAME,
        save_path   = os.path.join(config.OUTPUT_DIR, "joint_ground_track.png"),
    )
    plot_joint_keplerian(
        target_traj, chaser_traj,
        target_name = config.TARGET_NAME,
        chaser_name = config.CHASER_NAME,
        save_path   = os.path.join(config.OUTPUT_DIR, "joint_keplerian.png"),
    )

    # EKF estimation error plots
    _plot_ekf_errors(
        ekf_target, target_traj,
        label     = config.TARGET_NAME,
        save_path = os.path.join(config.OUTPUT_DIR, "ekf_errors_target.png"),
    )
    _plot_ekf_errors(
        ekf_chaser, chaser_traj,
        label     = config.CHASER_NAME,
        save_path = os.path.join(config.OUTPUT_DIR, "ekf_errors_chaser.png"),
    )

    print("\nSimulation complete.")
    plt.show()

    return target, chaser, target_traj, chaser_traj, stations, ekf_target, ekf_chaser


# ── EKF error plot ────────────────────────────────────────────────────────────

def _plot_ekf_errors(
    ekf        : EKF,
    truth_traj : list,
    label      : str,
    save_path  : str | None = None,
) -> None:
    """
    Plot position and velocity estimation errors with 3-sigma covariance bounds.

    The EKF history and truth trajectory must have the same length.
    If they differ (due to the EKF starting one step later) the shorter
    length is used.
    """
    n = min(len(ekf.state_history), len(truth_traj))

    t0   = truth_traj[0].epoch
    t_hr = np.array([
        (truth_traj[k].epoch - t0).total_seconds() / 3600.0
        for k in range(n)
    ])

    # Position and velocity errors (estimated - truth)
    pos_err = np.array([
        ekf.state_history[k].position - truth_traj[k].position
        for k in range(n)
    ])
    vel_err = np.array([
        ekf.state_history[k].velocity - truth_traj[k].velocity
        for k in range(n)
    ])

    # 3-sigma bounds from covariance diagonal
    sig_pos = np.array([
        3.0 * np.sqrt(np.maximum(np.diag(ekf.covariance_history[k])[:3], 0.0))
        for k in range(n)
    ])
    sig_vel = np.array([
        3.0 * np.sqrt(np.maximum(np.diag(ekf.covariance_history[k])[3:], 0.0))
        for k in range(n)
    ])

    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
    fig.suptitle(f"EKF Estimation Errors — {label}", fontsize=13)

    comp_labels = ["X", "Y", "Z"]
    err_data    = [pos_err, vel_err]
    sig_data    = [sig_pos, sig_vel]
    row_labels  = ["Position error (m)", "Velocity error (m/s)"]
    colors      = ["tomato", "steelblue", "mediumseagreen"]

    for row in range(2):
        for col in range(3):
            ax  = axes[row, col]
            err = err_data[row][:, col]
            sig = sig_data[row][:, col]

            ax.plot(t_hr, err, color=colors[col], linewidth=0.8,
                    label="Error")
            ax.fill_between(t_hr, -sig, sig,
                            color=colors[col], alpha=0.2, label="3σ")
            ax.axhline(0, color="black", linewidth=0.5, linestyle="--")
            ax.set_ylabel(f"{row_labels[row]} {comp_labels[col]}", fontsize=8)
            ax.grid(True, linestyle="--", alpha=0.4)
            if row == 0 and col == 0:
                ax.legend(fontsize=7)

    for ax in axes[-1]:
        ax.set_xlabel("Time (hours)", fontsize=9)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=config.FIGURE_DPI, bbox_inches="tight")
        print(f"  Saved → {save_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    target, chaser, target_traj, chaser_traj, stations, ekf_target, ekf_chaser = (
        run_simulation()
    )