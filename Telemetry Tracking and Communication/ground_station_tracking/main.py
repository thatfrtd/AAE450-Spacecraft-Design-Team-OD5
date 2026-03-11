"""
main.py

Entry point for the spacecraft tracking simulation.

Current stage
-------------
    1. Initialise target from TLE, chaser from Keplerian elements
    2. Propagate both truth trajectories over 24 hours
    3. Print summary statistics for both spacecraft
    4. Plot ground track and 3D orbit for both on the same axes

Chaser control law
------------------
    The chaser control law is stubbed as zero thrust.
    Replace the `chaser_control_law` function below with your own when ready.
    Its signature must be:
        def control_law(t: float, state: np.ndarray) -> np.ndarray
            # t     : elapsed seconds since propagation start
            # state : (6,) [rx, ry, rz, vx, vy, vz] ECI, m / m·s⁻¹
            # return: (3,) ECI thrust acceleration, m/s²

Run from the project root:
    python main.py
"""

import os
import numpy as np
import matplotlib.pyplot as plt

import config
from dynamics.spacecraft import Spacecraft
from dynamics.propagator  import Propagator
from analysis.plotting import (
    print_keplerian,
    print_trajectory_stats,
    plot_joint_ground_track,
    plot_joint_orbit_3d,
    plot_joint_keplerian,
)
from sensors.ground_station import load_stations_from_json

os.makedirs(config.OUTPUT_DIR, exist_ok=True)


# ── Chaser control law stub ───────────────────────────────────────────────────

def chaser_control_law(t: float, state: np.ndarray) -> np.ndarray:
    """
    Zero-thrust placeholder.  Replace with your control law when ready.

    Args:
        t     : elapsed seconds since propagation start
        state : (6,) ECI state [r, v], metres / m/s

    Returns:
        (3,) ECI thrust acceleration, m/s²
    """

    return np.zeros(3)


# ── Main ──────────────────────────────────────────────────────────────────────

def run_simulation():

    print("=" * 62)
    print("  Spacecraft Tracking Simulation — Truth Propagation")
    print("=" * 62)
    print(f"  Start : {config.SIM_START.isoformat()}")
    print(f"  End   : {config.SIM_END.isoformat()}")
    print(f"  dt    : {config.TRUTH_DT} s")
    print("=" * 62)

    propagator = Propagator(rtol=config.PROP_RTOL, atol=config.PROP_ATOL)

    # ── 1. Target ─────────────────────────────────────────────────────────────
    print("\n[1/4] Initialising target from TLE ...")

    target = Spacecraft.from_tle(
        tle_line1  = config.TARGET_TLE_LINE1,
        tle_line2  = config.TARGET_TLE_LINE2,
        epoch      = config.SIM_START,
        name       = config.TARGET_NAME,
        propagator = propagator,
    )

    kep_target = target.initial_state.to_keplerian()
    print_keplerian("Target initial elements (osculating, ECI J2000)", kep_target)

    # ── 2. Chaser ─────────────────────────────────────────────────────────────
    print("\n[2/4] Initialising chaser from Keplerian elements ...")

    # Pull RAAN, argp, nu directly from the target's osculating state so the
    # chaser elements match exactly at t=0 — no mean-element ambiguity.
    chaser_kep = np.array([
        config.CHASER_A,      # a    (m)    — user specified
        kep_target[1],        # e           — matched to target
        kep_target[2],        # i    (rad)  — matched to target
        kep_target[3],        # RAAN (rad)  — matched to target
        kep_target[4],        # argp (rad)  — matched to target
        kep_target[5],        # ν    (rad)  — matched to target
    ])

    chaser = Spacecraft.from_keplerian(
        kep        = chaser_kep,
        epoch      = config.SIM_START,
        name       = config.CHASER_NAME,
        control_law= chaser_control_law,
        propagator = propagator,
    )

    print_keplerian("Chaser initial elements", chaser_kep)
    print(f"\n  Note: e / i / RAAN / argp / ν all copied from target osculating state.")
    print(f"  Only semi-major axis differs: "
          f"{(config.CHASER_A - kep_target[0])/1e3:+.4f} km vs target")
    
    # -- Ground Station Locations ---------
    stations = load_stations_from_json("data/ground_stations.json")

    for state in target.truth_trajectory:
        for station in stations:
            obs = station.observe(state)
            if obs is not None:
                obs.range, obs.range_rate, obs.azimuth, obs.elevation


    # ── 3. Propagate ──────────────────────────────────────────────────────────
    print("\n[3/4] Propagating truth trajectories ...")

    target_traj = target.generate_truth_trajectory(
        config.SIM_START, config.SIM_END, config.TRUTH_DT
    )
    print(f"  Target : {len(target_traj)} states")

    chaser_traj = chaser.generate_truth_trajectory(
        config.SIM_START, config.SIM_END, config.TRUTH_DT
    )
    print(f"  Chaser : {len(chaser_traj)} states")

    print_trajectory_stats(config.TARGET_NAME, target_traj)
    print_trajectory_stats(config.CHASER_NAME, chaser_traj)

    # ── 4. Plots ──────────────────────────────────────────────────────────────
    print("\n[4/4] Generating plots ...")

    plot_joint_orbit_3d(
        target_traj, chaser_traj,
        target_name=config.TARGET_NAME,
        chaser_name=config.CHASER_NAME,
        save_path=os.path.join(config.OUTPUT_DIR, "joint_orbit_3d.png"),
    )

    plot_joint_ground_track(
        target_traj, chaser_traj,
        target_name=config.TARGET_NAME,
        chaser_name=config.CHASER_NAME,
        save_path=os.path.join(config.OUTPUT_DIR, "joint_ground_track.png"),
    )

    plot_joint_keplerian(
        target_traj, chaser_traj,
        target_name=config.TARGET_NAME,
        chaser_name=config.CHASER_NAME,
        save_path=os.path.join(config.OUTPUT_DIR, "joint_keplerian.png"),
    )

    print("\nSimulation complete.")
    plt.show()
    return target, chaser, target_traj, chaser_traj


if __name__ == "__main__":
    target, chaser, target_traj, chaser_traj = run_simulation()