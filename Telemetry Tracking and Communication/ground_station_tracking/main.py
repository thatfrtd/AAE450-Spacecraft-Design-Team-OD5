"""
main.py

Entry point for the spacecraft tracking simulation.

Current stage
-------------
    1. Initialise target from TLE, chaser from Keplerian elements
    2. Load ground stations from JSON
    3. Propagate both truth trajectories over 24 hours
    4. Print summary statistics for both spacecraft
    5. Plot ground track, 3D orbit, and Keplerian elements

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

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import matplotlib.pyplot as plt

import config
from dynamics.spacecraft import Spacecraft
from dynamics.propagator  import Propagator
from sensors.ground_station import load_stations_from_json
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
    Zero-thrust placeholder. Replace with your control law when ready.

    Args:
        t     : elapsed seconds since propagation start
        state : (6,) ECI state [r, v], metres / m·s⁻¹

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

    # All elements matched to the target's osculating state at t=0.
    # Only semi-major axis is user-specified.
    chaser_kep = np.array([
        config.CHASER_A,                       # a    (m)   — user specified
        kep_target[1],                         # e          — matched to target
        kep_target[2],                         # i    (rad) — matched to target
        kep_target[3],                         # RAAN (rad) — matched to target
        kep_target[4],                         # argp (rad) — matched to target
        kep_target[5] + np.deg2rad(30),        # ν    (rad) — offset 30 deg
    ])

    chaser = Spacecraft.from_keplerian(
        kep         = chaser_kep,
        epoch       = config.SIM_START,
        name        = config.CHASER_NAME,
        control_law = chaser_control_law,
        propagator  = propagator,
    )

    print_keplerian("Chaser initial elements", chaser_kep)
    print(f"\n  Note: e / i / RAAN / argp copied from target osculating state.")
    print(f"  Only semi-major axis differs: "
          f"{(config.CHASER_A - kep_target[0])/1e3:+.4f} km vs target")

    # ── 3. Ground stations ────────────────────────────────────────────────────
    print("\n[3/4] Loading ground stations ...")

    stations = load_stations_from_json(config.GROUND_STATIONS_JSON)
    for s in stations:
        print(f"  {s}")

    # ── 4. Propagate ──────────────────────────────────────────────────────────
    print("\n[4/4] Propagating truth trajectories ...")

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

    # ── 5. Plots ──────────────────────────────────────────────────────────────
    print("\n[5/5] Generating plots ...")

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

    print("\nSimulation complete.")
    plt.show()
    return target, chaser, target_traj, chaser_traj, stations


if __name__ == "__main__":
    target, chaser, target_traj, chaser_traj, stations = run_simulation()