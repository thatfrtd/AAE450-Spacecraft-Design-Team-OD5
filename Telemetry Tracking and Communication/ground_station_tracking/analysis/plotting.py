"""
analysis/plotting.py

Visualisation utilities for truth trajectory verification.

Functions
---------
plot_orbit_3d          : 3D ECI orbit with Earth sphere
plot_ground_track      : latitude / longitude ground track
plot_keplerian_elements: six Keplerian elements vs time
plot_cartesian_states  : position and velocity components vs time
"""

from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401 (registers 3d projection)
from typing import List

from dynamics.orbital_state import OrbitalState
from utils.coordinates import (
    RE_EARTH,
    MU_EARTH,
    eci_to_ecef,
    ecef_to_geodetic,
)


# ── Diagnostic utilities ──────────────────────────────────────────────────────

def print_keplerian(label: str, kep: np.ndarray) -> None:
    """Pretty-print a set of Keplerian elements."""
    print(f"\n  {label}:")
    print(f"    a    = {kep[0]/1e3:>12.4f} km")
    print(f"    e    = {kep[1]:>12.6f}")
    print(f"    i    = {np.rad2deg(kep[2]):>12.4f} deg")
    print(f"    RAAN = {np.rad2deg(kep[3]):>12.4f} deg")
    print(f"    ω    = {np.rad2deg(kep[4]):>12.4f} deg")
    print(f"    ν    = {np.rad2deg(kep[5]):>12.4f} deg")


def print_trajectory_stats(name: str, traj: List[OrbitalState]) -> None:
    """Print altitude range, speed range, and energy drift for a trajectory."""
    altitudes = np.array([(s.radius - RE_EARTH) / 1e3 for s in traj])
    speeds    = np.array([s.speed for s in traj])
    energies  = np.array([0.5 * s.speed**2 - MU_EARTH / s.radius for s in traj])

    print(f"\n  {name} trajectory stats:")
    print(f"    Altitude  : min = {altitudes.min():.2f} km, "
          f"max = {altitudes.max():.2f} km")
    print(f"    Speed     : min = {speeds.min():.2f} m/s, "
          f"max = {speeds.max():.2f} m/s")
    print(f"    Δ energy  : {energies[-1] - energies[0]:.4f} J/kg  "
          f"({abs((energies[-1] - energies[0]) / energies[0]) * 100:.4e} %)")


# ── Internal helpers ──────────────────────────────────────────────────────────

def _times_hours(states: List[OrbitalState]) -> np.ndarray:
    """Return elapsed time in hours from the first state's epoch."""
    t0 = states[0].epoch
    return np.array([(s.epoch - t0).total_seconds() / 3600.0 for s in states])


def _geodetic_track(states: List[OrbitalState]) -> tuple:
    """Return (lats_deg, lons_deg) arrays for a list of OrbitalState."""
    lats, lons = [], []
    for s in states:
        r_ecef = eci_to_ecef(s.position, s.epoch)
        lat, lon, _ = ecef_to_geodetic(r_ecef)
        lats.append(np.rad2deg(lat))
        lons.append(np.rad2deg(lon))
    return np.array(lats), np.array(lons)


def _earth_sphere(n: int = 50):
    """Return (x, y, z) arrays for an Earth-radius sphere."""
    u = np.linspace(0, 2 * np.pi, n)
    v = np.linspace(0,     np.pi, n)
    Re = RE_EARTH / 1e3   # km for plotting
    x = Re * np.outer(np.cos(u), np.sin(v))
    y = Re * np.outer(np.sin(u), np.sin(v))
    z = Re * np.outer(np.ones_like(u), np.cos(v))
    return x, y, z


# ── 3D orbit ──────────────────────────────────────────────────────────────────

def plot_orbit_3d(
    states: List[OrbitalState],
    title: str = "Truth Orbit (ECI)",
    save_path: str | None = None,
) -> None:
    """
    Plot the orbit in 3-D ECI coordinates alongside an Earth sphere.

    Args:
        states    : list of OrbitalState from the truth trajectory
        title     : figure title
        save_path : if given, save the figure to this path
    """
    pos_km = np.array([s.position / 1e3 for s in states])   # convert m → km

    fig = plt.figure(figsize=(9, 8))
    ax  = fig.add_subplot(111, projection="3d")

    # Earth
    xs, ys, zs = _earth_sphere()
    ax.plot_surface(xs, ys, zs, color="deepskyblue", alpha=0.35,
                    linewidth=0, zorder=0)

    # Orbit
    ax.plot(pos_km[:, 0], pos_km[:, 1], pos_km[:, 2],
            color="tomato", linewidth=1.0, label="Orbit", zorder=1)

    # Start / end markers
    ax.scatter(*pos_km[0],  color="lime",   s=40, zorder=2, label="Start")
    ax.scatter(*pos_km[-1], color="yellow", s=40, zorder=2, label="End")

    Re = RE_EARTH / 1e3
    ax.set_xlim(-1.5 * Re, 1.5 * Re)
    ax.set_ylim(-1.5 * Re, 1.5 * Re)
    ax.set_zlim(-1.5 * Re, 1.5 * Re)
    ax.set_xlabel("X ECI (km)")
    ax.set_ylabel("Y ECI (km)")
    ax.set_zlabel("Z ECI (km)")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=8)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")



# ── Ground track ──────────────────────────────────────────────────────────────

def plot_ground_track(
    states: List[OrbitalState],
    title: str = "Ground Track",
    save_path: str | None = None,
) -> None:
    """
    Plot geodetic latitude vs longitude ground track.

    Handles the 180° meridian wrap so lines don't cross the plot.

    Args:
        states    : list of OrbitalState from the truth trajectory
        title     : figure title
        save_path : if given, save the figure to this path
    """
    lats, lons = [], []
    for s in states:
        r_ecef = eci_to_ecef(s.position, s.epoch)
        lat, lon, _ = ecef_to_geodetic(r_ecef)
        lats.append(np.rad2deg(lat))
        lons.append(np.rad2deg(lon))

    lats = np.array(lats)
    lons = np.array(lons)

    fig, ax = plt.subplots(figsize=(12, 5))

    # Split at meridian crossings to avoid connecting lines across the map
    breaks = np.where(np.abs(np.diff(lons)) > 180)[0] + 1
    segs   = np.split(np.column_stack([lons, lats]), breaks)

    for k, seg in enumerate(segs):
        ax.plot(seg[:, 0], seg[:, 1],
                color="tomato", linewidth=0.8,
                label="Ground track" if k == 0 else None)

    # Mark start and end
    ax.scatter(lons[0],  lats[0],  color="lime",   s=60, zorder=5,
               label="Start", edgecolors="black", linewidths=0.5)
    ax.scatter(lons[-1], lats[-1], color="yellow", s=60, zorder=5,
               label="End",   edgecolors="black", linewidths=0.5)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 15))
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_title(title)
    ax.legend(fontsize=8)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")



# ── Keplerian elements ────────────────────────────────────────────────────────

def plot_keplerian_elements(
    states: List[OrbitalState],
    title: str = "Keplerian Elements vs Time",
    save_path: str | None = None,
) -> None:
    """
    Plot all six classical Keplerian elements over time.

    Shows J2 secular drift clearly in RAAN, argp, and mean motion.

    Args:
        states    : list of OrbitalState from the truth trajectory
        title     : figure title
        save_path : if given, save the figure to this path
    """
    t_hr = _times_hours(states)
    kep  = np.array([s.to_keplerian() for s in states])

    a    = kep[:, 0] / 1e3           # m → km
    e    = kep[:, 1]
    i    = np.rad2deg(kep[:, 2])
    raan = np.rad2deg(kep[:, 3])
    argp = np.rad2deg(kep[:, 4])
    nu   = np.rad2deg(kep[:, 5])

    labels = ["a (km)", "e", "i (deg)", "RAAN (deg)", "ω (deg)", "ν (deg)"]
    data   = [a, e, i, raan, argp, nu]
    colors = ["steelblue", "darkorange", "mediumseagreen",
              "mediumpurple", "tomato", "saddlebrown"]

    fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    fig.suptitle(title, fontsize=13)

    for ax, d, lbl, col in zip(axes.flat, data, labels, colors):
        ax.plot(t_hr, d, color=col, linewidth=0.9)
        ax.set_ylabel(lbl, fontsize=9)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.yaxis.set_major_formatter(
            plt.FuncFormatter(lambda x, _: f"{x:.4f}" if abs(x) < 1 else f"{x:.2f}")
        )

    for ax in axes[-1]:
        ax.set_xlabel("Time (hours)", fontsize=9)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")



# ── Cartesian state ───────────────────────────────────────────────────────────

def plot_cartesian_states(
    states: List[OrbitalState],
    title: str = "Position & Velocity vs Time",
    save_path: str | None = None,
) -> None:
    """
    Plot ECI position (km) and velocity (m/s) components over time.

    Args:
        states    : list of OrbitalState from the truth trajectory
        title     : figure title
        save_path : if given, save the figure to this path
    """
    t_hr = _times_hours(states)
    pos  = np.array([s.position / 1e3 for s in states])   # m → km
    vel  = np.array([s.velocity       for s in states])    # m/s

    fig  = plt.figure(figsize=(12, 7))
    fig.suptitle(title, fontsize=13)
    gs   = gridspec.GridSpec(2, 1, hspace=0.35)

    # Position
    ax_r = fig.add_subplot(gs[0])
    for k, (comp, col) in enumerate(zip(["X", "Y", "Z"],
                                         ["tomato", "steelblue", "mediumseagreen"])):
        ax_r.plot(t_hr, pos[:, k], color=col, linewidth=0.9, label=f"r_{comp}")
    ax_r.set_ylabel("Position (km)")
    ax_r.legend(fontsize=8, ncol=3)
    ax_r.grid(True, linestyle="--", alpha=0.4)

    # Velocity
    ax_v = fig.add_subplot(gs[1])
    for k, (comp, col) in enumerate(zip(["X", "Y", "Z"],
                                         ["tomato", "steelblue", "mediumseagreen"])):
        ax_v.plot(t_hr, vel[:, k], color=col, linewidth=0.9, label=f"v_{comp}")
    ax_v.set_ylabel("Velocity (m/s)")
    ax_v.set_xlabel("Time (hours)")
    ax_v.legend(fontsize=8, ncol=3)
    ax_v.grid(True, linestyle="--", alpha=0.4)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")


# ── Joint plots (target + chaser) ─────────────────────────────────────────────

def plot_joint_orbit_3d(
    target_traj: List[OrbitalState],
    chaser_traj: List[OrbitalState],
    target_name: str = "Target",
    chaser_name: str = "Chaser",
    save_path: str | None = None,
) -> None:
    """Plot target and chaser orbits in 3D ECI on the same axes."""
    fig = plt.figure(figsize=(9, 8))
    ax  = fig.add_subplot(111, projection="3d")
    fig.suptitle("Truth Orbits — Target & Chaser (ECI)", fontsize=13)

    Re = RE_EARTH / 1e3
    u  = np.linspace(0, 2 * np.pi, 50)
    v  = np.linspace(0,     np.pi, 50)
    ax.plot_surface(
        Re * np.outer(np.cos(u), np.sin(v)),
        Re * np.outer(np.sin(u), np.sin(v)),
        Re * np.outer(np.ones_like(u), np.cos(v)),
        color="deepskyblue", alpha=0.30, linewidth=0, zorder=0,
    )

    for traj, colour, label in [
        (target_traj, "tomato",     target_name),
        (chaser_traj, "dodgerblue", chaser_name),
    ]:
        pos_km = np.array([s.position / 1e3 for s in traj])
        ax.plot(pos_km[:, 0], pos_km[:, 1], pos_km[:, 2],
                color=colour, linewidth=0.9, label=label)
        ax.scatter(*pos_km[0],  color=colour, s=30, marker="o", zorder=2)
        ax.scatter(*pos_km[-1], color=colour, s=30, marker="s", zorder=2)

    lim = 1.5 * Re
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_xlabel("X ECI (km)")
    ax.set_ylabel("Y ECI (km)")
    ax.set_zlabel("Z ECI (km)")
    ax.legend(fontsize=9)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")


def plot_all_targets(
    target1_traj: List[OrbitalState],
    target2_traj: List[OrbitalState],
    target3_traj: List[OrbitalState],
    target1_name: str= "Target 1",
    target2_name: str= "Target 2",
    target3_name: str= "Target 3",
) -> None:
    pass
    

def plot_joint_ground_track(
    target_traj: List[OrbitalState],
    chaser_traj: List[OrbitalState],
    target_name: str = "Target",
    chaser_name: str = "Chaser",
    save_path: str | None = None,
) -> None:
    """Plot target and chaser ground tracks on the same axes."""
    fig, ax = plt.subplots(figsize=(14, 6))
    fig.suptitle("Ground Track — Target & Chaser", fontsize=13)

    for traj, colour, label in [
        (target_traj, "tomato",     target_name),
        (chaser_traj, "dodgerblue", chaser_name),
    ]:
        lats, lons = _geodetic_track(traj)
        breaks = np.where(np.abs(np.diff(lons)) > 180)[0] + 1
        segs   = np.split(np.column_stack([lons, lats]), breaks)
        for k, seg in enumerate(segs):
            ax.plot(seg[:, 0], seg[:, 1], color=colour, linewidth=0.8,
                    label=label if k == 0 else None)
        ax.scatter(lons[0],  lats[0],  color=colour, s=60, zorder=5,
                   edgecolors="black", linewidths=0.5, marker="o")
        ax.scatter(lons[-1], lats[-1], color=colour, s=60, zorder=5,
                   edgecolors="black", linewidths=0.5, marker="s")

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 30))
    ax.set_yticks(range(-90, 91, 15))
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.legend(fontsize=9)
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")


def plot_joint_keplerian(
    target_traj: List[OrbitalState],
    chaser_traj: List[OrbitalState],
    target_name: str = "Target",
    chaser_name: str = "Chaser",
    save_path: str | None = None,
) -> None:
    """Plot Keplerian elements for both spacecraft on the same axes."""
    t_tgt = _times_hours(target_traj)
    t_cha = _times_hours(chaser_traj)
    kep_tgt = np.array([s.to_keplerian() for s in target_traj])
    kep_cha = np.array([s.to_keplerian() for s in chaser_traj])

    labels = ["a (km)", "e", "i (deg)", "RAAN (deg)", "ω (deg)", "ν (deg)"]
    scales = [1/1e3, 1.0, np.rad2deg(1), np.rad2deg(1), np.rad2deg(1), np.rad2deg(1)]

    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    fig.suptitle("Keplerian Elements — Target & Chaser", fontsize=13)

    for ax, col_idx, lbl, sc in zip(axes.flat, range(6), labels, scales):
        ax.plot(t_tgt, kep_tgt[:, col_idx] * sc,
                color="tomato",     linewidth=0.9, label=target_name)
        ax.plot(t_cha, kep_cha[:, col_idx] * sc,
                color="dodgerblue", linewidth=0.9, label=chaser_name,
                linestyle="--")
        ax.set_ylabel(lbl, fontsize=9)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(fontsize=7)

    for ax in axes[-1]:
        ax.set_xlabel("Time (hours)", fontsize=9)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"  Saved → {save_path}")