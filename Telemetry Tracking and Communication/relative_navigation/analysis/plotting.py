import numpy as np
import matplotlib.pyplot as plt

def plot_relative_state_overlay(t, target_truth, chaser_truth,
                               target_est, chaser_est):

    # --- Truth relative ---
    rel_truth = chaser_truth[0:3, :] - target_truth[0:3, :]

    # --- Estimated relative ---
    rel_est = chaser_est[0:3, :] - target_est[0:3, :]

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Truth trajectory
    ax.plot3D(rel_truth[0, :], rel_truth[1, :], rel_truth[2, :],
              color='blue', label='Truth')

    # Estimated trajectory
    ax.plot3D(rel_est[0, :], rel_est[1, :], rel_est[2, :],
              color='red', linestyle='--', label='Estimated')

    # Target at origin
    ax.scatter(0, 0, 0, marker='*', s=150, color='gold', label='Target')

    # Start markers
    ax.scatter(*rel_truth[:, 0], color='green', s=60, label='Start (Truth)')
    ax.scatter(*rel_est[:, 0], color='lime',  s=60, label='Start (Est)')

    # End markers
    ax.scatter(*rel_truth[:, -1], color='black', s=60, label='End (Truth)')
    ax.scatter(*rel_est[:, -1], color='red',   s=60, label='End (Est)')

    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_zlabel("z [km]")
    ax.set_title("Relative Trajectory: Truth vs Estimated")

    ax.legend()
    plt.tight_layout()

def plot_relative_state(t: np.ndarray, target_states: np.ndarray, chaser_states: np.ndarray) -> None:

    # Extract positions
    target_r = target_states[0:3, :]
    chaser_r = chaser_states[0:3, :]

    # Relative position
    rel_pos = chaser_r - target_r

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    # Trajectory
    ax.plot3D(rel_pos[0, :], rel_pos[1, :], rel_pos[2, :], label="Relative Trajectory")

    # Target at origin
    ax.scatter(0, 0, 0, marker='*', s=150, color='gold', label="Target")

    # Start point
    ax.scatter(rel_pos[0, 0], rel_pos[1, 0], rel_pos[2, 0],
               color='green', s=60, label="Start")

    # End point
    ax.scatter(rel_pos[0, -1], rel_pos[1, -1], rel_pos[2, -1],
               color='red', s=60, label="End")

    # Labels
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_zlabel("z [km]")
    ax.set_title("Relative Position (Chaser w.r.t Target)")

    ax.legend()
    plt.tight_layout()

def plot_orbits_3d(t: np.ndarray, target_states: np.ndarray, chaser_states: np.ndarray) -> None:
    """
    Plot both target and chaser trajectories in 3D ECI frame with Earth.

    Args:
        t              : (N,) time array, seconds
        target_states  : (13, N) target state history
        chaser_states  : (13, N) chaser state history
    """
    fig = plt.figure(figsize=(9, 9))
    ax  = fig.add_subplot(111, projection="3d")

    # Earth sphere
    R_E = 6371.0
    u   = np.linspace(0, 2 * np.pi, 60)
    v   = np.linspace(0, np.pi, 60)
    xs  = R_E * np.outer(np.cos(u), np.sin(v))
    ys  = R_E * np.outer(np.sin(u), np.sin(v))
    zs  = R_E * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="royalblue", alpha=0.25, linewidth=0)

    # Target trajectory — solid line
    ax.plot(target_states[0, :], target_states[1, :], target_states[2, :],
            color="orangered", linewidth=1.2, label="Target")
    ax.scatter(*target_states[:3,  0], color="green", s=40, zorder=5, label="Target Start")
    ax.scatter(*target_states[:3, -1], color="red",   s=40, zorder=5, label="Target End")

    # Chaser trajectory — dotted markers instead of line
    ax.plot(chaser_states[0, ::20], chaser_states[1, ::20], chaser_states[2, ::20],
            color="gold", marker='+', markersize=6, linestyle='None', label="Chaser")
    ax.scatter(*chaser_states[:3,  0], color="limegreen", s=40, zorder=5, label="Chaser Start")
    ax.scatter(*chaser_states[:3, -1], color="yellow",    s=40, zorder=5, label="Chaser End")
    
    ax.set_xlabel("X [km]")
    ax.set_ylabel("Y [km]")
    ax.set_zlabel("Z [km]")
    ax.set_title("ECI Orbit — Target & Chaser")
    ax.legend()
    plt.tight_layout()

def plot_orbit_3d(t: np.ndarray, states: np.ndarray) -> None:
    """
    Plot spacecraft trajectory in 3D ECI frame with Earth surface sphere.

    Args:
        t      : (N,) time array, seconds  (unused, kept for consistent signature)
        states : (N, 13) state history
    """
    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111, projection="3d")

    # Earth sphere
    R_E = 6371.0   # km
    u   = np.linspace(0, 2 * np.pi, 60)
    v   = np.linspace(0, np.pi, 60)
    xs  = R_E * np.outer(np.cos(u), np.sin(v))
    ys  = R_E * np.outer(np.sin(u), np.sin(v))
    zs  = R_E * np.outer(np.ones_like(u), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="royalblue", alpha=0.3, linewidth=0)

    # Trajectory
    ax.plot(states[:, 0], states[:, 1], states[:, 2],
            color="orangered", linewidth=1.2, label="Trajectory")
    ax.scatter(*states[0,  :3], color="green", s=40, zorder=5, label="Start")
    ax.scatter(*states[-1, :3], color="red",   s=40, zorder=5, label="End")

    ax.set_xlabel("X [km]")
    ax.set_ylabel("Y [km]")
    ax.set_zlabel("Z [km]")
    ax.set_title("ECI Orbit")
    ax.legend()
    plt.tight_layout()
    


def plot_velocity(t: np.ndarray, states: np.ndarray) -> None:
    """
    Plot ECI velocity components over time.

    Args:
        t      : (N,) time array, seconds
        states : (N, 13) state history  — v at indices 3:6
    """
    fig, ax = plt.subplots(figsize=(9, 4))

    labels = ["$v_x$", "$v_y$", "$v_z$"]
    for k, lbl in enumerate(labels):
        ax.plot(t, states[:, 3 + k], label=lbl)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Velocity [km/s]")
    ax.set_title("ECI Velocity")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    


def plot_quaternion(t: np.ndarray, states: np.ndarray) -> None:
    """
    Plot unit quaternion components over time.

    Args:
        t      : (N,) time array, seconds
        states : (N, 13) state history  — q at indices 6:10  [qx, qy, qz, qw]
    """
    fig, ax = plt.subplots(figsize=(9, 4))

    labels = ["$q_x$", "$q_y$", "$q_z$", "$q_w$"]
    for k, lbl in enumerate(labels):
        ax.plot(t, states[:, 6 + k], label=lbl)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Quaternion component")
    ax.set_title("Attitude — Quaternion")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    


def plot_angular_velocity(t: np.ndarray, states: np.ndarray) -> None:
    """
    Plot body-frame angular velocity components over time.

    Args:
        t      : (N,) time array, seconds
        states : (N, 13) state history  — w at indices 10:13
    """
    fig, ax = plt.subplots(figsize=(9, 4))

    labels = ["$\\omega_x$", "$\\omega_y$", "$\\omega_z$"]
    for k, lbl in enumerate(labels):
        ax.plot(t, states[:, 10 + k], label=lbl)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Angular velocity [rad/s]")
    ax.set_title("Attitude — Angular Velocity")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
   


def plot_all(t: np.ndarray, states: np.ndarray) -> None:
    """
    Call all four plotting functions in sequence.

    Args:
        t      : (N,) time vector, seconds
        states : (N, 13) state history
    """
    plot_orbit_3d(t, states)
    plot_velocity(t, states)
    plot_quaternion(t, states)
    plot_angular_velocity(t, states)