import numpy as np
import matplotlib.pyplot as plt


def plot_relative_state(t: np.ndarray, target_states: np.ndarray, chaser_states: np.ndarray) -> None:

    # Find relative position
    target_r = target_states[0:3, :]
    chaser_r = chaser_states[0:3, :]

    rel_pos = chaser_r - target_r

    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    ax.plot3D(rel_pos[0, :], rel_pos[1, :], rel_pos[2, :])
    ax.plot3D(0, 0, 0, marker='*', markersize=15, color='gold')
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_zlabel("z [km]")


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