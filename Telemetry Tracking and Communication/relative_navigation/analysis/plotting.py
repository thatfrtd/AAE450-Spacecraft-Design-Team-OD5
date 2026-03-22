import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# ─────────────────────────────────────────────────────────────────────────────
# Existing plots (unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def plot_relative_state_overlay(t, target_truth, chaser_truth,
                                target_est, chaser_est):
    rel_truth = chaser_truth[0:3, :] - target_truth[0:3, :]
    rel_est   = chaser_est[0:3, :]   - target_est[0:3, :]

    fig = plt.figure()
    ax  = fig.add_subplot(projection="3d")
    ax.plot3D(*rel_truth, color='blue',  label='Truth')
    ax.plot3D(*rel_est,   color='red', linestyle='--', label='Estimated')
    ax.scatter(0, 0, 0, marker='*', s=150, color='gold',  label='Target')
    ax.scatter(*rel_truth[:, 0],  color='green', s=60, label='Start (Truth)')
    ax.scatter(*rel_est[:, 0],    color='lime',  s=60, label='Start (Est)')
    ax.scatter(*rel_truth[:, -1], color='black', s=60, label='End (Truth)')
    ax.scatter(*rel_est[:, -1],   color='red',   s=60, label='End (Est)')
    ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]"); ax.set_zlabel("z [km]")
    ax.set_title("Relative Trajectory: Truth vs Estimated")
    ax.legend(); plt.tight_layout()


def plot_relative_state(t, target_states, chaser_states):
    rel_pos = chaser_states[0:3, :] - target_states[0:3, :]
    fig = plt.figure()
    ax  = fig.add_subplot(projection="3d")
    ax.plot3D(*rel_pos, label="Relative Trajectory")
    ax.scatter(0, 0, 0, marker='*', s=150, color='gold', label="Target")
    ax.scatter(*rel_pos[:, 0],  color='green', s=60, label="Start")
    ax.scatter(*rel_pos[:, -1], color='red',   s=60, label="End")
    ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]"); ax.set_zlabel("z [km]")
    ax.set_title("Relative Position (Chaser w.r.t Target)")
    ax.legend(); plt.tight_layout()


def plot_orbits_3d(t, target_states, chaser_states):
    fig = plt.figure(figsize=(9, 9))
    ax  = fig.add_subplot(111, projection="3d")
    R_E = 6371.0
    u = np.linspace(0, 2*np.pi, 60); v = np.linspace(0, np.pi, 60)
    ax.plot_surface(R_E*np.outer(np.cos(u), np.sin(v)),
                    R_E*np.outer(np.sin(u), np.sin(v)),
                    R_E*np.outer(np.ones_like(u), np.cos(v)),
                    color="royalblue", alpha=0.25, linewidth=0)
    ax.plot(*target_states[:3, :], color="orangered", lw=1.2, label="Target")
    ax.scatter(*target_states[:3,  0], color="green",  s=40, label="Target Start")
    ax.scatter(*target_states[:3, -1], color="red",    s=40, label="Target End")
    ax.plot(chaser_states[0, ::20], chaser_states[1, ::20], chaser_states[2, ::20],
            color="gold", marker='+', markersize=6, linestyle='None', label="Chaser")
    ax.scatter(*chaser_states[:3,  0], color="limegreen", s=40, label="Chaser Start")
    ax.scatter(*chaser_states[:3, -1], color="yellow",    s=40, label="Chaser End")
    ax.set_xlabel("X [km]"); ax.set_ylabel("Y [km]"); ax.set_zlabel("Z [km]")
    ax.set_title("ECI Orbit — Target & Chaser")
    ax.legend(); plt.tight_layout()


def plot_orbit_3d(t, states):
    fig = plt.figure(figsize=(8, 8))
    ax  = fig.add_subplot(111, projection="3d")
    R_E = 6371.0
    u = np.linspace(0, 2*np.pi, 60); v = np.linspace(0, np.pi, 60)
    ax.plot_surface(R_E*np.outer(np.cos(u), np.sin(v)),
                    R_E*np.outer(np.sin(u), np.sin(v)),
                    R_E*np.outer(np.ones_like(u), np.cos(v)),
                    color="royalblue", alpha=0.3, linewidth=0)
    ax.plot(states[:, 0], states[:, 1], states[:, 2],
            color="orangered", lw=1.2, label="Trajectory")
    ax.scatter(*states[0,  :3], color="green", s=40, label="Start")
    ax.scatter(*states[-1, :3], color="red",   s=40, label="End")
    ax.set_xlabel("X [km]"); ax.set_ylabel("Y [km]"); ax.set_zlabel("Z [km]")
    ax.set_title("ECI Orbit"); ax.legend(); plt.tight_layout()


def plot_velocity(t, states):
    fig, ax = plt.subplots(figsize=(9, 4))
    for k, lbl in enumerate(["$v_x$", "$v_y$", "$v_z$"]):
        ax.plot(t, states[:, 3+k], label=lbl)
    ax.set_xlabel("Time [s]"); ax.set_ylabel("Velocity [km/s]")
    ax.set_title("ECI Velocity"); ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()


def plot_quaternion(t, states):
    fig, ax = plt.subplots(figsize=(9, 4))
    for k, lbl in enumerate(["$q_x$", "$q_y$", "$q_z$", "$q_w$"]):
        ax.plot(t, states[:, 6+k], label=lbl)
    ax.set_xlabel("Time [s]"); ax.set_ylabel("Quaternion component")
    ax.set_title("Attitude — Quaternion"); ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()


def plot_angular_velocity(t, states):
    fig, ax = plt.subplots(figsize=(9, 4))
    for k, lbl in enumerate(["$\\omega_x$", "$\\omega_y$", "$\\omega_z$"]):
        ax.plot(t, states[:, 10+k], label=lbl)
    ax.set_xlabel("Time [s]"); ax.set_ylabel("Angular velocity [rad/s]")
    ax.set_title("Attitude — Angular Velocity"); ax.legend(); ax.grid(True, alpha=0.3)
    plt.tight_layout()


def plot_all(t, states):
    plot_orbit_3d(t, states)
    plot_velocity(t, states)
    plot_quaternion(t, states)
    plot_angular_velocity(t, states)


# ─────────────────────────────────────────────────────────────────────────────
# NEW: Chaser quaternion — truth vs estimated
# ─────────────────────────────────────────────────────────────────────────────

def plot_chaser_quaternion_truth_vs_est(t, chaser_truth, chaser_est):
    """
    Compare chaser attitude quaternion truth vs EKF estimate over time.

    Parameters
    ----------
    t            : (N,)   time vector [s]
    chaser_truth : (10,N) truth state rows  — q at rows 6:10
    chaser_est   : (10,N) estimated state rows — q at rows 6:10
    """
    labels  = ["$q_x$", "$q_y$", "$q_z$", "$q_w$"]
    colors  = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Chaser Attitude Quaternion — Truth vs Estimated")

    for i, (lbl, col) in enumerate(zip(labels, colors)):
        ax = axes[i]
        ax.plot(t, chaser_truth[6+i, :], color=col,   lw=1.5, label=f"Truth {lbl}")
        ax.plot(t, chaser_est[6+i, :],   color=col, lw=1.0,
                linestyle='--', alpha=0.8, label=f"Est {lbl}")
        ax.set_ylabel(lbl)
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time [s]")
    plt.tight_layout()


# ─────────────────────────────────────────────────────────────────────────────
# NEW: Chaser angular velocity — truth vs estimated
# ─────────────────────────────────────────────────────────────────────────────

def plot_chaser_omega_truth_vs_est(t, chaser_truth, x_est_hist):
    """
    Compare chaser angular velocity truth vs gyro measurement over time.

    Parameters
    ----------
    t            : (N,)   time vector [s]
    chaser_truth : (10,N) truth chaser state — no omega stored (only 10 states)
                          so truth omega is taken from x_true at [23:26] which
                          is NOT in hist_chaser_truth. Pass None to skip truth line.
    x_est_hist   : (N,32) full estimated state history — gyro bias at [23:26]
                          Note: estimated chaser omega is NOT in x_est (gyro
                          drives attitude directly). We plot the gyro-corrected
                          omega = w_est - b_w_est as a proxy.
    """
    # Estimated gyro bias
    b_w_est = x_est_hist[:, 23:26]        # (N, 3)  bias estimate

    labels = ["$\\omega_x$", "$\\omega_y$", "$\\omega_z$"]
    colors = ["tab:blue", "tab:orange", "tab:green"]

    fig, axes = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
    fig.suptitle("Chaser Angular Velocity — Estimated Gyro Bias")

    for i, (lbl, col) in enumerate(zip(labels, colors)):
        ax = axes[i]
        ax.plot(t, b_w_est[:, i], color=col, lw=1.2, label=f"Bias est {lbl}")
        ax.set_ylabel(lbl + " bias [rad/s]")
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time [s]")
    plt.tight_layout()


# ─────────────────────────────────────────────────────────────────────────────
# NEW: Relative position with thrust vectors and chaser attitude arrows
# ─────────────────────────────────────────────────────────────────────────────

def _quat2dcm(q):
    """q = [qx,qy,qz,qw] → 3×3 DCM  v_body = C @ v_inertial."""
    qx, qy, qz, qw = q
    return np.array([
        [1-2*qy**2-2*qz**2,  2*(qx*qy+qw*qz),   2*(qx*qz-qw*qy)],
        [2*(qx*qy-qw*qz),    1-2*qx**2-2*qz**2,  2*(qy*qz+qw*qx)],
        [2*(qx*qz+qw*qy),    2*(qy*qz-qw*qx),    1-2*qx**2-2*qy**2],
    ])


def plot_relative_with_arrows(t, target_states, chaser_states, x_est_hist,
                               u_applied_hist, arrow_stride=50,
                               thrust_scale=5e3, att_scale=0.015):
    """
    3D relative trajectory with:
      - Blue arrows : applied thrust direction (body frame → inertial)
      - Red arrows  : chaser x-body-axis (pointing direction)

    Parameters
    ----------
    t              : (N,)   time [s]
    target_states  : (13,N) target truth
    chaser_states  : (13,N) chaser truth (10 states padded to 13)
    x_est_hist     : (N,32) estimated state history
    u_applied_hist : (N,3)  applied translational thrust [km/s²]
    arrow_stride   : plot one arrow every this many steps
    thrust_scale   : scale factor for thrust arrows [s²] so km/s² → km
    att_scale      : length of attitude arrows [km]
    """
    rel_pos = chaser_states[0:3, :] - target_states[0:3, :]  # (3,N)

    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection="3d")

    # Trajectory
    ax.plot3D(*rel_pos, color='royalblue', lw=1.2, label='Relative trajectory')
    ax.scatter(0, 0, 0, marker='*', s=200, color='gold', zorder=10, label='Target')
    ax.scatter(*rel_pos[:, 0],  color='green', s=60, zorder=10, label='Start')
    ax.scatter(*rel_pos[:, -1], color='black', s=60, zorder=10, label='End')

    N = rel_pos.shape[1]
    indices = range(0, N, arrow_stride)

    for k in indices:
        origin = rel_pos[:, k]
        q_C    = chaser_states[6:10, k]          # truth chaser quaternion
        u_cmd  = u_applied_hist[k, :]             # applied thrust [km/s²]

        # ── Thrust arrow — normalize so all arrows same visual length
        u_norm = np.linalg.norm(u_cmd)
        if u_norm > 1e-15:
            C_CI   = _quat2dcm(np.array([-q_C[0],-q_C[1],-q_C[2],q_C[3]]))  # body→inertial
            u_dir  = C_CI @ (u_cmd / u_norm)      # unit thrust direction in inertial
            u_rel  = u_dir * att_scale             # same fixed length as attitude arrows
            ax.quiver(*origin, *u_rel,
                      color='crimson', alpha=0.7, linewidth=1.2,
                      arrow_length_ratio=0.3)

        # ── Attitude arrow (chaser x-body-axis in relative frame)
        C_IC  = _quat2dcm(q_C)                   # inertial→body
        x_hat = C_IC[0, :]                        # body x-axis in inertial
        ax.quiver(*origin, *(x_hat * att_scale),
                  color='darkorange', alpha=0.6, linewidth=1.0,
                  arrow_length_ratio=0.4)

    # Legend proxies for arrows
    from matplotlib.lines import Line2D
    legend_extras = [
        Line2D([0],[0], color='crimson',    lw=2, label='Thrust direction'),
        Line2D([0],[0], color='darkorange', lw=2, label='Chaser x-axis (attitude)'),
    ]
    handles, labels_leg = ax.get_legend_handles_labels()
    ax.legend(handles + legend_extras, labels_leg + [e.get_label() for e in legend_extras],
              fontsize=8)

    ax.set_xlabel("x [km]"); ax.set_ylabel("y [km]"); ax.set_zlabel("z [km]")
    ax.set_title("Relative Position with Thrust & Attitude Vectors")
    plt.tight_layout()


# ─────────────────────────────────────────────────────────────────────────────
# NEW: Target quaternion vs Chaser quaternion (truth overlay)
# ─────────────────────────────────────────────────────────────────────────────

def plot_target_vs_chaser_quaternion(t, target_states, chaser_states):
    """
    Overlay target and chaser quaternion truth histories.
    Shows how well the chaser attitude is matching the target.

    Parameters
    ----------
    t             : (N,)   time vector [s]
    target_states : (13,N) target truth — q at rows 6:10
    chaser_states : (13,N) chaser truth (padded) — q at rows 6:10
    """
    labels = ["$q_x$", "$q_y$", "$q_z$", "$q_w$"]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Target vs Chaser Attitude Quaternion (Truth)")

    for i, (lbl, col) in enumerate(zip(labels, colors)):
        ax = axes[i]
        ax.plot(t, target_states[6+i, :], color=col, lw=1.5,
                label=f"Target {lbl}")
        ax.plot(t, chaser_states[6+i, :], color=col, lw=1.0,
                linestyle='--', alpha=0.8, label=f"Chaser {lbl}")
        ax.set_ylabel(lbl)
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time [s]")
    plt.tight_layout()


# ─────────────────────────────────────────────────────────────────────────────
# NEW: Target quaternion — truth vs EKF estimate
# ─────────────────────────────────────────────────────────────────────────────

def plot_target_quaternion_truth_vs_est(t, target_states, x_est_hist):
    """
    Compare target attitude quaternion truth vs EKF estimate over time.

    Parameters
    ----------
    t             : (N,)   time vector [s]
    target_states : (13,N) target truth — q at rows 6:10
    x_est_hist    : (N,32) full estimated state history — target q at cols 6:10
    """
    labels = ["$q_x$", "$q_y$", "$q_z$", "$q_w$"]
    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]

    fig, axes = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
    fig.suptitle("Target Attitude Quaternion — Truth vs Estimated")

    for i, (lbl, col) in enumerate(zip(labels, colors)):
        ax = axes[i]
        ax.plot(t, target_states[6+i, :], color=col, lw=1.5,
                label=f"Truth {lbl}")
        ax.plot(t, x_est_hist[:, 6+i],   color=col, lw=1.0,
                linestyle='--', alpha=0.8, label=f"Est {lbl}")
        ax.set_ylabel(lbl)
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(True, alpha=0.3)

    axes[-1].set_xlabel("Time [s]")
    plt.tight_layout()


# ─────────────────────────────────────────────────────────────────────────────
# 3D Rendezvous Animation — cylinder target, box chaser, principal axes
# ─────────────────────────────────────────────────────────────────────────────

def _make_cylinder(radius, length, n=30):
    """
    Return surface arrays (X,Y,Z) for a cylinder centred at origin,
    axis along body-x, expressed in body frame.
    """
    theta = np.linspace(0, 2*np.pi, n)
    z     = np.linspace(-length/2, length/2, 2)
    theta, z = np.meshgrid(theta, z)
    x = z
    y = radius * np.cos(theta)
    z = radius * np.sin(theta)
    return x, y, z


def _make_box(lx, ly, lz):
    """
    Return a list of (X,Y,Z) face arrays for a box centred at origin
    with half-extents (lx/2, ly/2, lz/2), expressed in body frame.
    """
    hx, hy, hz = lx/2, ly/2, lz/2
    faces = []
    for sign in [-1, 1]:
        # ±X faces
        Y, Z = np.meshgrid([-hy, hy], [-hz, hz])
        X    = sign * hx * np.ones_like(Y)
        faces.append((X, Y, Z))
        # ±Y faces
        X, Z = np.meshgrid([-hx, hx], [-hz, hz])
        Y    = sign * hy * np.ones_like(X)
        faces.append((X, Y, Z))
        # ±Z faces
        X, Y = np.meshgrid([-hx, hx], [-hy, hy])
        Z    = sign * hz * np.ones_like(X)
        faces.append((X, Y, Z))
    return faces


def _transform_surface(X, Y, Z, C_BI, pos):
    """
    Rotate body-frame surface (X,Y,Z) by DCM C_BI (body→inertial)
    and translate to pos.  Returns inertial (X',Y',Z').
    """
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()])  # (3, N)
    pts_I = C_BI @ pts + pos[:, None]
    return (pts_I[0].reshape(X.shape),
            pts_I[1].reshape(X.shape),
            pts_I[2].reshape(X.shape))


def animate_rendezvous(target_hist, chaser_hist, x_est_hist,
                       stride=10,
                       target_radius=1.85e-3, target_length=30e-3,
                       chaser_size=(2e-3, 2e-3, 2e-3),
                       axis_scale=8e-3,
                       save_path=None):
    """
    Animate the rendezvous in the relative (target-centred) frame.

    The target is drawn as a cylinder, the chaser as a rectangular box.
    Three principal-axis arrows are drawn for each spacecraft.

    Parameters
    ----------
    target_hist    : (N,13) target truth history
    chaser_hist    : (N,10) chaser truth history (padded OK)
    x_est_hist     : (N,32) EKF estimate history
    stride         : render every `stride` timesteps (speed vs smoothness)
    target_radius  : cylinder radius [km]
    target_length  : cylinder length [km]
    chaser_size    : (lx,ly,lz) box half-widths [km]
    axis_scale     : length of principal-axis arrows [km]
    save_path      : if given, save to this .mp4 / .gif path instead of showing
    """
    from matplotlib.animation import FuncAnimation
    import matplotlib.colors as mcolors

    # Geometry
    cyl_X, cyl_Y, cyl_Z = _make_cylinder(target_radius, target_length)
    box_faces            = _make_box(*chaser_size)
    axis_colors          = ['red', 'green', 'blue']   # x, y, z body axes
    axis_labels          = ['x', 'y', 'z']

    frames = range(0, min(len(target_hist), len(chaser_hist)), stride)

    fig = plt.figure(figsize=(11, 9), facecolor='#0a0a14')
    ax  = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('#0a0a14')
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
        pane.set_edgecolor('#2a2a3a')

    # Pre-compute relative position range for axis limits
    rel_all = chaser_hist[:, 0:3] - target_hist[:, 0:3]
    rng     = np.max(np.abs(rel_all)) * 1.4
    rng     = max(rng, target_length * 2)

    def _draw_axes(ax, pos, q, scale, label_prefix):
        """Draw three principal-axis quivers at pos with attitude q."""
        C_BI = _quat2dcm(q).T          # body→inertial = C_IB^T = C_BI
        for j, (col, lbl) in enumerate(zip(axis_colors, axis_labels)):
            axis_vec = C_BI[:, j] * scale
            ax.quiver(*pos, *axis_vec,
                      color=col, linewidth=1.5, alpha=0.9,
                      arrow_length_ratio=0.25)

    # ── Initial draw ──────────────────────────────────────────────────────────
    k0      = list(frames)[0]
    r_T0    = target_hist[k0, 0:3]
    r_C0    = chaser_hist[k0, 0:3]
    q_T0    = target_hist[k0, 6:10]
    q_C0    = chaser_hist[k0, 6:10]
    rel0    = r_C0 - r_T0

    C_T0    = _quat2dcm(q_T0).T
    C_C0    = _quat2dcm(q_C0).T

    # Target cylinder surfaces
    cyl_surfs = []
    X, Y, Z = _transform_surface(cyl_X, cyl_Y, cyl_Z, C_T0, np.zeros(3))
    cyl_surfs.append(
        ax.plot_surface(X, Y, Z, color='#c0392b', alpha=0.55, linewidth=0))

    # Chaser box surfaces
    box_surfs = []
    for (bX, bY, bZ) in box_faces:
        X, Y, Z = _transform_surface(bX, bY, bZ, C_C0, rel0)
        box_surfs.append(
            ax.plot_surface(X, Y, Z, color='#2980b9', alpha=0.65, linewidth=0))

    # Trajectory tail
    tail_truth, = ax.plot([], [], [], color='#3498db', lw=0.8, alpha=0.6)
    tail_est,   = ax.plot([], [], [], color='#e67e22', lw=0.8,
                          linestyle='--', alpha=0.5)

    # Time label
    time_text = ax.text2D(0.02, 0.96, '', transform=ax.transAxes,
                          color='white', fontsize=9,
                          fontfamily='monospace')

    ax.set_xlim(-rng, rng); ax.set_ylim(-rng, rng); ax.set_zlim(-rng, rng)
    ax.set_xlabel('x [km]', color='#aaaacc', labelpad=6)
    ax.set_ylabel('y [km]', color='#aaaacc', labelpad=6)
    ax.set_zlabel('z [km]', color='#aaaacc', labelpad=6)
    ax.tick_params(colors='#7777aa', labelsize=7)
    ax.set_title('Rendezvous Animation — Target (red) · Chaser (blue)',
                 color='white', pad=10)

    # Legend proxies
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    legend_els = [
        Patch(facecolor='#c0392b', alpha=0.7, label='Target (cylinder)'),
        Patch(facecolor='#2980b9', alpha=0.8, label='Chaser (box)'),
        Line2D([0],[0], color='red',   lw=1.5, label='Target x-axis'),
        Line2D([0],[0], color='green', lw=1.5, label='Target y-axis'),
        Line2D([0],[0], color='blue',  lw=1.5, label='Target z-axis'),
        Line2D([0],[0], color='#3498db', lw=1, label='Chaser truth path'),
        Line2D([0],[0], color='#e67e22', lw=1,
               linestyle='--', label='Chaser est path'),
    ]
    ax.legend(handles=legend_els, loc='upper right', fontsize=7,
              facecolor='#1a1a2e', edgecolor='#3a3a5e', labelcolor='white')

    frame_list = list(frames)

    def update(fi):
        k   = frame_list[fi]
        r_T = target_hist[k, 0:3]
        r_C = chaser_hist[k, 0:3]
        q_T = target_hist[k, 6:10]
        q_C = chaser_hist[k, 6:10]
        rel = r_C - r_T

        # Estimated relative position
        r_T_est = x_est_hist[k, 0:3]
        r_C_est = x_est_hist[k, 13:16]
        rel_est = r_C_est - r_T_est

        # Remove old surfaces
        for s in cyl_surfs + box_surfs:
            s.remove()
        cyl_surfs.clear(); box_surfs.clear()

        # Remove old quivers (stored as collections added since last frame)
        for coll in list(ax.collections):
            if hasattr(coll, '_is_axis_quiver'):
                coll.remove()

        C_T = _quat2dcm(q_T).T
        C_C = _quat2dcm(q_C).T

        # Redraw target cylinder at origin
        X, Y, Z = _transform_surface(cyl_X, cyl_Y, cyl_Z, C_T, np.zeros(3))
        cyl_surfs.append(
            ax.plot_surface(X, Y, Z, color='#c0392b', alpha=0.55, linewidth=0))

        # Redraw chaser box at relative position
        for (bX, bY, bZ) in box_faces:
            X, Y, Z = _transform_surface(bX, bY, bZ, C_C, rel)
            box_surfs.append(
                ax.plot_surface(X, Y, Z, color='#2980b9', alpha=0.65, linewidth=0))

        # Draw principal axes for target
        for j, col in enumerate(axis_colors):
            ax.quiver(0, 0, 0, *(C_T[:, j] * axis_scale),
                      color=col, lw=1.5, alpha=0.9, arrow_length_ratio=0.25)

        # Draw principal axes for chaser
        for j, col in enumerate(axis_colors):
            ax.quiver(*rel, *(C_C[:, j] * axis_scale * 0.7),
                      color=col, lw=1.0, alpha=0.7, arrow_length_ratio=0.3)

        # Update trajectory tails
        hist_start = max(0, k - 200*stride)
        rel_hist  = chaser_hist[hist_start:k:stride, 0:3] - \
                    target_hist[hist_start:k:stride, 0:3]
        rel_e_hist = x_est_hist[hist_start:k:stride, 13:16] - \
                     x_est_hist[hist_start:k:stride, 0:3]

        if len(rel_hist) > 1:
            tail_truth.set_data(rel_hist[:, 0], rel_hist[:, 1])
            tail_truth.set_3d_properties(rel_hist[:, 2])
        if len(rel_e_hist) > 1:
            tail_est.set_data(rel_e_hist[:, 0], rel_e_hist[:, 1])
            tail_est.set_3d_properties(rel_e_hist[:, 2])

        sep_m = np.linalg.norm(rel) * 1e3
        time_text.set_text(f't = {k:.0f} s  |  sep = {sep_m:.1f} m')

        return cyl_surfs + box_surfs + [tail_truth, tail_est, time_text]

    anim = FuncAnimation(fig, update, frames=len(frame_list),
                         interval=40, blit=False)

    if save_path:
        anim.save(save_path, dpi=120,
                  writer='ffmpeg' if save_path.endswith('.mp4') else 'pillow')
        print(f"Animation saved to {save_path}")
    else:
        plt.show()

    return anim


# ─────────────────────────────────────────────────────────────────────────────
# Estimated states with 3-sigma covariance bounds
# ─────────────────────────────────────────────────────────────────────────────

def plot_estimated_states_with_covariance(t, x_est_hist, P_diag_hist,
                                          target_hist=None, chaser_hist=None):
    """
    Plot all 30 estimated error-states with ±3σ bounds from P diagonal.

    Organized into 7 logical groups, one figure each:
      1. Target position      δr^T   (3)
      2. Target velocity      δv^T   (3)
      3. Target attitude      δθ^T   (3)  + target angular velocity δω^T (3)
      4. Chaser position      δr^C   (3)
      5. Chaser velocity      δv^C   (3)
      6. Chaser attitude      δθ^C   (3)
      7. Parameters           δb_ω, δS_s, δO_o (9)

    If truth histories are provided, the true error (truth − estimate) is
    overlaid so you can see whether the 3σ envelope actually contains the error.

    Parameters
    ----------
    t             : (N,)   time vector [s]
    x_est_hist    : (N,32) estimated state history
    P_diag_hist   : (N,30) diagonal of P at each timestep
    target_hist   : (N,13) optional truth target state (transposed NOT needed)
    chaser_hist   : (N,10) optional truth chaser state
    """

    # ── State and covariance index mapping ───────────────────────────────────
    # x_est cols:   rT(0:3) vT(3:6) qT(6:10) wT(10:13)
    #               rC(13:16) vC(16:19) qC(19:23)
    #               bw(23:26) Ss(26:29) Oo(29:32)
    #
    # P_diag cols:  δrT(0:3) δvT(3:6) δθT(6:9) δwT(9:12)
    #               δrC(12:15) δvC(15:18) δθC(18:21)
    #               δbw(21:24) δSs(24:27) δOo(27:30)

    sigma = np.sqrt(np.abs(P_diag_hist))   # (N,30) — 1σ for each error-state

    def _truth_error(x_est_col, truth_col):
        """Return truth - estimate if truth provided, else None."""
        if truth_col is not None:
            return truth_col - x_est_col
        return None

    def _plot_group(title, groups, fig_size=(12, 8)):
        """
        groups: list of dicts with keys:
          label       : subplot title string
          est_col     : column index into x_est_hist  (or slice)
          p_col       : column index into P_diag / sigma (error-state index)
          truth_vals  : (N,) array of truth values or None
          units       : y-axis unit string
        """
        n = len(groups)
        fig, axes = plt.subplots(n, 1, figsize=fig_size, sharex=True)
        if n == 1:
            axes = [axes]
        fig.suptitle(title, fontsize=11)

        for ax, g in zip(axes, groups):
            est  = x_est_hist[:, g['est_col']]
            sig3 = 3.0 * sigma[:, g['p_col']]

            ax.plot(t, est, color='tab:blue', lw=1.2, label='Estimate')
            ax.fill_between(t, est - sig3, est + sig3,
                            color='tab:blue', alpha=0.18, label='±3σ')
            ax.plot(t, est + sig3, color='tab:blue', lw=0.5, linestyle='--', alpha=0.5)
            ax.plot(t, est - sig3, color='tab:blue', lw=0.5, linestyle='--', alpha=0.5)

            if g['truth_vals'] is not None:
                err = g['truth_vals'] - est
                ax.plot(t, err + est, color='tab:orange', lw=1.0,
                        linestyle=':', label='Truth')

            ax.set_ylabel(f"{g['label']}\n[{g['units']}]", fontsize=8)
            ax.grid(True, alpha=0.25)
            ax.legend(loc='upper right', fontsize=7, ncol=3)

        axes[-1].set_xlabel("Time [s]")
        plt.tight_layout()

    # ── 1. Target Position ───────────────────────────────────────────────────
    _plot_group("Target Position Estimate ± 3σ", [
        dict(label='r^T_x', est_col=0, p_col=0,
             truth_vals=target_hist[:, 0] if target_hist is not None else None,
             units='km'),
        dict(label='r^T_y', est_col=1, p_col=1,
             truth_vals=target_hist[:, 1] if target_hist is not None else None,
             units='km'),
        dict(label='r^T_z', est_col=2, p_col=2,
             truth_vals=target_hist[:, 2] if target_hist is not None else None,
             units='km'),
    ])

    # ── 2. Target Velocity ───────────────────────────────────────────────────
    _plot_group("Target Velocity Estimate ± 3σ", [
        dict(label='v^T_x', est_col=3, p_col=3,
             truth_vals=target_hist[:, 3] if target_hist is not None else None,
             units='km/s'),
        dict(label='v^T_y', est_col=4, p_col=4,
             truth_vals=target_hist[:, 4] if target_hist is not None else None,
             units='km/s'),
        dict(label='v^T_z', est_col=5, p_col=5,
             truth_vals=target_hist[:, 5] if target_hist is not None else None,
             units='km/s'),
    ])

    # ── 3. Target Attitude + Angular Velocity ────────────────────────────────
    # Attitude: quaternion vector part (q stored in x_est at 6:10),
    # but P tracks rotation vector δθ at P_diag cols 6:8
    _plot_group("Target Attitude (δθ) & Angular Velocity ± 3σ", [
        dict(label='δθ^T_x', est_col=6,  p_col=6,
             truth_vals=None, units='rad'),
        dict(label='δθ^T_y', est_col=7,  p_col=7,
             truth_vals=None, units='rad'),
        dict(label='δθ^T_z', est_col=8,  p_col=8,
             truth_vals=None, units='rad'),
        dict(label='ω^T_x',  est_col=10, p_col=9,
             truth_vals=target_hist[:, 10] if target_hist is not None else None,
             units='rad/s'),
        dict(label='ω^T_y',  est_col=11, p_col=10,
             truth_vals=target_hist[:, 11] if target_hist is not None else None,
             units='rad/s'),
        dict(label='ω^T_z',  est_col=12, p_col=11,
             truth_vals=target_hist[:, 12] if target_hist is not None else None,
             units='rad/s'),
    ], fig_size=(12, 11))

    # ── 4. Chaser Position ───────────────────────────────────────────────────
    _plot_group("Chaser Position Estimate ± 3σ", [
        dict(label='r^C_x', est_col=13, p_col=12,
             truth_vals=chaser_hist[:, 0] if chaser_hist is not None else None,
             units='km'),
        dict(label='r^C_y', est_col=14, p_col=13,
             truth_vals=chaser_hist[:, 1] if chaser_hist is not None else None,
             units='km'),
        dict(label='r^C_z', est_col=15, p_col=14,
             truth_vals=chaser_hist[:, 2] if chaser_hist is not None else None,
             units='km'),
    ])

    # ── 5. Chaser Velocity ───────────────────────────────────────────────────
    _plot_group("Chaser Velocity Estimate ± 3σ", [
        dict(label='v^C_x', est_col=16, p_col=15,
             truth_vals=chaser_hist[:, 3] if chaser_hist is not None else None,
             units='km/s'),
        dict(label='v^C_y', est_col=17, p_col=16,
             truth_vals=chaser_hist[:, 4] if chaser_hist is not None else None,
             units='km/s'),
        dict(label='v^C_z', est_col=18, p_col=17,
             truth_vals=chaser_hist[:, 5] if chaser_hist is not None else None,
             units='km/s'),
    ])

    # ── 6. Chaser Attitude ───────────────────────────────────────────────────
    _plot_group("Chaser Attitude (δθ) ± 3σ", [
        dict(label='δθ^C_x', est_col=19, p_col=18,
             truth_vals=None, units='rad'),
        dict(label='δθ^C_y', est_col=20, p_col=19,
             truth_vals=None, units='rad'),
        dict(label='δθ^C_z', est_col=21, p_col=20,
             truth_vals=None, units='rad'),
    ])

    # ── 7. Parameters: gyro bias, scale factor, misalignment ─────────────────
    _plot_group("Navigation Parameters ± 3σ", [
        dict(label='b_ω_x', est_col=23, p_col=21,
             truth_vals=None, units='rad/s'),
        dict(label='b_ω_y', est_col=24, p_col=22,
             truth_vals=None, units='rad/s'),
        dict(label='b_ω_z', est_col=25, p_col=23,
             truth_vals=None, units='rad/s'),
        dict(label='S_s_x', est_col=26, p_col=24,
             truth_vals=None, units='—'),
        dict(label='S_s_y', est_col=27, p_col=25,
             truth_vals=None, units='—'),
        dict(label='S_s_z', est_col=28, p_col=26,
             truth_vals=None, units='—'),
        dict(label='O_o_x', est_col=29, p_col=27,
             truth_vals=None, units='rad'),
        dict(label='O_o_y', est_col=30, p_col=28,
             truth_vals=None, units='rad'),
        dict(label='O_o_z', est_col=31, p_col=29,
             truth_vals=None, units='rad'),
    ], fig_size=(12, 14))