import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad, solve_ivp
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt

def create_panels():
    # Creates panels by defining the centroids
    # Right now this is a simple cube
    centroids = np.array([
        # ===== +X FACE (Right) =====
        [ 1.5, -1.0, -1.0], [ 1.5, -1.0,  0.0], [ 1.5, -1.0,  1.0],
        [ 1.5,  0.0, -1.0], [ 1.5,  0.0,  0.0], [ 1.5,  0.0,  1.0],
        [ 1.5,  1.0, -1.0], [ 1.5,  1.0,  0.0], [ 1.5,  1.0,  1.0],
        # ===== -X FACE (Left) =====
        [-1.5, -1.0, -1.0], [-1.5, -1.0,  0.0], [-1.5, -1.0,  1.0],
        [-1.5,  0.0, -1.0], [-1.5,  0.0,  0.0], [-1.5,  0.0,  1.0],
        [-1.5,  1.0, -1.0], [-1.5,  1.0,  0.0], [-1.5,  1.0,  1.0],
        # ===== +Y FACE (Top) =====
        [-1.0,  1.5, -1.0], [-1.0,  1.5,  0.0], [-1.0,  1.5,  1.0],
        [ 0.0,  1.5, -1.0], [ 0.0,  1.5,  0.0], [ 0.0,  1.5,  1.0],
        [ 1.0,  1.5, -1.0], [ 1.0,  1.5,  0.0], [ 1.0,  1.5,  1.0],
        # ===== -Y FACE (Bottom) =====
        [-1.0, -1.5, -1.0], [-1.0, -1.5,  0.0], [-1.0, -1.5,  1.0],
        [ 0.0, -1.5, -1.0], [ 0.0, -1.5,  0.0], [ 0.0, -1.5,  1.0],
        [ 1.0, -1.5, -1.0], [ 1.0, -1.5,  0.0], [ 1.0, -1.5,  1.0],
        # ===== +Z FACE (Front) =====
        [-1.0, -1.0,  1.5], [-1.0,  0.0,  1.5], [-1.0,  1.0,  1.5],
        [ 0.0, -1.0,  1.5], [ 0.0,  0.0,  1.5], [ 0.0,  1.0,  1.5],
        [ 1.0, -1.0,  1.5], [ 1.0,  0.0,  1.5], [ 1.0,  1.0,  1.5],
        # ===== -Z FACE (Back) =====
        [-1.0, -1.0, -1.5], [-1.0,  0.0, -1.5], [-1.0,  1.0, -1.5],
        [ 0.0, -1.0, -1.5], [ 0.0,  0.0, -1.5], [ 0.0,  1.0, -1.5],
        [ 1.0, -1.0, -1.5], [ 1.0,  0.0, -1.5], [ 1.0,  1.0, -1.5]
    ])
    return centroids

def get_visible_panels(centroids, normals, areas, r_nozzle):
    """
    Filters out target panels that are facing away from the servicer's thruster.
    """
    # Calculate view vectors FROM every panel centroid TO the thruster nozzle
    view_vectors = r_nozzle - centroids
    
    # Compute dot product between the view vectors and the panel normals
    dot_products = np.sum(view_vectors * normals, axis=1)
    is_visible = dot_products > 0
    
    # cull the hidden panels
    visible_centroids = centroids[is_visible]
    visible_normals = normals[is_visible]
    visible_areas = areas[is_visible]
    
    return visible_centroids, visible_normals, visible_areas


def calculate_plume_density(r, theta, plume_params):
    """
    Calculates the gas density at specific distances and angles from the nozzle.
    """
    # Extract constants
    rho_star = plume_params['rho_star']
    r_star = plume_params['r_star']
    phi_0 = plume_params['phi_0']
    theta_0 = plume_params['theta_0']
    theta_1 = plume_params['theta_1']
    gamma = plume_params['gamma']
    beta = plume_params['beta']
    
    # pre-allocate the density array
    rho = np.zeros_like(r)
    
    # isentropic core mask (theta <= theta_0)
    core_mask = theta <= theta_0
    
    # boundary layer mask 
    boundary_mask = theta > theta_0
    
    # f(theta) for the isentropic core
    power = 2.0 / (gamma - 1.0)
    f_theta_core = (np.cos(np.pi * theta[core_mask] / (2.0 * theta_1))) ** power
    
    # f(theta) for the Boundary Layer
    f_theta_0 = (np.cos(np.pi * theta_0 / (2.0 * theta_1))) ** power
    f_theta_boundary = f_theta_0 * np.exp(-beta * (theta[boundary_mask] - theta_0))
    
    # f(theta) array
    f_theta = np.zeros_like(theta)
    f_theta[core_mask] = f_theta_core
    f_theta[boundary_mask] = f_theta_boundary
    
    # Apply the Simons Density Equation
    rho = phi_0 * rho_star * ( (r_star / r)**2 ) * f_theta
    
    return rho

def calculate_impingement_torque(vis_centroids, vis_normals, vis_areas, r_nozzle, thrust_dir, plume_params):
    """
    Calculates the total 3D detumbling torque applied to the debris.
    """
    if len(vis_centroids) == 0:
        return np.zeros(3)

    # rel flow vectors 
    # This is L in the paper's force equation
    flow_vectors = vis_centroids - r_nozzle
    r_distances = np.linalg.norm(flow_vectors, axis=1)
    
    # Unit flow vectors (L_hat)
    L_hat = flow_vectors / r_distances[:, np.newaxis]
    
    # Calc angles
    cos_thetas = np.dot(L_hat, thrust_dir)
    # Clip to prevent arccos domain errors from floating point imprecision
    thetas = np.arccos(np.clip(cos_thetas, -1.0, 1.0))
    rho = calculate_plume_density(r_distances, thetas, plume_params)
    
    # 4. Calculate Angle of Incidence (lambda)
    # lambda is the angle between the inward panel normal and the flow vector
    # Because vis_normals are outward, cos(lambda) = vis_normal dot (-L_hat)
    cos_lambdas = np.sum(vis_normals * -L_hat, axis=1)
    
    # Extract force equation constants
    V_lim = plume_params['V_lim']
    V_w = plume_params['V_w']
    c_n = plume_params['c_n']
    c_t = plume_params['c_t']
    
    # Hyper-Thermal Force Equation
    # Scalar multiplier term: rho * V_lim^2 * dA * cos(lambda)
    scalar_mult = rho * (V_lim**2) * vis_areas * cos_lambdas
    
    # Normal bracket term
    normal_bracket = ((2.0 - c_n - c_t) * cos_lambdas + c_n * (V_w / V_lim))[:, np.newaxis]
    
    # Assemble the force vector on each panel
    forces = scalar_mult[:, np.newaxis] * (normal_bracket * vis_normals + c_t * L_hat)
    
    # sum Torques
    torques = np.cross(vis_centroids, forces)
    total_torque = np.sum(torques, axis=0)
    
    return total_torque


def init_hypergolic_plume():
    """
    Initializes the Simons plume model parameters for a 22 N MMH/NTO 
    hypergolic bipropellant thruster. Runs once at simulation startup.
    """
    F = 22.0 # Thrust [N]
    Isp = 285.0 # Specific Impulse [s]
    g0 = 9.80665 # [m/s^2]
    gamma = 1.24  # Sp heat ratio
    R_gas = 378.0  # Sp gas constant [J/(kg K)]
    T_c = 3070.0 # Chamber temp [K]
    P_c = 1.5e6 # Chamber press [Pa]
    E_r = 150.0 # Expansion ratio (A_e / A*)
    delta_Re = 0.02 # Boundary layer ratio (delta / R_E) 
    alpha_e = np.radians(15.0) # Nozzle half-angle [rad]
    
    # Surface interaction constants 
    T_w = 300.0 # Wall temperature [K] 
    c_n = 1.0 # Normal accommodation 
    c_t = 0.97 # Tangential

    # Wall flow-limiting velocity 
    V_w = np.sqrt(2.0 * R_gas * T_w) 
    
    # Gas limiting velocity (Corrected from the paper's typo)
    V_lim = np.sqrt((2.0 * gamma / (gamma - 1.0)) * R_gas * T_c)
    
    # Throat density (rho_star) using choked isentropic flow
    rho_c = P_c / (R_gas * T_c)
    rho_star = rho_c * (2.0 / (gamma + 1.0))**(1.0 / (gamma - 1.0))
    
    # Mass flow rate and Throat Radius (R_star)
    m_dot = F / (Isp * g0)
    T_star = T_c * (2.0 / (gamma + 1.0))
    V_star = np.sqrt(gamma * R_gas * T_star)
    A_star = m_dot / (rho_star * V_star)
    r_star = np.sqrt(A_star / np.pi)

    # Solve for Exit Mach Number (M_e) using the Area-Mach relation
    def area_mach_relation(M):
        term1 = 1.0 / M
        term2 = (2.0 / (gamma + 1.0)) * (1.0 + (gamma - 1.0) / 2.0 * M**2)
        term3 = (gamma + 1.0) / (2.0 * (gamma - 1.0))
        return term1 * (term2 ** term3) - E_r

    M_e = fsolve(area_mach_relation, 4.0)[0] # Initial guess of Mach 4
    
    # Prandtl-Meyer function
    def prandtl_meyer(M):
        t1 = np.sqrt((gamma + 1.0)/(gamma - 1.0))
        t2 = np.arctan(np.sqrt((gamma - 1.0)/(gamma + 1.0) * (M**2 - 1.0)))
        t3 = np.arctan(np.sqrt(M**2 - 1.0))
        return t1 * t2 - t3
    
    nu_e = prandtl_meyer(M_e)
    nu_max = (np.pi / 2.0) * (np.sqrt((gamma + 1.0)/(gamma - 1.0)) - 1.0)
    
    # Limiting turning angle (theta_1) 
    theta_1 = nu_max - nu_e + alpha_e
    
    # Transition angle (theta_0) 
    theta_0 = theta_1 * (1.0 - (2.0 / np.pi) * (delta_Re)**((gamma - 1.0)/(gamma + 1.0)))
    
    # Solve beta and phi
    phi_0 = 1.0 # Initial guess
    tolerance = 1e-6
    error = 1.0
    
    power = 2.0 / (gamma - 1.0)
    f_theta_0 = (np.cos(np.pi * theta_0 / (2.0 * theta_1))) ** power
    
    while error > tolerance:
        # Update beta based on current phi_0 
        beta = phi_0 * np.sqrt((gamma + 1.0)/(gamma - 1.0)) * 1.5 * (delta_Re)**((gamma - 1.0)/(gamma + 1.0))
        
        # Define the angular distribution function f(theta) 
        def f_theta(theta):
            if theta <= theta_0:
                return (np.cos(np.pi * theta / (2.0 * theta_1))) ** power
            else:
                return f_theta_0 * np.exp(-beta * (theta - theta_0))
                
        # Integrand for continuity equation 
        def integrand(theta):
            return f_theta(theta) * np.sin(theta)
            
        # Integrate from 0 to theta_1
        integral_val, _ = quad(integrand, 0, theta_1)
        
        # Update phi_0 
        new_phi_0 = (0.5 * np.sqrt((gamma - 1.0)/(gamma + 1.0))) / integral_val
        
        error = abs(new_phi_0 - phi_0)
        phi_0 = new_phi_0

    plume_params = {
        'rho_star': rho_star,
        'r_star': r_star,
        'phi_0': phi_0,
        'theta_0': theta_0,
        'theta_1': theta_1,
        'gamma': gamma,
        'beta': beta,
        'V_lim': V_lim,
        'V_w': V_w,
        'c_n': c_n,
        'c_t': c_t
    }
    
    return plume_params

def dynamic_diff_equation(t, state, I, I_inv, r_nozzle_RTN, thrust_dir_RTN, centroids, normals, areas, plume_params):
    """
    Computes the state derivative [dq/dt, dw/dt] for the ODE solver.
    State vector: [q_x, q_y, q_z, q_w, w_x, w_y, w_z]
    """
    # Extract quaternion and angular velocity
    q = state[0:4]
    w = state[4:7]
    
    # Normalize quaternion to prevent numerical drift
    q = q / np.linalg.norm(q)
    
    # Create DCM (Rotation object from scipy expects [x, y, z, w])
    attitude = R.from_quat(q)
    
    r_nozzle_B = attitude.inv().apply(r_nozzle_RTN)
    thrust_dir_B = attitude.inv().apply(thrust_dir_RTN)
    
    # Calculate impingement torque in the Body Frame
    vis_centroids, vis_normals, vis_areas = get_visible_panels(
        centroids, normals, areas, r_nozzle_B
    )
    
    torque_B = calculate_impingement_torque(
        vis_centroids, vis_normals, vis_areas, r_nozzle_B, thrust_dir_B, plume_params
    )
    
    # Euler's rigid body
    w_cross_Iw = np.cross(w, np.dot(I, w))
    dw_dt = np.dot(I_inv, (torque_B - w_cross_Iw))
    
    # Quaternion kinematics
    wx, wy, wz = w
    qx, qy, qz, qw = q
    dq_dt = 0.5 * np.array([
         qw*wx + qy*wz - qz*wy,
         qw*wy - qx*wz + qz*wx,
         qw*wz + qx*wy - qy*wx,
        -qx*wx - qy*wy - qz*wz
    ])
    
    return np.concatenate((dq_dt, dw_dt))

############################
#
# Main 
#
############################

## Initialize parameters and initial states
hypergolic_config = init_hypergolic_plume()

# Target inertia matrix from the Borelli paper 
I_target = np.diag([45.0, 25.0, 50.0]) 
I_inv = np.linalg.inv(I_target)

# Servicer position and thrust direction in the global/RTN frame
r_nozzle_RTN = np.array([7.0, 1.0, 1.0])   
thrust_dir_RTN = np.array([-1.0, 0.0, 0.0]) 

areas = np.ones(54)  
centroids = create_panels()
normals = np.array([
    # ===== +X FACE (Right) =====
    [ 1.0,  0.0,  0.0], [ 1.0,  0.0,  0.0], [ 1.0,  0.0,  0.0],
    [ 1.0,  0.0,  0.0], [ 1.0,  0.0,  0.0], [ 1.0,  0.0,  0.0],
    [ 1.0,  0.0,  0.0], [ 1.0,  0.0,  0.0], [ 1.0,  0.0,  0.0],
    # ===== -X FACE (Left) =====
    [-1.0,  0.0,  0.0], [-1.0,  0.0,  0.0], [-1.0,  0.0,  0.0],
    [-1.0,  0.0,  0.0], [-1.0,  0.0,  0.0], [-1.0,  0.0,  0.0],
    [-1.0,  0.0,  0.0], [-1.0,  0.0,  0.0], [-1.0,  0.0,  0.0],
    # ===== +Y FACE (Top) =====
    [ 0.0,  1.0,  0.0], [ 0.0,  1.0,  0.0], [ 0.0,  1.0,  0.0],
    [ 0.0,  1.0,  0.0], [ 0.0,  1.0,  0.0], [ 0.0,  1.0,  0.0],
    [ 0.0,  1.0,  0.0], [ 0.0,  1.0,  0.0], [ 0.0,  1.0,  0.0],
    # ===== -Y FACE (Bottom) =====
    [ 0.0, -1.0,  0.0], [ 0.0, -1.0,  0.0], [ 0.0, -1.0,  0.0],
    [ 0.0, -1.0,  0.0], [ 0.0, -1.0,  0.0], [ 0.0, -1.0,  0.0],
    [ 0.0, -1.0,  0.0], [ 0.0, -1.0,  0.0], [ 0.0, -1.0,  0.0],
    # ===== +Z FACE (Front) =====
    [ 0.0,  0.0,  1.0], [ 0.0,  0.0,  1.0], [ 0.0,  0.0,  1.0],
    [ 0.0,  0.0,  1.0], [ 0.0,  0.0,  1.0], [ 0.0,  0.0,  1.0],
    [ 0.0,  0.0,  1.0], [ 0.0,  0.0,  1.0], [ 0.0,  0.0,  1.0],
    # ===== -Z FACE (Back) =====
    [ 0.0,  0.0, -1.0], [ 0.0,  0.0, -1.0], [ 0.0,  0.0, -1.0],
    [ 0.0,  0.0, -1.0], [ 0.0,  0.0, -1.0], [ 0.0,  0.0, -1.0],
    [ 0.0,  0.0, -1.0], [ 0.0,  0.0, -1.0], [ 0.0,  0.0, -1.0]
])

# Setup initial conditions
initial_quat = np.array([0.0, 0.0, 0.0, 1.0]) # [x, y, z, w]
initial_w = np.radians([0.0, 0.0, 0.0])       # Fast tumble in x-axis 
initial_state = np.concatenate((initial_quat, initial_w))

# Setup Time Span for 5 seconds
t_span = (0.0, 5.0)
t_eval = np.linspace(t_span[0], t_span[1], 100)

print("Starting 5-second numerical integration...")

# Propagate using solve_ivp (RK45 by default)
sol = solve_ivp(
    dynamic_diff_equation, 
    t_span, 
    initial_state, 
    t_eval=t_eval, 
    args=(I_target, I_inv, r_nozzle_RTN, thrust_dir_RTN, centroids, normals, areas, hypergolic_config),
    rtol=1e-6, 
    atol=1e-8
)

# Extract final state
final_quat = sol.y[0:4, -1]
final_w = sol.y[4:7, -1]

print(f"Integration complete. Success: {sol.success}")
print(f"Final Angular Velocity (deg/s): {np.degrees(final_w)}")

time = sol.t
w_history_deg = np.degrees(sol.y[4:7, :])

# Create the plot
w_history_deg_reversed = w_history_deg[:, ::-1]

# Create the plot
plt.figure(figsize=(10, 6))

# Plot the reversed data against the original forward time array
plt.plot(time, w_history_deg_reversed[0, :], label=r'$\omega_x$ (Roll)', linewidth=2, color='#1f77b4')
plt.plot(time, w_history_deg_reversed[1, :], label=r'$\omega_y$ (Pitch)', linewidth=2, color='#ff7f0e')
plt.plot(time, w_history_deg_reversed[2, :], label=r'$\omega_z$ (Yaw)', linewidth=2, color='#2ca02c')
# Formatting
plt.title('Target Debris Angular Velocity Decay (Hypergolic Plume Impingement)', fontsize=14, pad=15)
plt.xlabel('Time [seconds]', fontsize=12)
plt.ylabel('Angular Velocity [deg/s]', fontsize=12)
plt.axhline(0, color='black', linewidth=1, linestyle='--') # Zero reference line
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()

# Display
plt.show()

## Generate Density Plot


# 1. Generate the spatial grid (Restricted to 1-10m in both axes)
x = np.linspace(1.0, 10.0, 200) 
y = np.linspace(1.0, 10.0, 200)
X, Y = np.meshgrid(x, y)

# 2. Calculate distances (r) and off-axis angles (theta)
# r = sqrt(x^2 + y^2)
r_grid = np.sqrt(X**2 + Y**2)

# theta = arctan(|y| / x)
theta_grid = np.arctan2(np.abs(Y), X)

# 3. Vectorized Density Calculation
r_flat = r_grid.flatten()
theta_flat = theta_grid.flatten()

# Pass through our previously built density function
rho_flat = calculate_plume_density(r_flat, theta_flat, hypergolic_config)
rho_grid = rho_flat.reshape(X.shape)

# 4. Process for Logarithmic Plotting
rho_grid_safe = np.clip(rho_grid, 1e-30, None)
log_rho = np.log10(rho_grid_safe)

# 5. Generate the Figure
plt.figure(figsize=(9, 7)) # Adjusted slightly for the square 10x10 aspect ratio

# CRITICAL FIX: Dynamically stretch the color levels between the local min and max
min_density_log = np.min(log_rho)
max_density_log = np.max(log_rho)
levels = np.linspace(min_density_log, max_density_log, 150)

# Generate the filled contour plot
contour = plt.contourf(X, Y, log_rho, levels=levels, cmap='jet')

# Add the colorbar
cbar = plt.colorbar(contour)
cbar.set_label(r'$\log_{10}(\rho)$ $[kg/m^3]$', fontsize=14)

# Formatting
plt.title('Hypergolic Plume Density Field (1m to 10m off-axis)', fontsize=14, pad=15)
plt.xlabel(r'$x_{LOS}$ [m]', fontsize=14)
plt.ylabel(r'$y$ [m]', fontsize=14)

# Lock axis limits explicitly to 1-10
plt.xlim(1, 10)
plt.ylim(1, 10)

plt.tight_layout()
plt.show()