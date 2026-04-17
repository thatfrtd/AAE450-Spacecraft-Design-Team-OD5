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
        'c_t': c_t,
        'm_dot': m_dot 
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

def compute_guidance_and_control(w_B, q_B2RTN, r_nozzle_RTN, I_target, centroids, normals, areas, plume_params, ctrl):
    """
    Implements the target detumbling guidance and impingement firing logic.
    Returns the applied torque in the target's Body frame AND a firing flag (1 or 0).
    """
    # Target angular momentum in Body and RTN frames
    h_B = np.dot(I_target, w_B)
    h_norm = np.linalg.norm(h_B)
    
    if h_norm < 1e-6: # Has to be this value bc this is roughly 0.005 deg / s
        return np.zeros(3), 0.0 # Target is already detumbled, Thruster OFF

    # Convert quaternion to Rotation object (Body to RTN mapping)
    attitude = R.from_quat(q_B2RTN) 
    h_L = attitude.apply(h_B)

    # Impingement Torque Guidance 
    T_g_L = -ctrl['K_imp'] * (h_L / np.linalg.norm(h_L))

    # Chaser Pointing Guidance
    r_norm = np.linalg.norm(r_nozzle_RTN)
    r_hat = r_nozzle_RTN / r_norm

    # Project the desired guidance torque onto plane P (orthogonal to r_hat)
    T_g_L_P = T_g_L - np.dot(T_g_L, r_hat) * r_hat
    T_g_L_P_norm = np.linalg.norm(T_g_L_P)

    if T_g_L_P_norm < 1e-6:
        P_h_L = np.zeros(3)
    else:
        T_g_L_P_hat = T_g_L_P / T_g_L_P_norm
        P_h_L = np.cross(T_g_L_P_hat, r_hat)         
        P_h_L = P_h_L / np.linalg.norm(P_h_L)

    # Calculate final thruster Line of Sight (LOS)
    pointing_vector = r_nozzle_RTN - (ctrl['D_imp'] * P_h_L)
    thrust_dir_RTN = -pointing_vector / np.linalg.norm(pointing_vector)

    # Transform vectors back to Body Frame for plume surface physics
    r_nozzle_B = attitude.inv().apply(r_nozzle_RTN)
    thrust_dir_B = attitude.inv().apply(thrust_dir_RTN)

    # Evaluate Impingement Torque
    vis_centroids, vis_normals, vis_areas = get_visible_panels(
        centroids, normals, areas, r_nozzle_B
    )

    T_imp_B = calculate_impingement_torque(
        vis_centroids, vis_normals, vis_areas, r_nozzle_B, thrust_dir_B, plume_params
    )

    T_imp_L = attitude.apply(T_imp_B)
    T_imp_L_norm = np.linalg.norm(T_imp_L)
    
    # 6. PWM Firing Logic Thresholds
    if T_imp_L_norm >= ctrl['eps_m']:
        # Calculate angle between actual plume torque and desired guidance torque
        cos_theta = np.dot(T_imp_L, T_g_L) / (T_imp_L_norm * np.linalg.norm(T_g_L))
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        angle = np.arccos(cos_theta)

        # Only proceed if the torque is pointing in the correct stabilizing direction
        if angle <= ctrl['eps_theta']:
            
            # --- PWM / DUTY CYCLE APPROXIMATION ---
            # Calculate how much torque is actually needed to null the momentum in 1 second
            # T = dH/dt. To kill h_norm in 1 second, we need a torque magnitude of h_norm.
            required_torque = h_norm 
            
            # If the plume provides MORE torque than we need, we fractionally pulse it
            if T_imp_L_norm > required_torque:
                duty_cycle = required_torque / T_imp_L_norm
            else:
                duty_cycle = 1.0 # Fire continuously at 100%

            # Apply the hardware's Minimum Impulse Bit (MIB) limit
            # e.g., if the duty cycle is less than 1% (10ms pulse), the physical valve can't open that fast
            if duty_cycle < 0.01:
                return np.zeros(3), 0.0 

            # Return the scaled torque and the duty cycle for fuel tracking
            return (T_imp_B * duty_cycle), duty_cycle

    # If thresholds are not met, the thruster remains OFF
    return np.zeros(3), 0.0


def optimal_dynamic_diff_equation(t, state, I, I_inv, r_nozzle_RTN, centroids, normals, areas, plume_params, ctrl):
    """
    Computes the state derivative [dq/dt, dw/dt, dm/dt] for the ODE solver.
    State vector: [q_x, q_y, q_z, q_w, w_x, w_y, w_z, m_fuel]
    """
    q = state[0:4]
    w = state[4:7]
    q = q / np.linalg.norm(q)
    
    # Run the Guidance Law to get dynamic torque AND the firing state
    torque_B, is_firing = compute_guidance_and_control(
        w, q, r_nozzle_RTN, I, centroids, normals, areas, plume_params, ctrl
    )
    
    # Euler's rigid body equations
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
    
    # Integrate mass flow rate ONLY when the thruster is firing
    dm_dt = np.array([is_firing * plume_params['m_dot']])
    
    return np.concatenate((dq_dt, dw_dt, dm_dt))
def create_cylinder_geometry():
    """
    Generates the centroids, normals, and areas for a 12-sided 
    cylindrical spacecraft body.
    Dimensions: 25m length, 3m diameter (1.5m radius)
    """
    R = 1.5
    L = 25.0
    N_sides = 12
    N_segments = 25
    
    centroids = []
    normals = []
    areas = []
    
    d_theta = 2.0 * np.pi / N_sides
    dz = L / N_segments
    
    # --- 1. Lateral Faces (The sides of the cylinder) ---
    # The true distance to the center of the flat polygonal face
    R_c = R * np.cos(d_theta / 2.0) 
    
    # Area of one rectangular panel
    face_width = 2.0 * R * np.sin(d_theta / 2.0)
    lateral_area = face_width * dz
    
    # Calculate Z coordinates for the centers of each 1m segment
    z_centers = np.linspace(-L/2 + dz/2, L/2 - dz/2, N_segments)
    
    for z in z_centers:
        for i in range(N_sides):
            theta = i * d_theta
            
            # Normal vector points purely outward radially
            nx = np.cos(theta)
            ny = np.sin(theta)
            nz = 0.0
            normals.append([nx, ny, nz])
            
            # Centroid sits on the face
            cx = R_c * nx
            cy = R_c * ny
            cz = z
            centroids.append([cx, cy, cz])
            
            areas.append(lateral_area)
            
    # --- 2. End Caps (Top and Bottom) ---
    # Divides each circular cap into 12 triangular wedges
    cap_triangle_area = 0.5 * (R**2) * np.sin(d_theta)
    
    for i in range(N_sides):
        # Angular span of the current triangular wedge
        theta_left = i * d_theta - d_theta / 2.0
        theta_right = i * d_theta + d_theta / 2.0
        
        # --- Top Cap (+Z) ---
        V0_top = np.array([0.0, 0.0, L/2]) # Center point
        V1_top = np.array([R * np.cos(theta_left), R * np.sin(theta_left), L/2])
        V2_top = np.array([R * np.cos(theta_right), R * np.sin(theta_right), L/2])
        
        # Centroid of a triangle is the average of its vertices
        c_top = (V0_top + V1_top + V2_top) / 3.0
        centroids.append(c_top.tolist())
        normals.append([0.0, 0.0, 1.0])
        areas.append(cap_triangle_area)
        
        # --- Bottom Cap (-Z) ---
        V0_bot = np.array([0.0, 0.0, -L/2]) # Center point
        V1_bot = np.array([R * np.cos(theta_left), R * np.sin(theta_left), -L/2])
        V2_bot = np.array([R * np.cos(theta_right), R * np.sin(theta_right), -L/2])
        
        c_bottom = (V0_bot + V1_bot + V2_bot) / 3.0
        centroids.append(c_bottom.tolist())
        normals.append([0.0, 0.0, -1.0])
        areas.append(cap_triangle_area)
        
    return np.array(centroids), np.array(normals), np.array(areas)

def extract_aiming_history(sol, I_target, r_nozzle_RTN, ctrl):
    """
    Recalculates the chaser's exact aiming coordinates and Line of Sight  
    over the entire simulated time history.
    """
    num_steps = len(sol.t)
    
    # Init arrays to hold the 3D coordinates over time
    aiming_coords_RTN = np.zeros((3, num_steps))
    thrust_vectors_RTN = np.zeros((3, num_steps))

    for i in range(num_steps):
        q = sol.y[0:4, i]
        w = sol.y[4:7, i]
        q = q / np.linalg.norm(q)

        # Recalculate Guidance Law Parameters
        h_B = np.dot(I_target, w)
        h_norm = np.linalg.norm(h_B)

        # Check deadband 
        if h_norm < 1e-6:
            aiming_coords_RTN[:, i] = np.zeros(3)
            thrust_vectors_RTN[:, i] = np.zeros(3)
            continue

        attitude = R.from_quat(q)
        h_L = attitude.apply(h_B)
        T_g_L = -ctrl['K_imp'] * (h_L / h_norm)

        r_norm = np.linalg.norm(r_nozzle_RTN)
        r_hat = r_nozzle_RTN / r_norm

        T_g_L_P = T_g_L - np.dot(T_g_L, r_hat) * r_hat
        T_g_L_P_norm = np.linalg.norm(T_g_L_P)

        if T_g_L_P_norm < 1e-6:
            P_h_L = np.zeros(3)
        else:
            T_g_L_P_hat = T_g_L_P / T_g_L_P_norm
            P_h_L = np.cross(T_g_L_P_hat, r_hat)
            P_h_L = P_h_L / np.linalg.norm(P_h_L)

        # extract the Specific Aiming Coordinate
        
        aim_coord = ctrl['D_imp'] * P_h_L
        aiming_coords_RTN[:, i] = aim_coord

        # Extract the Unit Thrust Vector
        pointing_vector = r_nozzle_RTN - aim_coord
        thrust_vectors_RTN[:, i] = -pointing_vector / np.linalg.norm(pointing_vector)

    return aiming_coords_RTN, thrust_vectors_RTN

def cart2relative(x, v, a, n):
    '''
    Converts between Cartesian coordinates and Relative orbit dynamics.
    In literature referenced for this project, it is the T' matrix
    - x: Initial position of chaser
    - v: Initial velocity of chaser
    - a: Semi-major axis of chief
    - n: Mean Motion of chief
    '''
    inv_T = np.array([4, 0, 0, 0, 2/n, 0], [0, 1, 0, -2/n, 0, 0], [3*np.cos(n*t), 0, 0, np.sin(n*t)/n, 2*cos(n*t)/n, 0])


############################
#
# Main 
#
############################

## Initialize parameters and geometry
hypergolic_config = init_hypergolic_plume()

# Target inertia matrix 
I_target = np.diag([106400.0, 106400.0, 4500.0])
I_inv = np.linalg.inv(I_target)

# Create 25m Cylinder Geometry
centroids, normals, areas = create_cylinder_geometry()

# Setup initial conditions
initial_quat = np.array([0.0, 0.0, 0.0, 1.0]) 
# MUST BE NON-ZERO to trigger detumbling logic
initial_w = np.radians([0.0, 0.087, 0.037])       
# Append 0.0 for initial fuel consumed
#initial_CW_x = np.array([10, 0, 10, 1, 1, 1])
initial_state = np.concatenate((initial_quat, initial_w, [0.0])) 

# Servicer position (Hovering 14m out, aligned with target CoM for detumbling)
r_nozzle_RTN = np.array([10.0, 0.0, 10]) 

control_params = {
    'D_imp': 10.0,                   # INCREASED: Aim 10m down the cylinder to get a massive lever arm
    'eps_theta': np.radians(90.0),   # RELAXED: Allow firing even if curved surfaces scatter the torque direction up to 90 deg
    'eps_m': 1e-6,                   # RELAXED: Lower minimum threshold to catch smaller torques
    'K_imp': 1e-3                    # Gain of guidance torque 
}

# Setup Time Span for 7200 seconds (2 hours) to allow full momentum decay
t_span = (0.0, 14400.0)
t_eval = np.linspace(t_span[0], t_span[1], 5000)

print("Starting 1-hour numerical integration with ACTIVE Impingement Guidance...")

# Propagate using solve_ivp
sol = solve_ivp(
    optimal_dynamic_diff_equation, 
    t_span, 
    initial_state, 
    t_eval=t_eval, 
    args=(I_target, I_inv, r_nozzle_RTN, centroids, normals, areas, hypergolic_config, control_params),
    rtol=1e-6, 
    atol=1e-8
)

# Extract final state
final_quat = sol.y[0:4, -1]
final_w = sol.y[4:7, -1]
total_fuel_used = sol.y[7, -1] # The 8th state is our fuel tracker

print(f"Integration complete. Success: {sol.success}")
print(f"Final Angular Velocity (deg/s): {np.degrees(final_w)}")
print(f"Total Impingement Fuel Used: {total_fuel_used:.4f} kg")

# --- Plotting ---
time = sol.t
w_history_deg = np.degrees(sol.y[4:7, :])

plt.figure(figsize=(10, 6))
# Plotting normally (no reverse trick) since this is a real detumble
plt.plot(time, w_history_deg[0, :], label=r'$\omega_x$ (Roll)', linewidth=2, color='#1f77b4')
plt.plot(time, w_history_deg[1, :], label=r'$\omega_y$ (Pitch)', linewidth=2, color='#ff7f0e')
plt.plot(time, w_history_deg[2, :], label=r'$\omega_z$ (Yaw)', linewidth=2, color='#2ca02c')

plt.title('Target Debris Active Detumble (Closed-Loop Guidance)', fontsize=14, pad=15)
plt.xlabel('Time [seconds]', fontsize=12)
plt.ylabel('Angular Velocity [deg/s]', fontsize=12)
plt.axhline(0, color='black', linewidth=1, linestyle='--')
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()
plt.show()


## Generate Density Plot
# Spatial Grid
x = np.linspace(1.0, 20.0, 400) 
y = np.linspace(1.0, 20.0, 400)
X, Y = np.meshgrid(x, y)

# 2. Calculate distances (r) and off-axis angles (theta)
r_grid = np.sqrt(X**2 + Y**2)
theta_grid = np.arctan2(np.abs(Y), X)

# 3. Vectorized Density Calculation
r_flat = r_grid.flatten()
theta_flat = theta_grid.flatten()
rho_flat = calculate_plume_density(r_flat, theta_flat, hypergolic_config)
rho_grid = rho_flat.reshape(X.shape)

# 4. Process for Logarithmic Plotting
rho_grid_safe = np.clip(rho_grid, 1e-30, None)
log_rho = np.log10(rho_grid_safe)

# 5. Generate the Figure
plt.figure(figsize=(9, 7)) 

min_density_log = np.min(log_rho)
max_density_log = np.max(log_rho)
levels = np.linspace(min_density_log, max_density_log, 150)

contour = plt.contourf(X, Y, log_rho, levels=levels, cmap='jet')
cbar = plt.colorbar(contour)
cbar.set_label(r'$\log_{10}(\rho)$ $[kg/m^3]$', fontsize=14)

plt.title('Hypergolic Plume Density Field (1m to 10m off-axis)', fontsize=14, pad=15)
plt.xlabel(r'$x_{LOS}$ [m]', fontsize=14)
plt.ylabel(r'$y$ [m]', fontsize=14)
plt.xlim(1, 20)
plt.ylim(1, 20)
plt.tight_layout()
plt.show()


# --- Extract Aiming History ---
aiming_coords, thrust_vectors = extract_aiming_history(sol, I_target, r_nozzle_RTN, control_params)

# --- 3D Visualization of the Servicer's Aiming Strategy ---
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the Target Center of Mass (Origin)
ax.scatter(0, 0, 0, color='black', s=100, label='Target CoM (0,0,0)', marker='x')

# Plot the static Servicer position
ax.scatter(r_nozzle_RTN[0], r_nozzle_RTN[1], r_nozzle_RTN[2], color='red', s=150, label='Servicer Nozzle', marker='^')

# Plot the dynamic aiming coordinates (filtering out the [0,0,0] deadband points)
active_aiming_x = aiming_coords[0, aiming_coords[0,:] != 0]
active_aiming_y = aiming_coords[1, aiming_coords[1,:] != 0]
active_aiming_z = aiming_coords[2, aiming_coords[2,:] != 0]

ax.scatter(active_aiming_x, active_aiming_y, active_aiming_z, color='blue', s=10, alpha=0.3, label='Dynamic Aiming Points')

# Draw a sample Line of Sight from the servicer to the FIRST active aiming point
ax.plot([r_nozzle_RTN[0], active_aiming_x[0]], 
        [r_nozzle_RTN[1], active_aiming_y[0]], 
        [r_nozzle_RTN[2], active_aiming_z[0]], 
        color='red', linestyle='--', alpha=0.7, label='Initial Line of Sight')

# Formatting
ax.set_title('Servicer Plume Impingement Aiming Profile', fontsize=14, pad=15)
ax.set_xlabel('R (Radial) [m]')
ax.set_ylabel('T (Along-Track) [m]')
ax.set_zlabel('N (Cross-Track) [m]')

# Ensure axes are equally scaled so the geometry isn't distorted
max_range = np.array([active_aiming_x.max()-active_aiming_x.min(), 
                      active_aiming_y.max()-active_aiming_y.min(), 
                      active_aiming_z.max()-active_aiming_z.min()]).max() / 2.0
mid_x = (active_aiming_x.max()+active_aiming_x.min()) * 0.5
mid_y = (active_aiming_y.max()+active_aiming_y.min()) * 0.5
mid_z = (active_aiming_z.max()+active_aiming_z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

ax.legend()
plt.show()  