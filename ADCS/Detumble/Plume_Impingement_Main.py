import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import quad, solve_ivp
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt
from numba import njit
import pandas as pd

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

@njit(cache=True)
def fast_get_visible_panels(centroids, normals, areas, r_nozzle):
    """JIT-compiled hidden-surface removal."""
    n = len(centroids)
    mask = np.zeros(n, dtype=np.bool_)
    for i in range(n):
        # View vector from panel to nozzle
        vx = r_nozzle[0] - centroids[i, 0]
        vy = r_nozzle[1] - centroids[i, 1]
        vz = r_nozzle[2] - centroids[i, 2]
        
        # Dot product with outward normal
        dot = vx*normals[i, 0] + vy*normals[i, 1] + vz*normals[i, 2]
        mask[i] = dot > 0
        
    return centroids[mask], normals[mask], areas[mask]

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

@njit(cache=True)
def calculate_impingement_torque_fast(vis_centroids, vis_normals, vis_areas, r_nozzle, thrust_dir, plume_tuple):
    """
    JIT-compiled torque calculation. 
    Memory allocation is dropped to zero inside the loop for maximum speed.
    """
    # Unpack the tuple
# Unpack all 12 elements (added m_dot to the end)
    (rho_star, r_star, phi_0, theta_0, theta_1, gamma, beta, V_lim, V_w, c_n, c_t, m_dot) = plume_tuple    
    num_panels = len(vis_centroids)
    total_torque = np.zeros(3)
    
    if num_panels == 0:
        return total_torque
        
    # Pre-compute static constants for the Simons model
    power = 2.0 / (gamma - 1.0)
    f_theta_0 = (np.cos(np.pi * theta_0 / (2.0 * theta_1))) ** power
    
    for i in range(num_panels):
        # 1. Flow Vector (L) and Distance (r)
        dx = vis_centroids[i, 0] - r_nozzle[0]
        dy = vis_centroids[i, 1] - r_nozzle[1]
        dz = vis_centroids[i, 2] - r_nozzle[2]
        r_dist = np.sqrt(dx**2 + dy**2 + dz**2)
        
        # Unit Flow Vector (L_hat)
        L_x = dx / r_dist
        L_y = dy / r_dist
        L_z = dz / r_dist
        
        # 2. Plume Angle (theta)
        cos_theta = L_x * thrust_dir[0] + L_y * thrust_dir[1] + L_z * thrust_dir[2]
        # Clip to prevent domain errors
        if cos_theta > 1.0: cos_theta = 1.0
        if cos_theta < -1.0: cos_theta = -1.0
        theta = np.arccos(cos_theta)
        
        # 3. Inline Plume Density (rho)
        if theta <= theta_0:
            f_theta = (np.cos(np.pi * theta / (2.0 * theta_1))) ** power
        else:
            f_theta = f_theta_0 * np.exp(-beta * (theta - theta_0))
            
        rho = phi_0 * rho_star * ( (r_star / r_dist)**2 ) * f_theta
        
        # 4. Angle of Incidence (lambda)
        cos_lambda = -(vis_normals[i, 0] * L_x + vis_normals[i, 1] * L_y + vis_normals[i, 2] * L_z)
        
        # 5. Hyper-Thermal Force Calculation
        scalar_mult = rho * (V_lim**2) * vis_areas[i] * cos_lambda
        norm_bracket = (2.0 - c_n - c_t) * cos_lambda + c_n * (V_w / V_lim)
        
        fx = scalar_mult * (norm_bracket * vis_normals[i, 0] + c_t * L_x)
        fy = scalar_mult * (norm_bracket * vis_normals[i, 1] + c_t * L_y)
        fz = scalar_mult * (norm_bracket * vis_normals[i, 2] + c_t * L_z)
        
        # 6. Torque = Cross Product (vis_centroids x forces)
        tx = vis_centroids[i, 1] * fz - vis_centroids[i, 2] * fy
        ty = vis_centroids[i, 2] * fx - vis_centroids[i, 0] * fz
        tz = vis_centroids[i, 0] * fy - vis_centroids[i, 1] * fx
        
        # Accumulate totals
        total_torque[0] += tx
        total_torque[1] += ty
        total_torque[2] += tz
        
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
    vis_centroids, vis_normals, vis_areas = fast_get_visible_panels(
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
def guided_dynamic_diff_equation(t, state, I, I_inv, r_nozzle_RTN, centroids, normals, areas, plume_params, ctrl):
    """
    Computes the state derivative [dq/dt, dw/dt] for the ODE solver.
    State vector: [q_x, q_y, q_z, q_w, w_x, w_y, w_z]
    """
    q = state[0:4]
    w = state[4:7]
    q = q / np.linalg.norm(q)
    
    # Run the Guidance Law to get dynamic torque
    torque_B = compute_guidance_and_control(
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
    
    return np.concatenate((dq_dt, dw_dt))
def compute_guidance_and_control(w_B, q_B2RTN, r_nozzle_RTN, I_target, centroids, normals, areas, plume_params, ctrl):
    """
    Implements the target detumbling guidance and impingement firing logic.
    Returns the applied torque in the target's Body frame.
    """
    # Target angular momentum in Body and RTN frames
    h_B = np.dot(I_target, w_B)
    h_norm = np.linalg.norm(h_B)
    
    if h_norm < 1e-6:
        return np.zeros(3) # Target is already detumbled

    attitude = R.from_quat(q_B2RTN) 
    h_L = attitude.apply(h_B)

    T_g_L = -ctrl['K_imp'] * (h_L / np.linalg.norm(h_L))

    r_norm = np.linalg.norm(r_nozzle_RTN)
    r_hat = r_nozzle_RTN / r_norm

    # Project the desired guidance torque onto plane P (orthogonal to r_hat)
    T_g_L_P = T_g_L - np.dot(T_g_L, r_hat) * r_hat
    T_g_L_P_norm = np.linalg.norm(T_g_L_P)

    if T_g_L_P_norm < 1e-6:
        # Desired torque is parallel to position vector; cannot project.
        P_h_L = np.zeros(3)
    else:
        T_g_L_P_hat = T_g_L_P / T_g_L_P_norm
        # Firing line orthogonal to projection and position vector
        P_h_L = np.cross(r_hat, T_g_L_P_hat)
        P_h_L = P_h_L / np.linalg.norm(P_h_L)

    # Calculate final thruster Line of Sight (LOS) pointing at the offset
    # Vector points from offset target point TO chaser, so thrust is the negative
    pointing_vector = r_nozzle_RTN - (ctrl['D_imp'] * P_h_L)
    thrust_dir_RTN = -pointing_vector / np.linalg.norm(pointing_vector)

    # 4. Transform vectors back to Body Frame for plume surface physics
    r_nozzle_B = attitude.inv().apply(r_nozzle_RTN)
    thrust_dir_B = attitude.inv().apply(thrust_dir_RTN)

    # impingement Torque
    vis_centroids, vis_normals, vis_areas = get_visible_panels(
        centroids, normals, areas, r_nozzle_B
    )

    T_imp_B = calculate_impingement_torque(
        vis_centroids, vis_normals, vis_areas, r_nozzle_B, thrust_dir_B, plume_params
    )

    T_imp_L = attitude.apply(T_imp_B)

    # Firing Logic Thresholds
    T_imp_L_norm = np.linalg.norm(T_imp_L)
    
    # Is it strong enough?
    if T_imp_L_norm >= ctrl['eps_m']:
        # Calculate angle between actual torque and desired guidance torque
        cos_theta = np.dot(T_imp_L, T_g_L) / (T_imp_L_norm * np.linalg.norm(T_g_L))
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        angle = np.arccos(cos_theta)

        # Fire thruster if within angular tolerance
        if angle <= ctrl['eps_theta']:
            return T_imp_B 

def compute_guidance_and_control(w_B, q_B2RTN, r_nozzle_RTN, I_target, centroids, normals, areas, plume_tuple, ctrl):
    # 1. Check Deadband
    h_B = np.dot(I_target, w_B)
    h_norm = np.linalg.norm(h_B)
    
    if h_norm < 0.01: # Realistic 0.01 Nms deadband
        return np.zeros(3), 0.0

    # 2. Fast Coordinate Transforms
    R_B2RTN = quat_to_dcm(q_B2RTN)
    R_RTN2B = R_B2RTN.T 
    
    h_L = np.dot(R_B2RTN, h_B)

    # 3. Guidance Law (T_g_L)
    T_g_L = -ctrl['K_imp'] * (h_L / h_norm)

    r_norm = np.linalg.norm(r_nozzle_RTN)
    r_hat = r_nozzle_RTN / r_norm

    T_g_L_P = T_g_L - np.dot(T_g_L, r_hat) * r_hat
    T_g_L_P_norm = np.linalg.norm(T_g_L_P)

    if T_g_L_P_norm < 1e-6:
        P_h_L = np.zeros(3)
    else:
        P_h_L = np.cross((T_g_L_P / T_g_L_P_norm), r_hat)        
        P_h_L = P_h_L / np.linalg.norm(P_h_L)

    # Calculate final LOS pointing
    pointing_vector = r_nozzle_RTN - (ctrl['D_imp'] * P_h_L)
    thrust_dir_RTN = -pointing_vector / np.linalg.norm(pointing_vector)

    # 4. Transform to Body Frame
    r_nozzle_B = np.dot(R_RTN2B, r_nozzle_RTN)
    thrust_dir_B = np.dot(R_RTN2B, thrust_dir_RTN)

    # 5. Fast Plume Physics
    vis_centroids, vis_normals, vis_areas = fast_get_visible_panels(
        centroids, normals, areas, r_nozzle_B
    )

    T_imp_B = calculate_impingement_torque_fast(
        vis_centroids, vis_normals, vis_areas, r_nozzle_B, thrust_dir_B, plume_tuple
    )

    # 6. PWM Firing Logic
    T_imp_L = np.dot(R_B2RTN, T_imp_B)
    T_imp_L_norm = np.linalg.norm(T_imp_L)
    
    if T_imp_L_norm >= ctrl['eps_m']:
        cos_theta = np.dot(T_imp_L, T_g_L) / (T_imp_L_norm * np.linalg.norm(T_g_L))
        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        angle = np.arccos(cos_theta)

        if angle <= ctrl['eps_theta']:
            required_torque = h_norm 
            duty_cycle = required_torque / T_imp_L_norm if T_imp_L_norm > required_torque else 1.0

            if duty_cycle < 0.01:
                return np.zeros(3), 0.0 

            return (T_imp_B * duty_cycle), duty_cycle

    return np.zeros(3), 0.0


def optimal_dynamic_diff_equation(t, state, I, I_inv, centroids, normals, areas, plume_tuple, ctrl, omega):
    """14-Element State: [q_x, q_y, q_z, q_w, w_x, w_y, w_z, m_fuel, x, y, z, vx, vy, vz]"""
    
    # Extract Target Attitude & Spin
    q = state[0:4]
    w = state[4:7]
    
    # Extract Chaser Relative State
    r_nozzle_RTN = state[8:11] 
    v_chaser_RTN = state[11:14]
    
    # Run Guidance Law
    torque_B, duty_cycle = compute_guidance_and_control(
        w, q, r_nozzle_RTN, I, centroids, normals, areas, plume_tuple, ctrl
    )
    
    # Target Rigid Body Dynamics (Euler)
    w_cross_Iw = np.cross(w, np.dot(I, w))
    dw_dt = np.dot(I_inv, (torque_B - w_cross_Iw))
    
    # Target Kinematics (Quaternions - NO internal normalization!)
    wx, wy, wz = w
    qx, qy, qz, qw = q
    dq_dt = 0.5 * np.array([
         qw*wx + qy*wz - qz*wy,
         qw*wy - qx*wz + qz*wx,
         qw*wz + qx*wy - qy*wx,
        -qx*wx - qy*wy - qz*wz
    ])
    
    # Fuel Consumption (Index 11 in tuple is m_dot)
    dm_dt = np.array([duty_cycle * plume_tuple[11]])
    
    # Chaser Orbital Dynamics (CW Equations)
    n = omega
    ax = 3 * n**2 * r_nozzle_RTN[0] + 2 * n * v_chaser_RTN[1]
    ay = -2 * n * v_chaser_RTN[0]
    az = -n**2 * r_nozzle_RTN[2]
    
    dr_dt = v_chaser_RTN
    dv_dt = np.array([ax, ay, az])
    
    return np.concatenate((dq_dt, dw_dt, dm_dt, dr_dt, dv_dt))

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

def extract_aiming_history(sol, I_target, ctrl):
    """
    Extracts the 3D aiming coordinates and thrust vectors in the RTN frame 
    over the entire simulation using the dynamic chaser position.
    """
    num_steps = len(sol.t)
    
    # Pre-allocate arrays for speed
    aiming_coords = np.zeros((3, num_steps))
    thrust_vectors = np.zeros((3, num_steps))
    
    for i in range(num_steps):
        q_B2RTN = sol.y[0:4, i]
        w_B = sol.y[4:7, i]
        r_nozzle_RTN = sol.y[8:11, i] 
        
        h_B = np.dot(I_target, w_B)
        h_norm = np.linalg.norm(h_B)
        
        # If the debris is practically stopped, we aren't aiming
        if h_norm < 0.01:
            continue 
            
        R_B2RTN = quat_to_dcm(q_B2RTN)
        h_L = np.dot(R_B2RTN, h_B)
        
        T_g_L = -ctrl['K_imp'] * (h_L / h_norm)
        
        r_norm = np.linalg.norm(r_nozzle_RTN)
        r_hat = r_nozzle_RTN / r_norm
        
        T_g_L_P = T_g_L - np.dot(T_g_L, r_hat) * r_hat
        T_g_L_P_norm = np.linalg.norm(T_g_L_P)
        
        if T_g_L_P_norm < 1e-6:
            P_h_L = np.zeros(3)
        else:
            P_h_L = np.cross((T_g_L_P / T_g_L_P_norm), r_hat)        
            P_h_L = P_h_L / np.linalg.norm(P_h_L)
            
        aim_point = ctrl['D_imp'] * P_h_L
        pointing_vector = r_nozzle_RTN - aim_point
        
        if np.linalg.norm(pointing_vector) > 0:
            thrust_dir_RTN = -pointing_vector / np.linalg.norm(pointing_vector)
        else:
            thrust_dir_RTN = np.zeros(3)
            
        aiming_coords[:, i] = aim_point
        thrust_vectors[:, i] = thrust_dir_RTN
        
    return aiming_coords, thrust_vectors

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

def quat_to_dcm(q):
    """Fast, raw numpy quaternion to DCM conversion. q = [x, y, z, w]"""
    qx, qy, qz, qw = q
    x2, y2, z2 = qx**2, qy**2, qz**2
    w2 = qw**2
    xy, xz, xw = qx*qy, qx*qz, qx*qw
    yz, yw = qy*qz, qy*qw
    zw = qz*qw
    
    return np.array([
        [w2 + x2 - y2 - z2, 2.0 * (xy - zw), 2.0 * (xz + yw)],
        [2.0 * (xy + zw), w2 - x2 + y2 - z2, 2.0 * (yz - xw)],
        [2.0 * (xz - yw), 2.0 * (yz + xw), w2 - x2 - y2 + z2]
    ])

def calculate_mean_motion(altitude_km):
    """
    Calculates the orbital mean motion (omega) for a circular Earth orbit.
    """
    # Earth physical constants
    mu_earth = 3.986004418e14  # Earth's gravitational parameter [m^3/s^2]
    r_earth = 6378137.0        # Earth's equatorial radius [m]
    
    # Calculate the semi-major axis (orbital radius) in meters
    a = r_earth + (altitude_km * 1000.0)
    
    # Compute mean motion using Kepler's Third Law
    omega = np.sqrt(mu_earth / (a**3))
    
    return omega

def cart2rel(n, t):
    Gamma = np.array([[1, 0, -np.cos(n*t), -np.sin(n*t), 0, 0],
                      [0, 1, 2*np.sin(n*t), -2*np.cos(n*t), 0, 0],
                      [0, 0, 0, 0, np.sin(n*t), np.cos(n*t)],
                      [0, 0, n*np.sin(n*t), -n*np.cos(n*t), 0, 0],
                      [-3*n/2, 0, 2*n*np.cos(n*t), 2*n*np.sin(n*t), 0, 0],
                      [0, 0, 0, 0, n*np.cos(n*t), n*np.sin(n*t)]]) # negative sign at the end
    return Gamma
############################
#
# Main 
#
############################

## Initialize parameters and geometry
hypergolic_config = init_hypergolic_plume()

# Convert dictionary to a static tuple for the Numba compiler
hypergolic_tuple = (
    hypergolic_config['rho_star'], hypergolic_config['r_star'],
    hypergolic_config['phi_0'],    hypergolic_config['theta_0'],
    hypergolic_config['theta_1'],  hypergolic_config['gamma'],
    hypergolic_config['beta'],     hypergolic_config['V_lim'],
    hypergolic_config['V_w'],      hypergolic_config['c_n'],
    hypergolic_config['c_t'],      hypergolic_config['m_dot']  
)

# Target inertia matrix 
I_target = np.diag([106400.0, 106400.0, 4500.0])
I_inv = np.linalg.inv(I_target)

# Create 25m Cylinder Geometry
centroids, normals, areas = create_cylinder_geometry()

target_altitude_km = 800  # Adjust this to simulate different debris orbits
omega = calculate_mean_motion(target_altitude_km)

# Target Initial Conditions 
#initial_quat = np.array([0.0, 0.0, 0.0, 1.0]) 
initial_quat = np.array([0.57735026919, 0.57735026919, 0.57735026919, 0]) 
initial_w = np.radians([3, 3, 1])       
fuel_initial = [0.0]

# Chaser CW Initial Conditions 
h_B = np.dot(I_target, initial_w)

# Generate the rotation matrix from the initial quaternion
R_B2RTN = quat_to_dcm(initial_quat)

# Rotate the angular momentum vector into the RTN (orbit) frame
h_RTN = np.dot(R_B2RTN, h_B)

h_R = h_RTN[0] # Radial component
h_T = h_RTN[1] # Transverse (Along-track) component

# Find the phase angle of the target's tumble in the RTN frame
gamma_0 = np.arctan2(h_T, h_R)

# Design the Passively Safe Orbit 
phi_g = gamma_0 + (np.pi / 2.0)

ade = 9.0 # In-plane scale [m]
adi = 9.0 # Out-of-plane scale [m]

# Generate the Initial CW Cartesian State 
x_0 = ade * np.sin(phi_g)
y_0 = 2.0 * ade * np.cos(phi_g)
z_0 = adi * np.sin(phi_g)

vx_0 = omega * ade * np.cos(phi_g)
vy_0 = -2.0 * omega * ade * np.sin(phi_g)
vz_0 = omega * adi * np.cos(phi_g)

cw_initial = np.array([x_0, y_0, z_0, vx_0, vy_0, vz_0])
da = np.linalg.inv(cart2rel(omega, 0)) @ cw_initial / target_altitude_km

# --- 6. Build Unified 14-Element State Vector ---
initial_state = np.concatenate((initial_quat, initial_w, fuel_initial, cw_initial))

print(f"Target initial tumble phase (gamma_0): {np.degrees(gamma_0):.2f} deg")
print(f"Chaser optimal starting phase (phi_g): {np.degrees(phi_g):.2f} deg")
print(f"Chaser starting position [x,y,z]: [{x_0:.2f}, {y_0:.2f}, {z_0:.2f}] m")

control_params = {
    'D_imp': 10.0, # Aim 10m down the cylinder to get a massive lever arm
    'eps_theta': np.radians(20.0), # Allow firing even if curved surfaces scatter the torque direction up to 90 deg
    'eps_m': 1e-6, # Lower minimum threshold to catch smaller torques
    'K_imp': 1e-3 # Gain of guidance torque 
}

# Setup Time Span for 7200 seconds (2 hours) to allow full momentum decay
t_span = (0.0, 3600.0)
t_eval = np.linspace(t_span[0], t_span[1], 5000)

print("Starting 1-hour numerical integration with ACTIVE Impingement Guidance...")

# Propagate using solve_ivp
sol = solve_ivp(
    optimal_dynamic_diff_equation, 
    t_span, 
    initial_state, 
    t_eval=t_eval, 
    args=(I_target, I_inv, centroids, normals, areas, hypergolic_tuple, control_params, omega),
    rtol=1e-5, 
    atol=1e-7
)

# Remember to normalize the quaternions post-simulation
quats = sol.y[0:4, :]
sol.y[0:4, :] = quats / np.linalg.norm(quats, axis=0)
# --- 5. Extraction ---
final_w = sol.y[4:7, -1]
total_fuel_used = sol.y[7, -1]
chaser_trajectory_RTN = sol.y[8:11, :]

print(f"Integration complete. Success: {sol.success}")
print(f"Final Angular Velocity (deg/s): {np.degrees(final_w)}")
print(f"Total Impingement Fuel Used: {total_fuel_used:.4f} kg")

# --- Plotting ---
time = sol.t
w_history_deg = np.degrees(sol.y[4:7, :])

fuel_history = sol.y[7, :]
m_dot = hypergolic_tuple[11] # Assuming you used the tuple approach from earlier

# Reconstruct the firing history (Duty Cycle 0.0 to 1.0) using the derivative of fuel consumed
firing_history = np.gradient(fuel_history, time) / m_dot

# Create a figure with 2 subplots sharing the X-axis
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

# --- Top Subplot: Angular Velocity ---
ax1.plot(time, w_history_deg[0, :], label=r'$\omega_x$ (Pitch)', linewidth=2, color='#1f77b4')
ax1.plot(time, w_history_deg[1, :], label=r'$\omega_y$ (Yaw)', linewidth=2, color='#ff7f0e')
ax1.plot(time, w_history_deg[2, :], label=r'$\omega_z$ (Roll)', linewidth=2, color='#2ca02c')

ax1.set_title('Target Debris Active Detumble (Closed-Loop Guidance)', fontsize=14, pad=15)
ax1.set_ylabel('Angular Velocity [deg/s]', fontsize=12)
ax1.axhline(0, color='black', linewidth=1, linestyle='--')
ax1.grid(True, linestyle=':', alpha=0.7)
ax1.legend(loc='upper right', fontsize=11)

# --- Bottom Subplot: Firing & Fuel History ---
# Plot the Firing Command (Duty Cycle) on the primary Y-axis
color1 = 'purple'
ax2.fill_between(time, 0, firing_history, color=color1, alpha=0.3)
ax2.plot(time, firing_history, color=color1, linewidth=1.5, label='Thruster Firing (Duty Cycle)')
ax2.set_ylabel('Duty Cycle [0 to 1]', color=color1, fontsize=12)
ax2.tick_params(axis='y', labelcolor=color1)

# LOCK LEFT AXIS TO 0
ax2.set_ylim(0.0, 1.1) 
ax2.grid(True, linestyle=':', alpha=0.7)

# Create a twin Y-axis to overlay the total fuel consumed
ax3 = ax2.twinx()
color2 = 'red'
ax3.plot(time, fuel_history, color=color2, linewidth=2, linestyle='-.', label='Cumulative Fuel Used')
ax3.set_ylabel('Fuel Consumed [kg]', color=color2, fontsize=12)
ax3.tick_params(axis='y', labelcolor=color2)

# LOCK RIGHT AXIS TO 0 (with a 10% top margin so the line doesn't hit the ceiling)
max_fuel = np.max(fuel_history)
ax3.set_ylim(0.0, max_fuel * 1.1)

# Combine legends for the bottom subplot
lines_1, labels_1 = ax2.get_legend_handles_labels()
lines_2, labels_2 = ax3.get_legend_handles_labels()
ax2.legend(lines_1 + lines_2, labels_1 + labels_2, loc='center right', fontsize=11)

ax2.set_xlabel('Time [seconds]', fontsize=12)

plt.tight_layout()
#plt.show()


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
#plt.show()


# --- Extract Aiming History ---
aiming_coords, thrust_vectors = extract_aiming_history(sol, I_target, control_params)
# --- 3D Visualization of the Servicer's Aiming Strategy ---
chaser_x = sol.y[8, :]
chaser_y = sol.y[9, :]
chaser_z = sol.y[10, :]

# Find indices where duty cycle is greater than a small noise threshold (e.g., 0.01)
active_idx = np.where(firing_history > 0.01)[0]

# Downsample the active indices so the 3D arrows don't overlap into a solid block
step = 50  # Plot 1 arrow for every 50 firing steps (Adjust this to make it look best!)
plot_idx = active_idx[::step]

# Extract the specific locations and vectors at those downsampled firing times
X_fire = chaser_x[plot_idx]
Y_fire = chaser_y[plot_idx]
Z_fire = chaser_z[plot_idx]

# Assuming thrust_vectors is a (3, N) array from your extract_aiming_history function
U_plot = thrust_vectors[0, plot_idx]
V_plot = thrust_vectors[1, plot_idx]
W_plot = thrust_vectors[2, plot_idx]


# --- 2. 3D Visualization ---
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Plot the Target Center of Mass (Origin)
ax.scatter(0, 0, 0, color='black', s=150, label='Target CoM (0,0,0)', marker='X')

# Plot the continuous Chaser Orbit Path
ax.plot(chaser_x, chaser_y, chaser_z, color='gray', linestyle='--', alpha=0.5, label='Chaser Orbit Path (CW NMT)')

# Plot the Firing Locations (Dots)
ax.scatter(X_fire, Y_fire, Z_fire, color='red', s=20, label='Active Firing Locations')

# Plot the Thrust Vectors (3D Arrows)
# 'length' scales the arrows, 'normalize=True' ensures they all look uniform regardless of magnitude
ax.quiver(X_fire, Y_fire, Z_fire, U_plot, V_plot, W_plot, 
          color='red', length=2.5, normalize=True, alpha=0.8, 
          arrow_length_ratio=0.3, label='Plume Thrust Direction')

# Formatting
ax.set_title('Chaser Relative Orbit & Plume Impingement Trajectory', fontsize=14, pad=15)
ax.set_xlabel('R (Radial) [m]', labelpad=10)
ax.set_ylabel('T (Along-Track) [m]', labelpad=10)
ax.set_zlabel('N (Cross-Track) [m]', labelpad=10)

# Ensure axes are equally scaled so the elliptical orbit isn't distorted
max_range = np.array([chaser_x.max()-chaser_x.min(), 
                      chaser_y.max()-chaser_y.min(), 
                      chaser_z.max()-chaser_z.min()]).max() / 2.0
mid_x = (chaser_x.max()+chaser_x.min()) * 0.5
mid_y = (chaser_y.max()+chaser_y.min()) * 0.5
mid_z = (chaser_z.max()+chaser_z.min()) * 0.5

ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

# Adjust legend to be outside the plot so it doesn't cover the orbit
ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), fontsize=11)
plt.tight_layout()
plt.show()
thrust_output_array = np.column_stack((U_plot, V_plot, W_plot))

U_thrust = thrust_vectors[0, :]
V_thrust = thrust_vectors[1, :]
W_thrust = thrust_vectors[2, :]
#thrust_output_array.tofile('thrust_output.csv') # Idk why this isn't working
thrust_output = pd.DataFrame(thrust_output_array)
thrust_path = r'C:\Users\Athar\Documents\College\AAE450-Spacecraft-Design-Team-OD5\ADCS\Detumble\Thrust_Output.csv'
thrust_output.to_csv(thrust_path, 'a') # Idk why this isn't working
print("Finished Sim")

# Generating output for Animation CSV
time_array = sol.t
q_x, q_y, q_z, q_w = sol.y[0:4, :]

chaser_x, chaser_y, chaser_z = sol.y[8:11, :]

# Convert into a clean 1 (Firing) or 0 (Off) integer array based on your 1% deadband
is_thrusting = (firing_history > 0.01).astype(int)

# --- 5. Assemble the Master Array ---
master_output_array = np.column_stack((
    time_array, 
    q_x, q_y, q_z, q_w, 
    chaser_x, chaser_y, chaser_z, 
    U_thrust, V_thrust, W_thrust, 
    is_thrusting
))

# --- 6. Save to CSV ---
master_path = r'C:\Users\Athar\Documents\College\AAE450-Spacecraft-Design-Team-OD5\ADCS\Detumble\Detumble_History.csv'

# Create a clean, comma-separated header
header_str = "Time,q_x,q_y,q_z,q_w,chaser_x,chaser_y,chaser_z,U_thrust,V_thrust,W_thrust,is_thrusting"

np.savetxt(
    master_path, 
    master_output_array, 
    delimiter=",", 
    header=header_str, 
    comments="", 
    fmt="%.6f" # Formats all floats to 6 decimal places to keep the file size tight
)