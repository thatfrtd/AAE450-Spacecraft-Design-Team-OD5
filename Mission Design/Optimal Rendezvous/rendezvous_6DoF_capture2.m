%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Nonlinear Rendezvous Example
% Author: Travis Hastreiter 
% Created On: 28 February, 2026
% Description: Sequential Convex Programming Trajectory Optimization for 
% general rendezvous using relative orbital motion equations. Includes mass 
% in the state. You must have CVX installed. Multiple thrusters with
% different thrust levels and Isp but using the same fuel
% (resistojet/arcjet eprop with chemical RCS both using hydrazine).
% Last Modified On: 28 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

char_star = load_charecteristic_values_Earth();
nd_scalar = [ones([13, 1]); char_star.m];
g_0 = 9.81; % [km / s2]
R_E = char_star.l;

% Estimate spacecraft moment of inertia
cylinder_MoI = @(m, r, h) [1 / 2 * m * r ^ 2; ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2); ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2)];

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = [4155; 303]; % [s]
spacecraft_params.m_0 = 1500; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.r = 1; % [m]
spacecraft_params.h = 2; % [m]
spacecraft_params.I = cylinder_MoI(spacecraft_params.m_0, spacecraft_params.r, spacecraft_params.h);
spacecraft_params.F_max = [0.235, 22]; % [N]
% RCS (direction axis, position axis) x+y, x-y, x+z, x-z, y+z, y-z, y+x, y-x, z+x, z-x, z+y, z-y
spacecraft_params.thrust_directions = {[1; 0; 0], ... % Main thruster
                                       [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
                                        0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0;
                                        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]}; % RCS thruster directions
spacecraft_params.thrust_positions = {[1; 0; 0], ... % Main thruster
                                      [0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0;
                                       1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1;
                                       0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0] * spacecraft_params.r}; % RCS thruster positions
spacecraft_params.R_C_S = calculate_RCS_allocation_matrix(spacecraft_params.thrust_directions{2}, spacecraft_params.thrust_positions{2}); % RCS thruster force and torque allocation matrix
spacecraft_params.tau_max = 0.2; % [N m] Max reaction wheel torque
F_max_nd = spacecraft_params.F_max; % F_max in N, char_star.F in kN
tau_max_nd = spacecraft_params.tau_max; 
tau_max_RCS = 30; % [N m] Max allowable RCS torque (could make soft constraint lambda_RCS_tau * max(RCS_torque - tau_max_RCS, 0) )
max_angvel = deg2rad(5);

% Spacedebris Parameters
spacedebris_params = struct();
spacedebris_params.Isp = [4100; 232]; % [s]
spacedebris_params.m = 4000; % [kg]
spacedebris_params.r = 1.85; % [m]
spacedebris_params.h = 11.9; % [m]
spacedebris_params.I = cylinder_MoI(spacedebris_params.m, spacedebris_params.r, spacedebris_params.h);

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = R_E + 800; % [km] semi-major axis
e_c = 0.003; % [] eccentricity
i_c = deg2rad(10); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu0_c = deg2rad(0); % [rad] true anomaly at epoch
M0_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu0_c, e_c), e_c);
x_keplerian_c = [a_c; e_c; i_c; Omega_c; omega_c; M0_c];
n_c = sqrt(char_star.mu / a_c ^ 3);

% Rendezvous time
tf = 120; % [s] (nondimensionalized)
N = 60;
t_k_actual = linspace(0, tf, N);
tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

% Initial conditions for spacecraft - specify orbit instead?
data_table = readtable("Full_Detumble_History.csv");
rv_0 = [data_table.chaser_x(end); data_table.chaser_y(end); data_table.chaser_z(end); data_table.chaser_vx(end); data_table.chaser_vy(end); data_table.chaser_vz(end)] * 1e-3./ nd_scalar(1:6);
theta_0 = [deg2rad(0); deg2rad(45); deg2rad(45)]; % [rad]
R_0 = angle2dcm(theta_0(1), theta_0(2), theta_0(3));
r_0 = rv_0(1:3);
v_0 = rv_0(4:6);
%q_0 = qExp(RLog(R_0));
q_0 = point_at_vec(-r_0);
w_0 = deg2rad([0; 0; 0]); % [rad / s]
x_0 = [rv_0; q_0; w_0; spacecraft_params.m_0] ./ nd_scalar;

% Terminal conditions
seperation_distance = 5; % [m]
d_vec = [-(spacedebris_params.h / 2 + spacecraft_params.h / 2 + seperation_distance); 0; 0]*1e-3;
%w0_debris = deg2rad(-[0.0739732221634118; 0.0063002539406854; 0.0109314217396326]);
w0_debris = deg2rad([0.0739732221634118; 0.3002539406854; 0.2])*2;
q0_debris = [data_table.q_x(end); data_table.q_y(end); data_table.q_z(end); data_table.q_w(end)];
[x_f_interstellar, qw_debris_array] = interstellar_boundary_conditions(q0_debris, w0_debris, t_k, spacedebris_params.I, d_vec);
x_f = x_f_interstellar ./ nd_scalar(1:13);
q_f = x_f_interstellar(7:10);
r_f = x_f_interstellar(1:3);

%% Initialize

% Control discretization method
% Zero Order Hold (ZOH) - Piecewise constant
% First Order Hold (FOH) - Piecewise linear
u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

parser = "CVX"; % Can use CVXPY for more speed if needed (needs setup)
nx = 14; % Number of states (r, v, q, w, m)
nu = 16; % Number of controls (T_1_mag, T_2, tau_rw)
np = 0; % Number of parameters (tf, v_0, etc)

% PTR algorithm parameters
ptr_ops.iter_max = 20;
ptr_ops.iter_min = 1;
ptr_ops.Delta_min = 2e-2;
ptr_ops.w_vc = 5e5;
ptr_ops.w_tr = ones(1, N) * 5e-2;
ptr_ops.w_tr_p = 0;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 0e-2;
ptr_ops.alpha_p = 0;

% Scaling currently not helping...
scale = false;
% 
% scale_hint.x_max = [max(x_0(1:3)) * ones([3, 1]); max(x_0(4:6)) * ones([3, 1]); spacecraft_params.m_0 / char_star.m];
% scale_hint.x_min = [-max(x_0(1:3)) * ones([3, 1]); -max(x_0(4:6)) * ones([3, 1]); spacecraft_params.m_0 / char_star.m * 0.95];
% scale_hint.u_max = [F_max_nd * ones([6, 1])];
% scale_hint.u_min = [zeros([3, 1])];
% scale_hint.p_max = [];
% scale_hint.p_min = [];

%% Get Dynamics
f = @(t, x, u, p) relative_orbit_6DoF_twothruster_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp; char_star.mu; g_0; spacecraft_params.I; spacecraft_params.thrust_directions{1}; spacecraft_params.R_C_S(:)]);

%% Specify Constraints
mass_constraint = {1:N, @(t, x, u, p) spacecraft_params.m_dry / char_star.m - x(14)};
angular_velocity_constraint = {1:N, @(t, x, u, p) norm(x(11:13)) - max_angvel};
final_position_constraint = {N, @(t, x, u, p) norm(x(1:3)) - norm(r_f)};
state_convex_constraints = {mass_constraint, angular_velocity_constraint, final_position_constraint};

% Convex control constraints
max_thrust_constraint_1 = {1:Nu, @(t, x, u, p) u(1) - F_max_nd(1)};
min_thrust_constraint_1 = {1:Nu, @(t, x, u, p) 0 - u(1)};
max_RCS_thrust_constraint = {1:Nu, @(t, x, u, p) norm(u(2:13), Inf) - F_max_nd(2)};
max_RCS_torque_constraint = {1:Nu, @(t, x, u, p) norm(spacecraft_params.R_C_S * u(2:13), Inf) - tau_max_RCS};
max_reaction_wheel_torque_constraint = {1:Nu, @(t, x, u, p) norm(u(14:16), Inf) - tau_max_nd};
control_convex_constraints = {max_thrust_constraint_1, min_thrust_constraint_1, max_RCS_thrust_constraint, max_reaction_wheel_torque_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
% Keep out sphere
keep_out_distance = 0.010; % [km]
keep_out_sphere_constraint = @(t, x, u, p) keep_out_distance ^ 2 - (x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2);
keep_out_sphere_constraint_linearized = {1:N, linearize_constraint(keep_out_sphere_constraint, nx, nu, np, "x", 1:3)};
% Keepout Ellipsoid
keep_out_ellipsoid_form = ([spacedebris_params.h / 2; spacedebris_params.r; spacedebris_params.r] + seperation_distance)*1e-3 ./ nd_scalar(1:3); % [km] 
keep_out_ellipsoid_constraint = @(t, x, u, p) 1 - ((x(1) / keep_out_ellipsoid_form(1)) ^ 2 + (x(2) / keep_out_ellipsoid_form(2)) ^ 2 + (x(3) / keep_out_ellipsoid_form(3)) ^ 2);
keep_out_ellipsoid_constraint_linearized_func = linearize_constraint(keep_out_ellipsoid_constraint, nx, nu, np, "x", 1:3);
keep_out_ellipsoid_constraint_linearized_STC = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) keep_out_ellipsoid_constraint_linearized_func(t, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x(1:3), u, p, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x_ref(1:3, :), u_ref, p_ref, k)*1e5};
% Line of sight constraint
d_B_camera = [1; 0; 0]*1e-3; % Location of camera in body frame
p_B_camera = [1; 0; 0]; % Sensor boresight direction in body frame
angle_LOS_camera = deg2rad(10); % [rad]
camera_LOS_trigger_distance = 0.08; % [km]
camera_LOS_constraint = @(t, x, u, p) (quat_rot(q_conj(x(7:10)), x(1:3)) + d_B_camera).' * p_B_camera + dnorm(quat_rot(q_conj(x(7:10)), x(1:3)) + d_B_camera) * cos(angle_LOS_camera);
camera_LOS_constraint_linearized_func = linearize_constraint(camera_LOS_constraint, nx, nu, np, "x", 1:nx);
camera_LOS_constraint_linearized = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) camera_LOS_constraint_linearized_func(t, x, u, p, x_ref, u_ref, p_ref, k) * max(camera_LOS_trigger_distance - norm(x_ref(1:3, k)), 0)*1e12};
camera_LOS_constraint_relaxed = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) (quat_rot(q_conj(x_ref(7:10, k)), x(1:3)) + d_B_camera).' * p_B_camera + norm(x(1:3) + quat_rot(x_ref(7:10, k), d_B_camera)) * cos(angle_LOS_camera)};
% Combine nonconvex state constraints
state_nonconvex_constraints = {keep_out_ellipsoid_constraint_linearized_STC, camera_LOS_constraint_linearized};

% Nonconvex control constraints
% Plume impingement constraints
control_nonconvex_constraints = {};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Boundary conditions
initial_bc = @(x, p) [x - x_0];
terminal_bc = @(x, p, x_ref, p_ref) [x(1:13) - x_f(1:13); 0]; % Don't constrain final mass
% terminal_bc = @(x, p, x_ref, p_ref) [4 * x(1) * nd_scalar(1) + 2 / n_c * x(5) * nd_scalar(5); 
%                                          x(2) * nd_scalar(2) - 2 / n_c * x(4) * nd_scalar(4); 
%                                          x(7:13) - x_f(7:13); x(3); x(6); 0; 0; 0]; % Don't constrain final mass

%% Specify Objective
objective_min_fuel = @(x, u, p, x_ref, u_ref, p_ref) sum(norms(u(2:13, :), 1)) * delta_t / (spacecraft_params.Isp(2) * g_0) ...
                                                   + sum(norms(u(14:16, :), 1)) * delta_t * 1e-6 ...
                                                   + sum(norms(u(1, :), 1)) * delta_t / (spacecraft_params.Isp(1) * g_0)*1000;

%% Create Guess
% Straight Line Initial Guess
rv_int = linspace(0, 1, N) .* (x_f(1:6) - x_0(1:6)) + x_0(1:6);
[q_int, w_int_const] = q_interp(q_0, q_f, t_k);
w_int = repmat(w_int_const, 1, numel(t_k));
guess.x = [rv_int; q_int; w_int; x_0(14) * ones([1, N])];%linspace(0, 1, N) .* ([x_f; x_0(14)] - x_0) + x_0;
guess.u = ones([nu, Nu]) * 1e-6;
guess.p = [];

%% Construct Problem Object
problem = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, objective_min_fuel, scale = scale, nonconvex_constraints = nonconvex_constraints, initial_bc = initial_bc, terminal_bc = terminal_bc, integration_tolerance = 1e-12, discretization_method = "error", N_sub = 1, Name = "rendezvous_6DoF_twothruster");

[problem, Delta_disc] = problem.discretize(guess.x, guess.u, guess.p);

%% Test with known solution
% t_cont_sol = linspace(0, tf, 1000);
% [~, x_cont_sol_ck] = problem.cont_prop([zeros([1, Nu]); ptr_sol.u(:, :, i)], ptr_sol.p(:, i), tspan = t_cont_sol);
% 
% figure
% tiledlayout(3, 1)
% 
% nexttile
% plot(t_cont_sol, x_cont_sol_ck(7:10, :))
% xlabel("Time")
% ylabel("Quaternions")
% title("Quaternion History")
% grid on
% 
% nexttile
% plot(t_cont_sol, rad2deg(x_cont_sol_ck(11:13, :)))
% xlabel("Time")
% ylabel("Angular Velocities")
% title("Angular Velocity History")
% grid on
% 
% % Plot Control
% nexttile
% plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(13:15,:), LineWidth=1); hold on
% plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(13:15,:)), LineWidth=1)
% title("Reaction Wheel Torque")
% xlabel("Time")
% ylabel("Torque [N m]")
% legend("\hat{b}_1", "\hat{b}_2", "\hat{b}_3", "||\tau||", Interpreter="latex")
% grid on

%% Solve
ptr_sol = ptr(problem, ptr_ops, parser);

%% Plot Convergence
if ptr_sol.converged == false
    ptr_sol.converged_i = ptr_ops.iter_max;
end

figure
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [problem.objective(problem.guess.x, problem.guess.u, problem.guess.p), [ptr_sol.info.J]]); hold on
hold off
legend("PTR Iterations")
title("Objective vs Iteration")
grid on

nexttile
plot(ptr_sol.delta_xp)
title("Stopping Criteria vs Iteration")
yscale("log")
grid on

nexttile
plot(0:ptr_sol.converged_i, squeeze(sum(vecnorm(ptr_sol.Delta(:, :, 1:(ptr_sol.converged_i + 1)), 2, 1))))
yscale("log")
title("Defect Norm vs Iteration")
grid on

%% Extract Solution
i = ptr_sol.converged_i + 1;
x = ptr_sol.x(:, :, i) .* nd_scalar;
u = ptr_sol.u(:, :, i);

t_cont_sol = linspace(0, tf, 1000);
[~, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i), tspan = t_cont_sol);
x_cont_sol = x_cont_sol .* nd_scalar;

%% Plot Trajectory

u_main_vec = quat_rot_array(x(7:10, 1:Nu), spacecraft_params.thrust_directions{1} .* u(1, :)); % Rotate 
u_cont_sol_main_vec = quat_rot_array(x_cont_sol(7:10, 1:end - (N - Nu)), spacecraft_params.thrust_directions{1} .* u_cont_sol(1, :)); % Rotate 
u_RCS_vec = squeeze(quat_rot_array(x(7:10, 1:Nu), pagemtimes(spacecraft_params.R_C_S, reshape(u(2:13, :), 12, 1, [])))); % Rotate 
u_cont_sol_RCS_vec = squeeze(quat_rot_array(x_cont_sol(7:10, 1:end - (N - Nu)), pagemtimes(spacecraft_params.R_C_S, reshape(u_cont_sol(2:13, :), 12, 1, [])))); % Rotate 
for k = 1 : Nu
    u_norm = norm(u_main_vec(1:3, :));
    if u_norm < 1e-8
        u_main_vec(1:3, k) = 0;
    end
end

figure
scatter3(0, 0, 0, 60, "blue", "filled", "diamond"); hold on
plot_cartesian_orbit(x_cont_sol(1:3,:)', 'k', 0.4, 1); hold on
quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u_main_vec(1, :), u_main_vec(2, :), u_main_vec(3, :), 1, "filled", Color = "red")
quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u_RCS_vec(1, :), u_RCS_vec(2, :), u_RCS_vec(3, :), 2, "filled", Color = "m")
scatter3(r_0(1), r_0(2), r_0(3), 48, "green", "filled", "square"); hold on
scatter3(r_f(1), r_f(2), r_f(3), 48, "red", "x"); hold on
plot3(guess.x(1, :), guess.x(2, :), guess.x(3, :), Color = "green", LineStyle = "--");
for k = 1 :10: N
    [x_1, y_1, z_1] = ellipsoid(0,0,0,keep_out_ellipsoid_form(1),keep_out_ellipsoid_form(2),keep_out_ellipsoid_form(3));
    h = surf(x_1, y_1, z_1);
    tau = qLog(qw_debris_array(1:4, k));
    theta = norm(tau);
    if theta > 1e-10
        lambda = tau / theta;
    else
        lambda = [1; 0; 0];
    end
    rotate(h, lambda, -rad2deg(theta), [0; 0; 0]);
    set(h, 'FaceAlpha', 0.1)
end
keep_out_ellipsoid_form = ([spacedebris_params.h / 2; spacedebris_params.r; spacedebris_params.r])*1e-3 ./ nd_scalar(1:3); % [km] 
[x_1, y_1, z_1] = ellipsoid(0,0,0,keep_out_ellipsoid_form(1),keep_out_ellipsoid_form(2),keep_out_ellipsoid_form(3));
h2 = surf(x_1, y_1, z_1);
tau = qLog(q_f);
theta = norm(tau);
if theta > 1e-10
    lambda = tau / theta;
else
    lambda = [1; 0; 0];
end
rotate(h2, lambda, -rad2deg(theta), [0; 0; 0]);
set(h2, 'FaceAlpha', 0.5)
shading interp
title('Optimal Rendezvous')
xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
axis equal
legend("Target", 'Spacecraft', "", "Thrust 1", "Thrust 2", "Start", "End", "Guess", "Keepout Sphere", 'Location', 'northwest'); axis equal; grid on

%% Plot Control
figure
tiledlayout(3, 1)

nexttile
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol_main_vec, LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(1, :), LineWidth=1)
title("Main Thruster")
xlabel("Time")
ylabel("Force [N]")
legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
grid on

nexttile
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol_RCS_vec, LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol_RCS_vec), LineWidth=1)
title("Reaction Control Thrusters")
xlabel("Time")
ylabel("Force [N]")
legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
grid on

nexttile
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(14:16,:), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(14:16,:)), LineWidth=1)
title("Reaction Wheel Torque")
xlabel("Time")
ylabel("Torque [N m]")
legend("\hat{b}_1", "\hat{b}_2", "\hat{b}_3", "||\tau||", Interpreter="latex")
grid on

%% Plot RCS Thrusters
% RCS (direction axis, position axis) x+y, x-y, x+z, x-z, y+z, y-z, y+x, y-x, z+x, z-x, z+y, z-y
figure
tiledlayout(2, 3)

nexttile % xy
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(2, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(3, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("x+y", "x-y")
title("X Directed Y Positioned")
grid on
ylim(spacecraft_params.F_max(2) * [-1; 1])

nexttile % yz
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(6, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(7, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("y+z", "y-z")
title("Y Directed Z Positioned")
grid on
ylim(spacecraft_params.F_max(2) * [-1; 1])

nexttile % zx
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(10, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(11, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("z+x", "z-x")
title("Z Directed X Positioned")
grid on
ylim(spacecraft_params.F_max(2) * [-1; 1])

nexttile % xz
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(4, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(5, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("x+z", "x-z")
title("X Directed Z Positioned")
grid on
ylim(spacecraft_params.F_max(2) * [-1; 1])

nexttile % yx
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(8, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(9, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("y+x", "y-x")
title("Y Directed X Positioned")
grid on
ylim(spacecraft_params.F_max(2) * [-1; 1])

nexttile % zy
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(12, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(13, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("z+y", "z-y")
title("Z Directed Y Positioned")
grid on
ylim(spacecraft_params.F_max(2) * [-1; 1])

sgtitle("Reaction Control System Thrusters Control History")
%% Plot Camera Line of Sight
camera_LoS = acosd(dot(quat_rot_array(x_cont_sol(7:10, :), repmat([1; 0; 0], 1, numel(t_cont_sol))), -x_cont_sol(1:3, :) ./ vecnorm(x_cont_sol(1:3, :))));

figure
plot(t_cont_sol, camera_LoS)
yline(rad2deg(angle_LOS_camera))
xline(t_k(camera_LOS_constraint_linearized{1}(1)))
xlabel("Time [s]")
ylabel("Line of Sight Angle [deg]")
grid on
title("Camera Line of Sight")
legend("Solution", "Max", "Activation")
%%
for k = 1 : N
    (quat_rot(q_conj(x(7:10, k)), x(1:3, k)) + d_B_camera).' * p_B_camera + dnorm(quat_rot(q_conj(x(7:10, k)), x(1:3, k)) + d_B_camera) * cos(angle_LOS_camera)
end

%% Save to CSV for Blender Animation
[~, qw_debris_array_cont] = interstellar_boundary_conditions(q0_debris, w0_debris, t_cont_sol, spacedebris_params.I, d_vec);
output_array = [t_cont_sol(1:(end-1))', x_cont_sol(:, 1:(end-1))', qw_debris_array_cont(1:4, 1:(end - 1))', u_cont_sol'];
state_names = ["r_Hill_1", "r_Hill_2", "r_Hill_3", ...
               "v_Hill_1", "v_Hill_2", "v_Hill_3", ...
               "q_1", "q_2", "q_3", "q_4", ...
               "w_1", "w_2", "w_3", ...
               "mass", ...
               "q_1_d", "q_2_d", "q_3_d", "q_4_d"];
control_names = strings(1, nu);
control_names(1) = "u_main";
for i = 1 : size(spacecraft_params.thrust_directions{2}, 2)
    control_names(i + 1) = sprintf("u_RCS_%g", i);
end
for i = 1 : 3
    control_names(i + size(spacecraft_params.thrust_directions{2}, 2) + 1) = sprintf("u_W_%g", i);
end
output_names = ["Time", state_names, control_names];

output_table = array2table(output_array, VariableNames = output_names);
writetable(output_table,"./Animation/6DoF_capture.csv")

%% Iteration Plot
fig = scatter3(0, 0, 0, 60, "blue", "filled", "diamond"); hold on
C = cool(ptr_sol.converged_i);
t_iters = {};
x_iters = {};
u_iters = {};
linestyle = string(1:ptr_sol.converged_i);
for i = 1:ptr_sol.converged_i
    t_iters{i} = t_k;
    x_iters{i} = ptr_sol.x(:, :, i);
    u_iters{i} = ptr_sol.u(:, :, i);
    linestyle(i) = "-";
    patch([x_iters{i}(1,:) nan],[x_iters{i}(2,:) nan],[x_iters{i}(3,:) nan],[repmat(i / ptr_sol.converged_i, 1, numel(x_iters{i}(3,:))) nan],'FaceColor','none','EdgeColor','flat')
end
%comparison_plot_6DoF_trajectory(x_iters, "iter " + string(1:ptr_sol.converged_i), glideslope_angle_max, linestyle = linestyle, title = "Solution vs Iteration")

%quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(1, :), u(2, :), u(3, :), 1, "filled", Color = "red")

% quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(4, :), u(5, :), u(6, :), 2, "filled", Color = "m", LineWidth=2)
scatter3(x_0(1) * nd_scalar(1), x_0(2) * nd_scalar(2), x_0(3) * nd_scalar(3), 48, "green", "filled", "square"); hold on
scatter3(x_f(1) * nd_scalar(1), x_f(2) * nd_scalar(2), x_f(3) * nd_scalar(3), 48, "red", "x"); hold on
plot3(guess.x(1, :), guess.x(2, :), guess.x(3, :), Color = "green", LineStyle = "--");
C = jet(ptr_sol.converged_i);
colormap(C);
% for k = 1 : N
%     [x_1, y_1, z_1] = ellipsoid(0,0,0,keep_out_ellipsoid_form(1),keep_out_ellipsoid_form(2),keep_out_ellipsoid_form(3));
%     h = surf(x_1, y_1, z_1);
%     tau = qLog(qw_debris_array(1:4, k));
%     theta = norm(tau);
%     if theta > 1e-10
%         lambda = tau / theta;
%     else
%         lambda = [1; 0; 0];
%     end
%     rotate(h, lambda, -rad2deg(theta), [0; 0; 0]);
%     set(h, 'FaceAlpha', 0.1)
% end
keep_out_ellipsoid_form = ([spacedebris_params.h / 2; spacedebris_params.r; spacedebris_params.r])*1e-3 ./ nd_scalar(1:3); % [km] 
[x_1, y_1, z_1] = ellipsoid(0,0,0,keep_out_ellipsoid_form(1),keep_out_ellipsoid_form(2),keep_out_ellipsoid_form(3));
h2 = surf(x_1, y_1, z_1);
tau = qLog(q_f);
theta = norm(tau);
if theta > 1e-10
    lambda = tau / theta;
else
    lambda = [1; 0; 0];
end
rotate(h2, lambda, -rad2deg(theta), [0; 0; 0]);
set(h2, 'FaceAlpha', 1)
shading interp

title('SCP Iteration Trajectories for Capture Approach')

xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
axis equal

c = colorbar();
c.Label.String = 'Iteration';
c.TickLabels = string(0:ptr_sol.converged_i);


%% Get Net Force and Torque for Ansys Transient Modal Structural
u_cont_sol_RW = u_cont_sol(14:16,:);
RCS_net = pagemtimes(spacecraft_params.R_C_S, reshape(u_cont_sol(2:13, :), 12, 1, []));

F_net = RCS_net(1:3, :);
tau_net = RCS_net(4:6, :) + u_cont_sol_RW;

output_array = [t_cont_sol(1:(end-1))', F_net', tau_net'];
control_names = strings(1, 6);
for i = 1 : 3
    control_names(i) = sprintf("u_F_%g", i);
end
for i = 1 : 3
    control_names(i + 3) = sprintf("u_tau_%g", i);
end
output_names = ["Time", control_names];

output_table = array2table(output_array, VariableNames = output_names);

writetable(output_table,"./Structures/6DoF_capture_u_net.csv")

%% Helper Functions
function [q] = qExp(tau)
    N = size(tau, 2);
    q = zeros([4, N]);
    for k = 1 : N
        theta = norm(tau(:, k));
        if theta > 1e-10
            u = tau(:, k) / theta;
            q(:, k) = [u * sin(theta / 2); cos(theta / 2)];
        else
            q(:, k) = [0; 0; 0; 1];
        end
    end
end

function [tau] = qLog(q)
    N = size(q, 2);
    for k = 1 : N
        w = q(4, k);
        v = q(1:3, k);
        w = w * sign(w);
        v = v * sign(w);
        tau(:, k) = 2 * v * atan2(norm(v), w) / norm(v);
    end
end

function [tau] = RLog(R)
    theta = acos((trace(R) - 1) / 2);
    u = vee(R - R') / (2 * sin(theta));

    tau = theta * u;

    function [tau] = vee(tau_hat)
        tau = [tau_hat(3, 2); tau_hat(1, 3); tau_hat(2, 1)];
    end
end

function [q_interp, w_diff] = q_interp(q_0, q_f, tspan)
    q_diff = q_mul(q_0, q_conj(q_f));

    theta_diff = qLog(q_diff);

    w_diff = theta_diff ./ tspan(end);

    tau_interp = qExp(interp1([0, 1]', [zeros([3, 1]), theta_diff]', linspace(0, 1, numel(tspan)))');

    q_interp = q_mul_array(repmat(q_0, 1, numel(tspan)), tau_interp);
end

function [quat] = point_at_vec(T_vec)
    T_vec = T_vec ./ vecnorm(T_vec);
    w = 1 + T_vec(1, :);
    xyz = cross([ones([1, size(T_vec, 2)]); zeros([2, size(T_vec, 2)])], T_vec); % cross b_x with T_vec
    quat = [xyz; w];
    quat = quat ./ vecnorm(quat);
end