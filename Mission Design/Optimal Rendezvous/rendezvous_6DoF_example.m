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
nd_scalar = [char_star.l * ones([3, 1]); char_star.v * ones([3, 1]); ones([7, 1]); char_star.m];

% Estimate spacecraft moment of inertia
cylinder_MoI = @(m, r, h) [1 / 2 * m * r ^ 2; ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2); ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2)];

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = [300; 232]; % [s]
spacecraft_params.m_0 = 800; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.r = 1.85; % [m]
spacecraft_params.h = 4; % [m]
spacecraft_params.I = cylinder_MoI(spacecraft_params.m_0, spacecraft_params.r, spacecraft_params.h);
spacecraft_params.F_max = [0.8; 22]; % [N]
spacecraft_params.thrust_directions = {[1; 0; 0], ... % Main thruster
                                       {}}; % RCS thrusters
spacecraft_params.tau_max = spacecraft_params.F_max(2) * spacecraft_params.r * 4; % [N m]
F_max_nd = spacecraft_params.F_max / 1000 / char_star.F; % F_max in N, char_star.F in kN
tau_max_nd = spacecraft_params.tau_max;% / 1000 / char_star.F / char_star.l; 

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = 6728; % [km] semi-major axis
e_c = 0.01; % [] eccentricity
i_c = deg2rad(10); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu0_c = deg2rad(0); % [rad] true anomaly at epoch
M0_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu0_c, e_c), e_c);
x_keplerian_c = [a_c; e_c; i_c; Omega_c; omega_c; M0_c];

% Rendezvous time
tf = 15000 / char_star.t; % [s] (nondimensionalized)

% Initial conditions for spacecraft - specify orbit instead?
r_0 = [0.5; -0.5; 0.2]; % [km]
v_0 = [0.001; 1e-3; 0]; % [km / s]
theta_0 = [deg2rad(0); deg2rad(90); deg2rad(0)]; % [rad]
R_0 = angle2dcm(theta_0(1), theta_0(2), theta_0(3));
q_0 = qexp(RLog(R_0));
w_0 = deg2rad([0; 0; 0]); % [rad / s]
x_0 = [r_0; v_0; q_0; w_0; spacecraft_params.m_0] ./ nd_scalar;

% Terminal conditions
r_f = [0; 0.2; 0]; % [km]
v_f = [0e-3; 0; 0]; % [km / s]
theta_f = deg2rad([0; 0; 90]); % [rad]
R_f = angle2dcm(theta_f(1), theta_f(2), theta_f(3));
q_f = qexp(RLog(R_f));
w_f = deg2rad([0; 0; 0]); % [rad / s]
x_f = [r_f; v_f; q_f; w_f] ./ nd_scalar(1:13);

%% Initialize
N = 200;
t_k_actual = linspace(0, tf, N);
tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

% Control discretization method
% Zero Order Hold (ZOH) - Piecewise constant
% First Order Hold (FOH) - Piecewise linear
u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

parser = "CVX"; % Can use CVXPY for more speed if needed (needs setup)
nx = 14; % Number of states (r, v, q, w, m)
nu = 7; % Number of controls (T_1_mag, T_2, tau_rw)
np = 0; % Number of parameters (tf, v_0, etc)

% PTR algorithm parameters
ptr_ops.iter_max = 25;
ptr_ops.iter_min = 1;
ptr_ops.Delta_min = 5e-7;
ptr_ops.w_vc = 5e5;
ptr_ops.w_tr = ones(1, Nu) * 5e-2;
ptr_ops.w_tr_p = 0;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 0;
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
f = @(t, x, u, p) relative_orbit_6DoF_twothruster_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp; spacecraft_params.I; spacecraft_params.thrust_directions{1}]);

%% Specify Constraints
mass_constraint = {1:N, @(t, x, u, p) spacecraft_params.m_dry / char_star.m - x(14)};
%angular_velocity_constraint = {1:N, @(t, x, u, p) norm(x(11:13), Inf) - norm([w_0; deg2rad(20)], Inf)};
% sensor seeing target
state_convex_constraints = {mass_constraint};

% Convex control constraints
max_thrust_constraint_1 = {1:N, @(t, x, u, p) u(1) - F_max_nd(1)};
min_thrust_constraint_1 = {1:N, @(t, x, u, p) 0 - u(1)};
max_thrust_constraint_2 = {1:N, @(t, x, u, p) norm(u(2:4)) - F_max_nd(2)};
max_reaction_wheel_torque_constraint = {1:N, @(t, x, u, p) norm(u(5:7)) - tau_max_nd};
control_convex_constraints = {max_thrust_constraint_1, min_thrust_constraint_1, max_thrust_constraint_2, max_reaction_wheel_torque_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
keep_out_distance = 0.2; % [km]
keep_out_sphere_constraint = @(t, x, u, p) keep_out_distance ^ 2 - nd_scalar(1) ^ 2 * (x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2);
keep_out_sphere_constraint_linearized = {1:N, linearize_constraint(keep_out_sphere_constraint, nx, nu, np, "x", 1:3)};
state_nonconvex_constraints = {keep_out_sphere_constraint_linearized};

% Nonconvex control constraints
% Plume impingement constraints
control_nonconvex_constraints = {};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Boundary conditions
initial_bc = @(x, p) [x - x_0];
terminal_bc = @(x, p, x_ref, p_ref) [x([1:6, 11:13]') - x_f([1:6, 11:13]'); zeros([4, 1]); 0]; % Don't constrain final mass

%% Specify Objective
objective_min_fuel = @(x, u, p, x_ref, u_ref, p_ref) sum(u(1, :)) * delta_t * char_star.F / (spacecraft_params.Isp(1) * 9.81e-3) * 1000 ...
                                                   + sum(norms(u(2:4, :))) * delta_t * char_star.F / (spacecraft_params.Isp(2) * 9.81e-3) * 1000 ...
                                                   + sum(norms(u(5:7, :), 1))*0;

%% Create Guess
% Straight Line Initial Guess - Definitely need to actually slerp for
% proper quaternion interpolation
guess.x = linspace(0, 1, N) .* ([x_f; x_0(14)] - x_0) + x_0;
guess.u = ones([nu, Nu]) * 1e-6;
guess.p = [];

%% Construct Problem Object
problem = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, objective_min_fuel, scale = scale, nonconvex_constraints = nonconvex_constraints, initial_bc = initial_bc, terminal_bc = terminal_bc, integration_tolerance = 1e-12, discretization_method = "error", N_sub = 1, Name = "rendezvous_6DoF_twothruster");

[problem, Delta_disc] = problem.discretize(guess.x, guess.u, guess.p);

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
u = ptr_sol.u(:, :, i) * char_star.F * 1000;

[t_cont_sol, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));
t_cont_sol = t_cont_sol * char_star.t;
x_cont_sol = x_cont_sol .* nd_scalar;
u_cont_sol = u_cont_sol * char_star.F * 1000;

%% Plot Trajectory

u_main_vec = quat_rot_array(x(7:10, 1:Nu), spacecraft_params.thrust_directions{1} .* u(1, :)); % Rotate 
u_cont_sol_main_vec = quat_rot_array(x_cont_sol(7:10, 1:end - (N - Nu)), spacecraft_params.thrust_directions{1} .* u_cont_sol(1, :)); % Rotate 
u_RCS_vec = quat_rot_array(x(7:10, 1:Nu), u(2:4, :)); % Rotate 

figure
scatter3(0, 0, 0, 60, "blue", "filled", "diamond"); hold on
plot_cartesian_orbit(x_cont_sol(1:3,:)', 'k', 0.4, 1); hold on
quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u_main_vec(1, :), u_main_vec(2, :), u_main_vec(3, :), 1, "filled", Color = "red")
quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u_RCS_vec(1, :), u_RCS_vec(2, :), u_RCS_vec(3, :), 2, "filled", Color = "m")
scatter3(r_0(1), r_0(2), r_0(3), 48, "green", "filled", "square"); hold on
scatter3(r_f(1), r_f(2), r_f(3), 48, "red", "x"); hold on
plot3(guess.x(1, :) * nd_scalar(1), guess.x(2, :) * nd_scalar(2), guess.x(3, :) * nd_scalar(3), Color = "green", LineStyle = "--");
[s_x, s_y, s_z] = sphere(128);
h = surfl(s_x * keep_out_distance, s_y * keep_out_distance, s_z * keep_out_distance); 
set(h, 'FaceAlpha', 0.5)
shading interp
title('Optimal Rendezvous')
xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
axis equal
legend("Target", 'Spacecraft', "", "Thrust 1", "Thrust 2", "Start", "End", "Guess", "Keepout Sphere", 'Location', 'northwest'); axis equal; grid on

figure
tiledlayout(2, 1)

nexttile
plot(t_cont_sol, x_cont_sol(7:10, :))
xlabel("Time")
ylabel("Quaternions")
title("Quaternion History")
grid on

nexttile
plot(t_cont_sol, rad2deg(x_cont_sol(11:13, :)))
xlabel("Time")
ylabel("Angular Velocities")
title("Angular Velocity History")
grid on

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
plot(t_cont_sol(1:end - (N - Nu)), quat_rot_array(x_cont_sol(7:10, 1:end - (N - Nu)), u_cont_sol(2:4,:)), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(2:4,:)), LineWidth=1)
title("Reaction Control Thrusters")
xlabel("Time")
ylabel("Force [N]")
legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
grid on

nexttile
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(5:7,:), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(5:7,:)), LineWidth=1)
title("Reaction Wheel Torque")
xlabel("Time")
ylabel("Torque [N m]")
legend("\hat{b}_1", "\hat{b}_2", "\hat{b}_3", "||\tau||", Interpreter="latex")
grid on

%% Helper Functions

function [q] = qexp(tau)
    theta = norm(tau);
    u = tau / theta;
    q = [u * sin(theta / 2); cos(theta / 2)];
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
