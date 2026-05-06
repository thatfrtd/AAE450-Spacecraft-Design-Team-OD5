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
nd_scalar = [ones([6, 1]); char_star.m];
g_0 = 9.81; % [m / s2]

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = [4100; 300]; % [s]
spacecraft_params.m_0 = 1500; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = [0.25; 88]; % [N]
spacecraft_params.r = 1; % [m]
spacecraft_params.thrust_directions = {[1; 0; 0], ... % Main thruster
                                       [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
                                        0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0;
                                        0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]}; % RCS thruster directions
spacecraft_params.thrust_positions = {[1; 0; 0], ... % Main thruster
                                      [0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0;
                                       1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1;
                                       0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0] * spacecraft_params.r}; % RCS thruster positions
spacecraft_params.R_C_S = calculate_RCS_allocation_matrix(spacecraft_params.thrust_directions{2}, spacecraft_params.thrust_positions{2}); % RCS thruster force and torque allocation matrix
spacecraft_params.tau_max = 0; % [N m] Max reaction wheel torque
F_max_nd = spacecraft_params.F_max; % F_max in N, char_star.F in kN

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = char_star.l + 800; % [km] semi-major axis
e_c = 0.003; % [] eccentricity
i_c = deg2rad(98.6); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu0_c = deg2rad(0); % [rad] true anomaly at epoch
M0_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu0_c, e_c), e_c);
x_keplerian_c = [a_c; e_c; i_c; Omega_c; omega_c; M0_c];
n_c = sqrt(char_star.mu / a_c ^ 3);

% Passive Safety Parameters
P_c = 2 * pi * sqrt(a_c ^ 3 / char_star.mu); % Chief orbital period
T = 2 * P_c; % [s] Safety horizon length
N_safe = 200; 

% Rendezvous time
tf = P_c; % [s] (nondimensionalized)

% Initial conditions for spacecraft - specify orbit instead?
b = sqrt(40^2/2);
c = sqrt(20^2/2);
aROE_0 = [0; % [km] delta semimajor axis
          0; % [km] delta lambda
         -b*1e-3; % [km] delta e_x
          b*1e-3; % [km] delta e_y
         -c*1e-3; % [km] delta i_x 
          c*1e-3]; % [km] delta i_y

x_0 = [ROE_to_cart_matrix(n_c, 0 * char_star.t) * aROE_0; spacecraft_params.m_0] ./ nd_scalar;

% Terminal conditions
% aROE_f = [0; % [km] delta semimajor axis
%          0; % [km] delta lambda
%         -0.007954951288348661; % [km] delta e_x
%          0.007954951288348661; % [km] delta e_y
%         -0.007954951288348661; % [km] delta i_x 
%          0.007954951288348661]*1.5; % [km] delta i_y
% 
% x_f = ROE_to_cart_matrix(n_c, 0 * char_star.t) * aROE_f ./ nd_scalar(1:6);
data_table = readtable("Full_Detumble_History.csv");
x_f = [data_table.chaser_x(1); data_table.chaser_y(1); data_table.chaser_z(1); data_table.chaser_vx(1); data_table.chaser_vy(1); data_table.chaser_vz(1)] * 1e-3./ nd_scalar(1:6);

%% Initialize
N = 100;
t_k_actual = linspace(0, tf, N);
tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

% Control discretization method
% Zero Order Hold (ZOH) - Piecewise constant
% First Order Hold (FOH) - Piecewise linear (DOESN'T LIKE THIS PROBLEM... PROBABLY SCALING ISSUE :( )
u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

parser = "CVX"; % Can use CVXPY for more speed if needed (needs setup)
nx = 7; % Number of states
nu = 6; % Number of controls
np = 0; % Number of parameters (tf, v_0, etc)

% PTR algorithm parameters
ptr_ops.iter_max = 25;
ptr_ops.iter_min = 4;
ptr_ops.Delta_min = 1e-7;
ptr_ops.w_vc = 5e5;
ptr_ops.w_tr = ones(1, N) * 5e-2;
ptr_ops.w_tr_p = 0;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 0;
ptr_ops.alpha_p = 0;

% Scaling currently not helping...
scale = false;
% 
scale_hint.x_max = [ones([6, 1]); spacecraft_params.m_0 / char_star.m];
scale_hint.x_min = [zeros([7, 1])];
scale_hint.u_max = [ones([6, 1])];
scale_hint.u_min = [zeros([6, 1])];
scale_hint.p_max = [];
scale_hint.p_min = [];

%% Get Dynamics
f_nonlinear = @(t, x, u, p) nonlinear_relative_orbit_EoM_twothruster_dim(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp; char_star.mu; g_0]);

f_opt = f_nonlinear; % Dynamics to use for optimization
f_eval = f_nonlinear; % Dynamics to use for propagation and plotting

A_func = dynamics_jacobian(f_opt, nx, nu, np);

%% Specify Constraints
state_convex_constraints = {};

% Convex control constraints
max_thrust_constraint_1 = {1:N, @(t, x, u, p) norm(u(1:3)) - F_max_nd(1)};
max_thrust_constraint_2 = {1:N, @(t, x, u, p) norm(u(4:6)) - F_max_nd(2)};
control_convex_constraints = {max_thrust_constraint_1, max_thrust_constraint_2};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
keep_out_distance = 0.015; % [km]
keep_out_sphere_constraint = @(t, x, u, p) keep_out_distance ^ 2 - (x(1) ^ 2 + x(2) ^ 2 + x(3) ^ 2);
keep_out_sphere_constraint_linearized_func = linearize_constraint(keep_out_sphere_constraint, nx, nu, np, "x", 1:3);
keep_out_sphere_constraint_linearized = {1:N, linearize_constraint(keep_out_sphere_constraint, nx, nu, np, "x", 1:3)};
passive_safety_constraint = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) construct_passive_safety_constraint(x, x_ref(:, k), @(t_prop, x_prop) f_opt(t_prop + t, x_prop, zeros([nu, 1]), p_ref), @(t_prop, x_prop) A_func(t_prop + t, x_prop, zeros([nu, 1]), p_ref), @(t_prop, x_prop) keep_out_sphere_constraint(t_prop + t, x_prop, zeros([nu, 1]), p_ref), @(t_prop, x_prop, x_ref_safe) keep_out_sphere_constraint_linearized_func(t_prop + t, x_prop, zeros([nu, 1]), p_ref, x_ref_safe, zeros([nu, 1]), p_ref, 1), T, N_safe, 1e-5, -1e-5)};
state_nonconvex_constraints = {keep_out_sphere_constraint_linearized, passive_safety_constraint};

% Nonconvex control constraints
control_nonconvex_constraints = {};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Boundary conditions
initial_bc = @(x, p) [x - x_0];
terminal_bc = @(x, p, x_ref, p_ref) [x(1:6) - x_f; 0]; % Don't constrain final mass

%% Specify Objective
objective_min_fuel = @(x, u, p, x_ref, u_ref, p_ref) sum(norms(u(1:3, :))) * delta_t / (spacecraft_params.Isp(1) * g_0)*1000 ...
                                                   + sum(norms(u(4:6, :))) * delta_t / (spacecraft_params.Isp(2) * g_0);

%% Create Guess
% Straight Line Initial Guess
guess.x = linspace(0, 1, N) .* ([x_f; x_0(7)] - x_0) + x_0;
guess.u = ones([6, Nu]) * 1e-6;
guess.p = [];

%% Construct Problem Object
problem = DeterministicProblem(x_0, x_f, N, u_hold, tf, f_opt, guess, convex_constraints, objective_min_fuel, scale = scale, scale_hint = scale_hint, nonconvex_constraints = nonconvex_constraints, initial_bc = initial_bc, terminal_bc = terminal_bc, integration_tolerance = 1e-12, discretization_method = "error", N_sub = 1, Name = "nonlinear_rendezvous_twothruster");

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
u = ptr_sol.u(:, :, i);

problem.cont.f = f_eval;
[t_cont_sol, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));
x_cont_sol = x_cont_sol .* nd_scalar;
u_cont_sol = u_cont_sol;

%% Get Trajectories Checked for Passive Safety
N_safe_ck = 200;
min_safety = zeros([N, 1]);
x_safety_ck = zeros([nx, N_safe_ck, N]);
for k = 1 : N
    [~, min_safety(k), x_safety_ck(:, :, k)] = construct_passive_safety_constraint(ptr_sol.x(:, k, i), ptr_sol.x(:, k, i), @(t_prop, x_prop) f_opt(t_prop + t_k(k), x_prop, zeros([nu, 1]), []), @(t_prop, x_prop) A_func(t_prop + t_k(k), x_prop, zeros([nu, 1]), []), @(t_prop, x_prop) keep_out_sphere_constraint(t_prop + t_k(k), x_prop, zeros([nu, 1]), []), @(t_prop, x_prop, x_ref_safe) keep_out_sphere_constraint_linearized_func(t_prop + t_k(k), x_prop, zeros([nu, 1]), [], x_ref_safe, zeros([nu, 1]), [], 1), T, N_safe_ck, 1e-6);
end
x_safety_ck = x_safety_ck .* nd_scalar;

%% Recover RCS thruster inputs
q_body = zeros([4, N]);
u_RCS = zeros([size(spacecraft_params.R_C_S, 2), N]);
for k = 1 : N
    q_body(:, k) = point_at_vec(-x(1:3, k));
    if k <= Nu
        F_body = quat_rot(q_conj(q_body(:, k)), u(4:6, k));
        u_RCS(:, k) = pinv(spacecraft_params.R_C_S) * [F_body; 0; 0; 0];
    end
end

%% Save Drift Trajectories
output_array = [reshape(x_safety_ck(1:3, :, :), [], 1)];
state_names = ["drift_r"];
output_names = [state_names];

output_table = array2table(output_array, VariableNames = output_names);
writetable(output_table,"./Animation/transfer_inspection_to_detumble_drifts.csv")

output_array = [q_body', u_RCS'];
state_names = ["q_1", "q_2", "q_3", "q_4"];
control_names = "";
for i = 1 : size(spacecraft_params.thrust_directions{2}, 2)
    control_names(i) = sprintf("u_RCS_%g", i);
end
output_names = [state_names, control_names];
output_table = array2table(output_array, VariableNames = output_names);
writetable(output_table,"./Animation/transfer_inspection_to_detumble_attitude_controls.csv")


%% Plot Trajectory
figure;
fig = scatter3(0, 0, 0, 60, "blue", "filled", "diamond"); hold on
lim = max(abs(x(1:3, :)), [], "all");
for k = (N - N + 1):(N - 0)
    if k == 1
        handvis = "on";
    else
        handvis = "off";
    end
    x_safety_ck(:, vecnorm(x_safety_ck(1:3, :, k)) > lim, k) = nan;
    plot3(x_safety_ck(1, :, k), x_safety_ck(2, :, k), x_safety_ck(3, :, k), Color="b", HandleVisibility=handvis)
end
for k = 1 : Nu
    u_norm = norm(u(1:3, :));
    if u_norm < 1e-10
        u(1:3, k) = 0;
    end
end

plot_cartesian_orbit(x_cont_sol(1:3,:)', 'k', 0.4, 1, 1.5); hold on
%quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(1, :), u(2, :), u(3, :), 1, "filled", Color = "red")

quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(4, :), u(5, :), u(6, :), 2, "filled", Color = "m", LineWidth=2)
scatter3(x_0(1) * nd_scalar(1), x_0(2) * nd_scalar(2), x_0(3) * nd_scalar(3), 48, "green", "filled", "square"); hold on
scatter3(x_f(1) * nd_scalar(1), x_f(2) * nd_scalar(2), x_f(3) * nd_scalar(3), 48, "red", "x"); hold on
plot3(guess.x(1, :), guess.x(2, :), guess.x(3, :), Color = "green", LineStyle = "--");
[s_x, s_y, s_z] = sphere(128);
h = surfl(s_x * keep_out_distance, s_y * keep_out_distance, s_z * keep_out_distance); 
set(h, 'FaceAlpha', 0.5)
shading interp
plot_cartesian_orbit(r_cont_sol_inspection', 'g', 0.4, 1, 1); hold on
plot_cartesian_orbit([data_table.chaser_x(:), data_table.chaser_y(:), data_table.chaser_z(:)]*1e-3, 'r', 0.4, 1, 1); hold on

title('Optimal Passively Safe Relative Transfer')

xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
axis equal

legend("Target", "Drift", 'Spacecraft', "", "RCS Thrust", "Start", "End", "Guess", "Keepout Sphere", "", "Inspection eNMT", "", "Detumble eNMT", 'Location', 'northwest'); axis equal; grid on


%% Plot Control
figure
tiledlayout(2, 1)

nexttile
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(1:3,:), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(1:3,:)), LineWidth=1)
title("Control Thruster 1")
xlabel("Time")
ylabel("Force [N]")
legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
grid on

nexttile
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(4:6,:), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(4:6,:)), LineWidth=1)
title("Control Thruster 2")
xlabel("Time")
ylabel("Force [N]")
legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
grid on

%%
figure
plot(t_cont_sol(1:end - (N - Nu)) / 60 / 60, u_cont_sol(4:6,:), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)) / 60 / 60, vecnorm(u_cont_sol(4:6,:)), LineWidth=1)
xlim([0, tf / 60 / 60])
title("RCS Thrust for Safe Relative Orbit Transfer")
xlabel("Time [hr]")
ylabel("Force [N]")
legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
grid on

%% Constraint Satisfaction
figure
plot(t_cont_sol / 60 / 60, vecnorm(x_cont_sol(1:3, :)))
yline(keep_out_distance)
xlabel("Time [hr]")
ylabel("Distance [km]")
title("Distance from Target vs Time")
legend("Solution", "Minimum")
grid on

%% Save Trajectory
r_cont_sol = x(1:3, :);

output_array = [t_k', r_cont_sol'];
state_names = ["r_1", "r_2", "r_3"];
output_names = ["Time", state_names];

output_table = array2table(output_array, VariableNames = output_names);
writetable(output_table,"./Animation/transfer_inspection_to_detumble.csv")

%% Plot Iterations
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
[s_x, s_y, s_z] = sphere(128);
h = surfl(s_x * keep_out_distance, s_y * keep_out_distance, s_z * keep_out_distance); 
set(h, 'FaceAlpha', 0.5)
shading interp

title('SCP Iteration Trajectories for Passive Rendezvous')

xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
axis equal

c = colorbar();
c.Label.String = 'Iteration';
c.TickLabels = string(0:ptr_sol.converged_i);

%% Helper Functions
function [A] = dynamics_jacobian(f, nx, nu, np)

    t_sym = sym("t");
    x_sym = sym("x", [nx, 1]);
    u_sym = sym("u", [nu, 1]);
    p_sym = sym("p", [np, 1]);
    
    % Linearize Dynamics
    A = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
end


function [quat] = point_at_vec(T_vec)
    T_vec = T_vec ./ vecnorm(T_vec);
    w = 1 + T_vec(1, :);
    xyz = cross([ones([1, size(T_vec, 2)]); zeros([2, size(T_vec, 2)])], T_vec); % cross b_x with T_vec
    quat = [xyz; w];
    quat = quat ./ vecnorm(quat);
end