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

% Estimate spacecraft moment of inertia
cylinder_MoI = @(m, r, h) [1 / 2 * m * r ^ 2; ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2); ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2)];

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 232; % [s]
spacecraft_params.m_0 = 1500; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.r = 1.85; % [m]
spacecraft_params.h = 4; % [m]
spacecraft_params.I = cylinder_MoI(spacecraft_params.m_0, spacecraft_params.r, spacecraft_params.h);
spacecraft_params.F_max = 5; % [N]
% RCS (direction axis, position axis) x+y, x-y, x+z, x-z, y+z, y-z, y+x, y-x, z+x, z-x, z+y, z-y
spacecraft_params.thrust_directions = [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0;
                                       0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0;
                                       0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]; % RCS thruster directions
spacecraft_params.thrust_positions = [0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0;
                                      1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1;
                                      0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 0, 0] * spacecraft_params.r; % RCS thruster positions
spacecraft_params.R_C_S = calculate_RCS_allocation_matrix(spacecraft_params.thrust_directions, spacecraft_params.thrust_positions); % RCS thruster force and torque allocation matrix
spacecraft_params.tau_max = 15; % [N m] Max reaction wheel torque
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
tf = 35; % [s] (nondimensionalized)

% Initial conditions for spacecraft - specify orbit instead?
theta_0 = deg2rad([0; 45; 45]); % [rad]
R_0 = angle2dcm(theta_0(1), theta_0(2), theta_0(3));
q_0 = qExp(RLog(R_0));
w_0 = deg2rad([0; 0; 0]); % [rad / s]
x_0 = [q_0; w_0];

% Terminal conditions
theta_f = deg2rad([0; 0; 0]); % [rad]
R_f = angle2dcm(theta_f(1), theta_f(2), theta_f(3));
q_f = qExp(RLog(R_f));
w_f = deg2rad([0; 10; 0]); % [rad / s]
x_f = [q_f; w_f];

%% Initialize
N = 100;
t_k_actual = linspace(0, tf, N);
tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

% Control discretization method
% Zero Order Hold (ZOH) - Piecewise constant
% First Order Hold (FOH) - Piecewise linear
u_hold = "FOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

parser = "CVX"; % Can use CVXPY for more speed if needed (needs setup)
nx = 7; % Number of states (q, w)
nu = 15; % Number of controls (tau_rw)
np = 0; % Number of parameters (tf, v_0, etc)

% PTR algorithm parameters
ptr_ops.iter_max = 25;
ptr_ops.iter_min = 7;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 5e5;
ptr_ops.w_tr = ones(1, Nu) * 5e-3;
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
% scale_hint.x_max = 2 * pi * ones([7, 1]);
% scale_hint.x_min = -2 * pi * ones([7, 1]);
% scale_hint.u_max = [tau_max_nd * ones([3, 1])];
% scale_hint.u_min = [-tau_max_nd * zeros([3, 1])];
% scale_hint.p_max = [];
% scale_hint.p_min = [];

%% Get Dynamics
f = @(t, x, u, p) rotational_3DoF_RCS_EoM(t, x, u, p, [spacecraft_params.I; spacecraft_params.R_C_S(:)]);

%% Specify Constraints
max_ang_vel = deg2rad(50);
angular_velocity_constraint = {1:N, @(t, x, u, p) norm(x(5:7), Inf) - norm([w_0; max_ang_vel], Inf)};
% sensor seeing target
state_convex_constraints = {angular_velocity_constraint};

% Convex control constraints
max_RCS_thrust_constraint = {1:N, @(t, x, u, p) norm(u(1:12), Inf) - spacecraft_params.F_max};
max_reaction_wheel_torque_constraint = {1:N, @(t, x, u, p) norm(u(13:15), Inf) - tau_max_nd};
control_convex_constraints = {max_RCS_thrust_constraint, max_reaction_wheel_torque_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
state_nonconvex_constraints = {};

% Nonconvex control constraints
% Plume impingement constraints
control_nonconvex_constraints = {};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Boundary conditions
initial_bc = @(x, p) [x - x_0];
terminal_bc = @(x, p, x_ref, p_ref) [x - x_f]; % Don't constrain final mass
%terminal_bc = @(x, p, x_ref, p_ref) [x([4:7]') - x_f([4:7]'); zeros([4, 1])]; % Don't constrain final mass

%% Specify Objective
objective_min_fuel = @(x, u, p, x_ref, u_ref, p_ref) sum(norms(u(1:12, :), 1)) * delta_t + sum(norms(u(13:15, :), 1)) * delta_t * 1e-3;

%% Create Guess
% Straight Line Initial Guess - Definitely need to actually slerp for
% proper quaternion interpolation
[q_int, w_int_const] = q_interp(q_0, q_f, t_k);
w_int = repmat(w_int_const, 1, numel(t_k));

guess.x = [q_int; w_int];
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
x = ptr_sol.x(:, :, i);
u = ptr_sol.u(:, :, i);

t_cont_sol = linspace(0, tf, 1000);
[~, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i), tspan = t_cont_sol);
%t_cont_sol = t_cont_sol;
%x_cont_sol = x_cont_sol;
%u_cont_sol = u_cont_sol;

%% Plot Trajectory
figure
tiledlayout(3, 1)

nexttile
plot(t_cont_sol, x_cont_sol(1:4, :))
xlabel("Time")
ylabel("Quaternions")
title("Quaternion History")
grid on

nexttile
plot(t_cont_sol, rad2deg(x_cont_sol(5:7, :)))
xlabel("Time")
ylabel("Angular Velocities")
title("Angular Velocity History")
grid on

% Plot Control
nexttile
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(13:15,:), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(13:15,:)), LineWidth=1)
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
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(1, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(2, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("x+y", "x-y")
title("X Directed Y Positioned")
grid on
ylim(spacecraft_params.F_max * [-1; 1])

nexttile % yz
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(5, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(6, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("y+z", "y-z")
title("Y Directed Z Positioned")
grid on
ylim(spacecraft_params.F_max * [-1; 1])

nexttile % zx
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(9, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(10, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("z+x", "z-x")
title("Z Directed X Positioned")
grid on
ylim(spacecraft_params.F_max * [-1; 1])

nexttile % xz
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(3, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(4, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("x+z", "x-z")
title("X Directed Z Positioned")
grid on
ylim(spacecraft_params.F_max * [-1; 1])

nexttile % yx
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(7, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(8, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("y+x", "y-x")
title("Y Directed X Positioned")
grid on
ylim(spacecraft_params.F_max * [-1; 1])

nexttile % zy
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(11, :), Color = "r", LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(12, :), Color = "b", LineWidth=1);
xlabel("Time [s]")
ylabel("Force [N]")
legend("z+y", "z-y")
title("Z Directed Y Positioned")
grid on
ylim(spacecraft_params.F_max * [-1; 1])

sgtitle("Reaction Control System Thrusters Control History")

%% Save to CSV for Blender Animation
output_array = [t_cont_sol', [zeros([3, size(x_cont_sol, 2)]); x_cont_sol]', u_cont_sol'];
state_names = ["r_Hill_1", "r_Hill_2", "r_Hill_3", ...
               "q_1", "q_2", "q_3", "q_4", ...
               "w_1", "w_2", "w_3"];
control_names = strings(1, nu);
for i = 1 : size(spacecraft_params.thrust_directions, 2)
    control_names(i) = sprintf("u_RCS_%g", i);
end
for i = 1 : 3
    control_names(i + size(spacecraft_params.thrust_directions, 2)) = sprintf("u_W_%g", i);
end
output_names = ["Time", state_names, control_names];

output_table = array2table(output_array, VariableNames = output_names);
writetable(output_table,"./Animation/3DoF_test_animation_output.csv")

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
    tau = zeros([3, N]);
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