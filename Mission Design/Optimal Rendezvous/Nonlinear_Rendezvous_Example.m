%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Nonlinear Rendezvous Example
% Author: Travis Hastreiter 
% Created On: 11 February, 2026
% Description: Sequential Convex Programming Trajectory Optimization for 
% general rendezvous using relative orbital motion equations. Includes mass 
% in the state. You must have CVX installed.
% Created On: 11 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

char_star = load_charecteristic_values_Earth();

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 1000; % [s]
spacecraft_params.m_0 = 800; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.5; % [N]
F_max_nd = spacecraft_params.F_max / 1000 / char_star.F; % F_max in N, char_star.F in kN

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = 6728; % [km] semi-major axis
e_c = 0; % [] eccentricity
i_c = deg2rad(10); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu0_c = deg2rad(0); % [rad] true anomaly at epoch
M0_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu0_c, e_c), e_c);
x_keplerian_c = [a_c; e_c; i_c; Omega_c; omega_c; M0_c];

% Initial conditions for spacecraft - specify orbit instead?
r_0 = [-1.5; -0.5; 0.2]; % [km]
v_0 = [3e-3; 0; 0]; % [km / s]
x_0 = [r_0 / char_star.l; v_0 / char_star.v; spacecraft_params.m_0 / char_star.m];

% Terminal conditions
r_f = [-0.100; 0; 0]; % [km]
v_f = [0.2e-3; 0; 0]; % [km / s]
x_f = [r_f / char_star.l; v_f / char_star.v];

%% Initialize
tf = 3600*5 / char_star.t;
N = 50;
t_k_actual = linspace(0, tf, N);
tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

parser = "CVX";
nx = 7;
nu = 3;
np = 0;

% PTR algorithm parameters
ptr_ops.iter_max = 10;
ptr_ops.iter_min = 1;
ptr_ops.Delta_min = 5e-10;
ptr_ops.w_vc = 5e4;
ptr_ops.w_tr = ones(1, Nu) * 5e-2;
ptr_ops.w_tr_p = 0;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = false;

scale_hint.x_max = [max(x_0(1:3)) * ones([3, 1]); max(x_0(4:6)) * ones([3, 1]); spacecraft_params.m_0 / char_star.m];
scale_hint.x_min = [-max(x_0(1:3)) * ones([3, 1]); -max(x_0(4:6)) * ones([3, 1]); spacecraft_params.m_0 / char_star.m * 0.95];
scale_hint.u_max = [F_max_nd * ones([3, 1])];
scale_hint.u_min = [zeros([3, 1])];
scale_hint.p_max = [];
scale_hint.p_min = [];

%% Get Dynamics
%f = @(t, x, u, p) nonlinear_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
%f = @(t, x, u, p) linearized_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
f = @(t, x, u, p) CWH_relative_orbit_EoM(t, x, u, p, [a_c; spacecraft_params.Isp]);

%% Specify Constraints
state_convex_constraints = {};

% Convex control constraints
% min_periapsis_constraint = 
max_thrust_constraint = {1:N, @(t, x, u, p) norm(u) - F_max_nd};
control_convex_constraints = {max_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
state_nonconvex_constraints = {};

% Nonconvex control constraints
control_nonconvex_constraints = {};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Boundary conditions
initial_bc = @(x, p) x - x_0;
terminal_bc = @(x, p, x_ref, p_ref) [x(1:6) - x_f; 0]; % Don't constrain final mass

%% Specify Objective
objective_min_fuel = @(x, u, p, x_ref, u_ref, p_ref) sum(norms(u)) * delta_t * char_star.F / (spacecraft_params.Isp * 9.81e-3) * 1000;

%% Create Guess
% Straight Line Initial Guess - Lambert better?
guess.x = linspace(0, 1, N) .* (x_0 - [x_f; x_0(7)]) + x_0;
guess.u = ones([3, Nu]) * 1e-6;
guess.p = [];

%% Construct Problem Object
problem = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, objective_min_fuel, scale = scale, nonconvex_constraints = nonconvex_constraints, initial_bc = initial_bc, terminal_bc = terminal_bc, integration_tolerance = 1e-12, discretization_method = "error", N_sub = 1, Name = "nonlinear_rendezvous");

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
x = ptr_sol.x(:, :, i) .* [char_star.l * ones(3, 1); char_star.v * ones(3, 1); char_star.m];
u = ptr_sol.u(:, :, i) * char_star.F * 1000;

[t_cont_sol, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));
t_cont_sol = t_cont_sol * char_star.t;
x_cont_sol = x_cont_sol .* [char_star.l * ones(3, 1); char_star.v * ones(3, 1); char_star.m];
u_cont_sol = u_cont_sol * char_star.F * 1000;

%% Plot Trajectory
figure
plot_cartesian_orbit(x_cont_sol(1:3,:)', 'k', 0.4, 1); hold on
quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(1, :), u(2, :), u(3, :), 1, "filled", Color = "red")
title('Optimal Rendezvous')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]')
legend('Spacecraft', "", "Thrust", 'Location', 'northwest'); axis equal; grid on

%% Plot Control
figure

plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(1:3,:), LineWidth=1); hold on
plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(1:3,:)), LineWidth=1)
title("Control")
xlabel("Time")
ylabel("Force [N]")
legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
grid on