function [rendezvous_solution, ptr_sol] = nonlinear_rendezvous_func(x_0_hill, x_f_hill, ToF, x_keplerian_c, spacecraft_params, options)
arguments
    x_0_hill
    x_f_hill
    ToF
    x_keplerian_c
    spacecraft_params
    options.N = 100
    options.u_hold = "FOH"
    options.dynamics = "Nonlinear" % "CWH", "Linearized", "Nonlinear"
    options.plot_results = true
    options.plot_convergence = false
end
char_star = load_charecteristic_values_Earth();
nd_scalar = [char_star.l * ones([3, 1]); char_star.v * ones([3, 1]); char_star.m];

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
F_max_nd = spacecraft_params.F_max / 1000 / char_star.F; % F_max in N, char_star.F in kN

% Rendezvous time
tf = ToF / char_star.t; % [s] (nondimensionalized)

% Initial conditions for spacecraft - specify orbit instead?
x_0_nd = x_0_hill ./ nd_scalar;

% Terminal conditions
x_f_nd = x_f_hill ./ nd_scalar(1:6);

%% Initialize
N = options.N;
t_k_actual = linspace(0, tf, N);
tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

% Control discretization method
% Zero Order Hold (ZOH) - Piecewise constant
% First Order Hold (FOH) - Piecewise linear
u_hold = options.u_hold;
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

parser = "CVX"; % Can use CVXPY for more speed if needed (needs setup)
nx = 7; % Number of states
nu = 3; % Number of controls
np = 0; % Number of parameters (tf, v_0, etc)

% PTR algorithm parameters
ptr_ops.iter_max = 5;
ptr_ops.iter_min = 1;
ptr_ops.Delta_min = 1e-8;
ptr_ops.w_vc = 5e5;
ptr_ops.w_tr = ones(1, Nu) * 5e-2;
ptr_ops.w_tr_p = 0;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

% Scaling currently not helping...
scale = false;

scale_hint.x_max = [max(x_0_nd(1:3)) * ones([3, 1]); max(x_0_nd(4:6)) * ones([3, 1]); spacecraft_params.m_0 / char_star.m];
scale_hint.x_min = [-max(x_0_nd(1:3)) * ones([3, 1]); -max(x_0_nd(4:6)) * ones([3, 1]); spacecraft_params.m_0 / char_star.m * 0.95];
scale_hint.u_max = [F_max_nd * ones([3, 1])];
scale_hint.u_min = [zeros([3, 1])];
scale_hint.p_max = [];
scale_hint.p_min = [];

%% Get Dynamics
f_nonlinear = @(t, x, u, p) nonlinear_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
f_linearized = @(t, x, u, p) linearized_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
f_CWH = @(t, x, u, p) CWH_relative_orbit_EoM(t, x, u, p, [x_keplerian_c(1); spacecraft_params.Isp]);

if options.dynamics == "CWH"
    f_opt = f_CWH;
elseif options.dynamics == "Linearized"
    f_opt = f_linearized;
elseif options.dynamics == "Nonlinear"
    f_opt = f_nonlinear;
end
f_eval = f_nonlinear; % Dynamics to use for propagation and plotting

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
initial_bc = @(x, p) [x - x_0_nd];
terminal_bc = @(x, p, x_ref, p_ref) [x(1:6) - x_f_nd; 0]; % Don't constrain final mass

%% Specify Objective
objective_min_fuel = @(x, u, p, x_ref, u_ref, p_ref) sum(norms(u)) * delta_t * char_star.F / (spacecraft_params.Isp * 9.81e-3) * 1000;

%% Create Guess
% Straight Line Initial Guess - Lambert better?
guess.x = linspace(0, 1, N) .* (x_0_nd - [x_f_nd; x_0_nd(7)]) + x_0_nd;
guess.u = ones([3, Nu]) * 1e-6;
guess.p = [];

%% Construct Problem Object
problem = DeterministicProblem(x_0_nd, x_f_nd, N, u_hold, tf, f_opt, guess, convex_constraints, objective_min_fuel, scale = scale, nonconvex_constraints = nonconvex_constraints, initial_bc = initial_bc, terminal_bc = terminal_bc, integration_tolerance = 1e-12, discretization_method = "error", N_sub = 1, Name = "nonlinear_rendezvous");

[problem, Delta_disc] = problem.discretize(guess.x, guess.u, guess.p);

%% Solve
ptr_sol = ptr(problem, ptr_ops, parser);

%% Plot Convergence
if ptr_sol.converged == false
    ptr_sol.converged_i = ptr_ops.iter_max;
end

if options.plot_convergence
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
end

%% Extract Solution
i = ptr_sol.converged_i + 1;
x = ptr_sol.x(:, :, i) .* nd_scalar;
u = ptr_sol.u(:, :, i) * char_star.F * 1000;

problem.cont.f = f_eval;
[t_cont_sol, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));
t_cont_sol = t_cont_sol * char_star.t;
x_cont_sol = x_cont_sol .* nd_scalar;
u_cont_sol = u_cont_sol * char_star.F * 1000;

%% Package Output
rendezvous_solution.t = t_cont_sol;
rendezvous_solution.x = x_cont_sol;
rendezvous_solution.u = u_cont_sol;

%% Plot
if options.plot_results
    % Plot Trajectory
    figure
    scatter3(0, 0, 0, 60, "blue", "filled", "diamond"); hold on
    plot_cartesian_orbit(x_cont_sol(1:3,:)', 'k', 0.4, 1); hold on
    quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(1, :), u(2, :), u(3, :), 1, "filled", Color = "red")
    scatter3(x_0_hill(1), x_0_hill(2), x_0_hill(3), 48, "green", "filled", "square"); hold on
    scatter3(x_f_hill(1), x_f_hill(2), x_f_hill(3), 48, "red", "x"); hold off
    title('Optimal Rendezvous')
    xlabel("r [km]")
    ylabel("\theta [km]")
    zlabel("n [km]")
    legend("Target", 'Spacecraft', "", "Thrust", "Start", "End", 'Location', 'northwest'); axis equal; grid on
    
    % Plot Control
    figure
    
    plot(t_cont_sol(1:end - (N - Nu)), u_cont_sol(1:3,:), LineWidth=1); hold on
    plot(t_cont_sol(1:end - (N - Nu)), vecnorm(u_cont_sol(1:3,:)), LineWidth=1)
    title("Control")
    xlabel("Time")
    ylabel("Force [N]")
    legend("\hat{r}", "\hat{\theta}", "\hat{h}", "||u||", Interpreter="latex")
    grid on
end
end