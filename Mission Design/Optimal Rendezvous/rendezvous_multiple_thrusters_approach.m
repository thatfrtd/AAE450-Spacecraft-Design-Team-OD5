%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Nonlinear Rendezvous Example
% Author: Travis Hastreiter 
% Created On: 3 April, 2026
% Description: Sequential Convex Programming Trajectory Optimization for 
% general rendezvous using relative orbital motion equations. Includes mass 
% in the state. You must have CVX installed. Multiple thrusters with
% different thrust levels and Isp but using the same fuel
% (resistojet/arcjet eprop with chemical RCS both using hydrazine).
% Last Modified On: 3 April, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

char_star = load_charecteristic_values_Earth();
nd_scalar = [char_star.l * ones([3, 1]); char_star.v * ones([3, 1]); char_star.m];
g_0 = 9.81e-3; % [km / s2]

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.r = 1; % [m]
spacecraft_params.h = 2; % [m]
spacecraft_params.Isp = [4155; 303]; % [s]
spacecraft_params.m_0 = 1500; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = [0; 88]; % [N]
F_max_nd = spacecraft_params.F_max / 1000 / char_star.F; % F_max in N, char_star.F in kN

% Estimate spacecraft moment of inertia
cylinder_MoI = @(m, r, h) [1 / 2 * m * r ^ 2; ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2); ...
                           1 / 12 * m * (3 * r ^ 2 + h ^ 2)];

% Spacedebris Parameters
spacedebris_params = struct();
spacedebris_params.Isp = [4100; 232]; % [s]
spacedebris_params.m = 4000; % [kg]
spacedebris_params.r = 1.85; % [m]
spacedebris_params.h = 11.9; % [m]
spacedebris_params.I = cylinder_MoI(spacedebris_params.m, spacedebris_params.r, spacedebris_params.h);

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
tf = 60 / char_star.t; % [s] (nondimensionalized)
N = 100;
t_k_actual = linspace(0, tf, N);
tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

% Initial conditions for spacecraft - specify orbit instead?
data_table = readtable("Full_Detumble_History.csv");
rv_0 = [data_table.chaser_x(end); data_table.chaser_y(end); data_table.chaser_z(end); data_table.chaser_vx(end); data_table.chaser_vy(end); data_table.chaser_vz(end)] * 1e-3;
r_0 = rv_0(1:3);
v_0 = rv_0(4:6);
x_0 = [rv_0; spacecraft_params.m_0] ./ nd_scalar;

% Terminal conditions
seperation_distance = 5; % [m]
d_vec = [-(spacedebris_params.h / 2 + spacecraft_params.h / 2 + seperation_distance); 0; 0]*1e-3;
w0_debris = -deg2rad([0.0739732221634118; 0.0063002539406854; 0.0109314217396326]);
q0_debris = [data_table.q_x(end); data_table.q_y(end); data_table.q_z(end); data_table.q_w(end)];
[x_f_interstellar, qw_debris_array] = interstellar_boundary_conditions(q0_debris, w0_debris, t_k, spacedebris_params.I, d_vec);
x_f = x_f_interstellar(1:6) ./ nd_scalar(1:6);
q_f = x_f_interstellar(7:10);
r_f = x_f_interstellar(1:3);

%% Initialize
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
ptr_ops.iter_max = 15;
ptr_ops.iter_min = 4;
ptr_ops.Delta_min = 1e-7;
ptr_ops.w_vc = 5e5;
ptr_ops.w_tr = ones(1, N) * 5e2;
ptr_ops.w_tr_p = 0;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
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
f_nonlinear = @(t, x, u, p) nonlinear_relative_orbit_EoM_twothruster(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp; 1; g_0]);

f_opt = f_nonlinear; % Dynamics to use for optimization
f_eval = f_nonlinear; % Dynamics to use for propagation and plotting

%% Specify Constraints
% Convex state constraints
state_convex_constraints = {};

% Convex control constraints
max_thrust_constraint_1 = {1:N, @(t, x, u, p) norm(u(1:3)) - F_max_nd(1)};
max_thrust_constraint_2 = {1:N, @(t, x, u, p) norm(u(4:6)) - F_max_nd(2)};
control_convex_constraints = {max_thrust_constraint_1, max_thrust_constraint_2};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
approach_cone_angle = deg2rad(80);
d_capture = 11e-3 ./ char_star.l; % [m] Distance of capture point (where s/c CoM will be) from debris CoM
approach_cone_constraint = @(t, x, u, p) cos(approach_cone_angle) * norm(x(1:3) - [-d_capture; 0; 0]) + (x(1) + d_capture);
keep_out_ellipsoid_form = [18; 10; 10]*1e-3 ./ nd_scalar(1:3); % [km] 
keep_out_ellipsoid_constraint = @(t, x, u, p) 1 - ((x(1) / keep_out_ellipsoid_form(1)) ^ 2 + (x(2) / keep_out_ellipsoid_form(2)) ^ 2 + (x(3) / keep_out_ellipsoid_form(3)) ^ 2);
keep_out_ellipsoid_constraint_linearized_func = linearize_constraint(keep_out_ellipsoid_constraint, nx, nu, np, "x", 1:3);
keep_out_ellipsoid_constraint_linearized_STC = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) keep_out_ellipsoid_constraint_linearized_func(t, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x(1:3), u, p, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x_ref(1:3, :), u_ref, p_ref, k) * max(approach_cone_constraint(t, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x_ref(1:3, k), u_ref(:, k), p_ref), 0)*1e5};
%keep_out_ellipsoid_constraint_linearized_STC = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) keep_out_ellipsoid_constraint_linearized_func(t, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x(1:3), u, p, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x_ref(1:3, :), u_ref, p_ref, k)*1e5};
approach_cone_STC = {1:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) approach_cone_constraint(t, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x(1:3), u, p) * max(keep_out_ellipsoid_constraint(t, quat_rotmatrix(q_conj(qw_debris_array(1:4, k))) * x_ref(1:3, k), u_ref(:, k), p_ref), 0)*1e3};% * (keep_out_ellipsoid_constraint(t, x_ref(:, k), u_ref(:, k), p_ref) > 0)};
state_nonconvex_constraints = {keep_out_ellipsoid_constraint_linearized_STC, approach_cone_STC};

%%
approach_cone_constraint([], quat_rotmatrix(q_conj(qw_debris_array(1:4, end))) * x_f(1:3), [], [])
keep_out_ellipsoid_constraint_linearized_STC{2}(0, x_f .* [1.2 * ones([3, 1]); ones([3, 1])], [], [], x_f, 0, [], 1)
approach_cone_STC{2}(0, x_f .* [1.2 * ones([3, 1]); ones([3, 1])], [], [], x_f, 0, [], 1)
%%

% Nonconvex control constraints
control_nonconvex_constraints = {};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Boundary conditions
initial_bc = @(x, p) [x - x_0];
terminal_bc = @(x, p, x_ref, p_ref) [x(1:6) - x_f; 0] .* nd_scalar; % Don't constrain final mass

%% Specify Objective
objective_min_fuel = @(x, u, p, x_ref, u_ref, p_ref) sum(norms(u(1:3, :))) * delta_t * char_star.F / (spacecraft_params.Isp(1) * g_0) * 1000 ...
                                                   + sum(norms(u(4:6, :))) * delta_t * char_star.F / (spacecraft_params.Isp(2) * g_0) * 1000;

%% Create Guess
% Straight Line Initial Guess
guess.x = linspace(0, 1, N) .* ([x_f; x_0(7)] - x_0) + x_0;
guess.u = ones([nu, Nu]) * 1e-6;
guess.p = [];

%% Construct Problem Object
problem = DeterministicProblem(x_0, x_f, N, u_hold, tf, f_opt, guess, convex_constraints, objective_min_fuel, scale = scale, nonconvex_constraints = nonconvex_constraints, initial_bc = initial_bc, terminal_bc = terminal_bc, integration_tolerance = 1e-12, discretization_method = "error", N_sub = 1, Name = "nonlinear_rendezvous_twothruster");

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

problem.cont.f = f_eval;
[t_cont_sol, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));
t_cont_sol = t_cont_sol * char_star.t;
x_cont_sol = x_cont_sol .* nd_scalar;
u_cont_sol = u_cont_sol * char_star.F * 1000;

%% Plot Trajectory
figure
for k = 1 : Nu
    u_norm = norm(u(1:3, :));
    if u_norm < 1e-10
        u(1:3, k) = 0;
    end
end
scatter3(0, 0, 0, 60, "blue", "filled", "diamond"); hold on
plot_cartesian_orbit(x_cont_sol(1:3,:)', 'k', 0.4, 1); hold on
quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(1, :), u(2, :), u(3, :), 1, "filled", Color = "red")
quiver3(x(1, 1:Nu), x(2, 1:Nu), x(3, 1:Nu), u(4, :), u(5, :), u(6, :), 2, "filled", Color = "m")
scatter3(r_0(1), r_0(2), r_0(3), 48, "green", "filled", "square"); hold on
scatter3(r_f(1), r_f(2), r_f(3), 48, "red", "x"); hold on
plot3(guess.x(1, :) * nd_scalar(1), guess.x(2, :) * nd_scalar(2), guess.x(3, :) * nd_scalar(3), Color = "green", LineStyle = "--");
[x_1, y_1, z_1] = ellipsoid(0,0,0,keep_out_ellipsoid_form(1) * char_star.l,keep_out_ellipsoid_form(2) * char_star.l,keep_out_ellipsoid_form(3) * char_star.l);
h = surf(x_1, y_1, z_1);
tau = qLog(q_f);
theta = norm(tau);
if theta > 1e-10
    lambda = tau / theta;
else
    lambda = [1; 0; 0];
end
rotate(h, lambda, -rad2deg(theta), [0; 0; 0]);
set(h, 'FaceAlpha', 0.4)
[X,Y,Z]=cylinder([0 0.5],50 );
axis([0 1,-1 1,-.5 .5])
M=makehgtform('translate',[0,0,0],'xrotate',pi/4,'yrotate',pi/2);
h=surf(X*30e-3,Y*30e-3,Z*15e-3 + d_capture * char_star.l,'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.4);
view([30,35])
grid on
light
shading interp
title('Optimal Rendezvous')
xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
axis equal
legend("Target", 'Spacecraft', "", "Thrust 1", "Thrust 2", "Start", "End", "Guess", "Keepout Sphere", 'Location', 'northwest'); axis equal; grid on



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
