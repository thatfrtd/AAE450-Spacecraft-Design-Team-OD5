%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% CWH Rendezvous Example
% Author: Travis Hastreiter 
% Created On: 11 February, 2026
% Description: CWH Convex trajectory optimization with ZOH control. Computes 
% a fuel optimal relative orbit transfer. You must have CVX installed
% Created On: 12 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
% Physical parameters
mu = 398600.4415; % [km2 / s3] Earth gravitational parameter
u_max = 5e-6; % [km / s^2] Max control acceleration

% Initial conditions
r_0 = [-1.5; -0.5; 0.2]; % [km]
v_0 = [3e-3; 0; 0]; % [km / s]
x_0 = [r_0; v_0];

% Terminal conditions
r_f = [-0.25; 0; 0]; % [km]
v_f = [0.2e-3; 0; 0]; % [km / s]
x_f = [r_f; v_f];

% Time for relative orbit transfer
t_f = 3600; % [s]

t_span = [0, t_f];

% Problem Parameters
N = 50; % Number of time discretization points
Nu = N - 1; % Number of control steps to find

nx = 6; % Number of states
nu = 3; % Number of controls

t_k = linspace(t_span(1), t_span(2), N);
delta_t = t_k(2) - t_k(1);

% Algorithm parameters
default_tolerance = 1e-12;
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

% Chief orbit parameters
r_c = 6728; % [km] ISS Keplerian circular orbit
n = sqrt(mu / r_c ^ 3); % Mean motion of orbit

%% Linearization
f = @(t, x, u, p) CWH_EoM(t, x, u, mu, r_c);
A = @(t, x, u, p) [zeros(3), eye(3); [3 * n ^ 2, 0, 0; 0, 0, 0; 0, 0, -n ^ 2], [0, 2 * n, 0; -2 * n, 0, 0; 0, 0, 0]];
B = @(t, x, u, p) [zeros(nu); eye(nu)];
c = @(t, x, u, p) f(t, x, u, p) - A(t, x, u, p) * x - B(t, x, u, p) * u;

% Reference trajectory (not needed because the dynamics are linear)
x_ref = x_0 .* ones([nx, N]);
u_ref = zeros([nu, Nu]);

%% Discretization
% Discretize time by integrating the effect of the dynamics and control
% over each timestep to get transition matrices
[A_k, B_k, E_k, c_k, Delta] = discretize_dynamics_ZOH(f, A, B, @(t, x, u, p) zeros([nx, 0]), c, N, t_k, x_ref, u_ref, [], tolerances);

%% Solve with CVX
% Create optimization problem (CVX automatically solves it after cvx_end)
cvx_begin
    variable X_sol(nx, N) % State history
    variable U_sol(nu, Nu) % Control history
    minimize( sum(norms(U_sol, 2, 1) * delta_t) ) % Minimize total control use
    subject to
        norms(U_sol, 2, 1) <= u_max; % Maximum control acceleration constraint
        X_sol(:, 1) == x_0; % Initial condition
        X_sol(:, end) == x_f; % Terminal condition
        for k = 1:Nu % Loop for dynamics
            X_sol(:, k + 1) == A_k(:, :, k) * X_sol(:, k) + B_k(:, :, k) * U_sol(:, k) + c_k(:, :, k);
        end
cvx_end % End definition and solve problem

%% Plot Solution
figure
plot3(X_sol(1, :), X_sol(2, :), X_sol(3, :)); hold on
quiver3(X_sol(1, 1:Nu), X_sol(2, 1:Nu), X_sol(3, 1:Nu), U_sol(1, :), U_sol(2, :), U_sol(3, :)); hold on
scatter3(x_0(1), x_0(2), x_0(3), 48, "green", "filled", "square"); hold on
scatter3(x_f(1), x_f(2), x_f(3), 48, "red", "x"); hold off
title("Convex Minimum Fuel Control of Relative Orbit Transfer Near ISS")
xlabel("r [km]")
ylabel("\theta [km]")
zlabel("n [km]")
legend("Trajectory", "Thrust", "Start", "End", location = "eastoutside")
grid on

figure
tiledlayout(1, 3, "TileSpacing","compact")

nexttile
plot(t_k, X_sol(1, :), DisplayName="r"); hold on
plot(t_k, X_sol(2, :), DisplayName="\theta"); hold on
plot(t_k, X_sol(3, :), DisplayName="n"); hold off
title("RTN Frame Position History")
xlabel("Time [s]")
ylabel("Position [km]")
legend(Location="southoutside", Orientation="horizontal")
xlim(t_span)
grid on

nexttile
plot(t_k, X_sol(4, :), DisplayName="rdot"); hold on
plot(t_k, X_sol(5, :), DisplayName="\theta dot"); hold on
plot(t_k, X_sol(6, :), DisplayName="ndot"); hold off
title("RTN Frame Velocity History")
xlabel("Time [s]")
ylabel("Velocity [km / s]")
legend(Location="southoutside", Orientation="horizontal")
xlim(t_span)
grid on

nexttile

u_mag = vecnorm(U_sol, 2, 1);

stairs(t_k(1:Nu), U_sol(1, :) * 1e6, DisplayName="u_r"); hold on
stairs(t_k(1:Nu), U_sol(2, :) * 1e6, DisplayName="u_\theta"); hold on
stairs(t_k(1:Nu), U_sol(3, :) * 1e6, DisplayName="u_n"); hold on
stairs(t_k(1:Nu), u_mag * 1e6, DisplayName="||u||"); hold off
title("RTN Frame Control History")
xlabel("Time [s]")
ylabel("Control [mm / s2]")
legend(Location="southoutside", Orientation="horizontal")
xlim(t_span)
grid on

sgtitle("State and Control Time Histories for Relative Orbit Transfer Around ISS")

%% Helper Functions
function [xdot] = CWH_EoM(t, x, u, mu, r_c)
    r = x(1:3);
    v = x(4:6);

    n = sqrt(mu / r_c ^ 3);

    rdot = v;
    vdot = [3 * n ^ 2, 0, 0; 0, 0, 0; 0, 0, -n ^ 2]  * r ...
        +  [0, 2 * n, 0; -2 * n, 0, 0; 0, 0, 0] * v ...
        + u;

    xdot = [rdot; vdot];
end

function [xdot] = linearized_relative_orbit_EoM(t, x, u, mu, r_c, thetadot, thetaddot)
    r = x(1:3);
    v = x(4:6);

    rdot = v;
    vdot = [thetadot ^ 2 + 2 * mu / r_c ^ 3, thetaddot, 0; 
           -thetaddot, thetadot ^ 2 - mu / r_c ^ 3, 0; 
            0, 0, -mu / r_c ^ 3] * r ...
         + [0, 2 * thetadot, 0; -2 * thetadot, 0, 0; 0, 0, 0] * v ...
         + u;

    xdot = [rdot; vdot];
end

function [xdot] = nonlinear_relative_orbit_EoM(t, x, u, mu, r_c, r_cdot, nudot)
    r = x(1:3);
    v = x(4:6);

    r_d = norm([r_c; 0; 0] + r);

    rdot = v;
    vdot = [2 * nudot * (v(2) - r(2) * r_cdot / r_c) + r(1) * nudot ^ 2 + mu / r_c ^ 2 - mu / r_d ^ 3 * (r_c + r(1));
           -2 * nudot * (v(1) - r(1) * r_cdot / r_c) + r(2) * nudot ^ 2 - mu / r_d ^ 3 * r(2);
           -mu / r_d ^ 3 * r(3)] ...
         + u;

    xdot = [rdot; vdot];
end