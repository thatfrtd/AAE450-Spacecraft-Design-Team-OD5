%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% HW8 Q1c
% Author: Travis Hastreiter 
% Created On: 11 April, 2025
% Description: Lyapunov orbit transfer with minimum periapsis constraint
% Most Recent Change: 11 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = 10000; % [km] semi-major axis
e_c = 1e-3; % [] eccentricity
i_c = deg2rad(1e-3); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu_c = deg2rad(0); % [rad] true anomaly at epoch

M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);

% Initial conditions for spacecraft
a_d = 2 * (R_E + 500); % [km] semi-major axis
e_d = 0.3; % [] eccentricity
i_d = deg2rad(45); % [rad] inclination
Omega_d = deg2rad(0); % [rad] right ascension of ascending node
omega_d = deg2rad(0); % [rad] argument of periapsis
nu_d = deg2rad(0); % [rad] true anomaly at epoch

M_d = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_d, e_d), e_d);
x0_d_keplerian = [a_d; e_d; i_d; Omega_d; omega_d; M_d];
x0_d_cartesian = keplerian_to_cartesian(x0_d_keplerian, nu_d, mu_E);

[l_star, t_star] = nondimensionalized_quantities(a_c, mu_E);

% NEED TO ADD MASS IN STATE
Isp = 3000; % [s]
mass = 800; % [kg]
u_max = 0.001e-3; % [km / s2]

u_max_star = u_max / l_star * t_star ^ 2;
g_0 = 9.81; % [m / s2]

% Propagation Time 
orbits = 1200;
tspan = linspace(0, orbits * period(a_c, mu_E), 1e4);
t_orbits = linspace(0, orbits, numel(tspan));
t_hr = tspan / 60 / 60;

% Integration error tolerance
default_tolerance = 1e-6;

%% Penalized Lyapunov - periapsis constraint
a_d_0 = @(t, x) zeros([3, 1]);

tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

% Integrate target orbit
[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_keplerian(x, 1), B_keplerian(x, 1), a_d_0(t,x)), tspan / t_star, x0_c_keplerian .* [1 / l_star, ones([1, 5])]', tolerances);
x_keplerian_cartesian_c = keplerian_to_cartesian_array(x_keplerian_c ./ [1 / l_star, ones([1, 5])], [], mu_E);

% Define Lyapunov feedback controller
K = eye(5) * sqrt(2); % Gain matrix
type = "quadratic"; % Penalty type
k = 100; % Smoothing parameter
w = 100; % Penalty weight
r_p_min = R_E + 400; % [km] min periapsis
g = @(x_slow) r_p_min / l_star - x_slow(1) * (1 - x_slow(2)); % Min periapsis constraint function r_p_d < r_p_c
g_partial_xslow = @(x_slow) [x_slow(2) - 1, x_slow(1), 0, 0, 0]; % Constraint jacobian
epsilon = 1e-8;
% Control acceleration function
a_d_pkl = @(t, x) penalized_keplerian_lyapunov_bounded(x, interp1(t_keplerian_c, x_keplerian_c, t, "linear","extrap")', K, 1, w, g(x(1:5)), g_partial_xslow(x(1:5)), k, epsilon, type, u_max_star);

% Integrate spacecraft orbit with Lyapunov feedback control
[t_keplerian_d, x_keplerian_d] = ode45(@(t,x) gauss_planetary_eqn(f0_keplerian(x, 1), B_keplerian(x, 1), a_d_pkl(t,x)), tspan / t_star, x0_d_keplerian .* [1 / l_star, ones([1, 5])]', tolerances);
x_keplerian_cartesian_d = keplerian_to_cartesian_array(x_keplerian_d ./ [1 / l_star, ones([1, 5])], [], mu_E);

% Convert control acceleration to thrust history
u_cl_crit_penalty = zeros([numel(tspan), 3]);
for i = 1:numel(tspan)
    u_cl_crit_penalty(i, :) = a_d_pkl(t_keplerian_d(i), x_keplerian_d(i, :)') * l_star / t_star ^ 2;
end

%% Plot
figure
tiledlayout(1,3,"TileSpacing","compact")
nexttile
earthy(R_E, "Earth", 0.5, [0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d, 'b', 0, 1); hold off

title("Penalized Keplerian Lyapunov Orbit Transfer")
legend("", "Target", "", "Spacecraft")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")


plot_orbit_transfer_histories(t_hr, x_keplerian_c, x_keplerian_d, u_cl_crit_penalty)
sgtitle("Keplerian Lyapunov Orbit Transfer with Eccentricity Constraint Results")

%%
plot(t_hr, x_keplerian_d(:, 1) .* (1 - x_keplerian_d(:, 2)) * l_star - R_E); hold on
yline(r_p_min - R_E); hold off
xlabel("Time [hr]")
ylabel("Periapsis [km]")
legend("Spacecraft", "Target Orbit")
title("Lyapunonv Orbit Transfer Eccentricity Constraint")
grid on
%%
figure
tiledlayout(1,1,"TileSpacing","compact")
nexttile
earthy(R_E, "Earth", 0.5, [0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d, 'b', 0, 1); hold off

title("Penalized Keplerian Lyapunov Orbit Transfer")
legend("", "Target", "", "Spacecraft")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")

figure
plot_orbit_transfer_histories(t_hr, x_keplerian_c, x_keplerian_d, u_cl_crit_penalty)

%% Calculate final mass
dV_spiral = sqrt(mu_E / a_d) - sqrt(mu_E / a_c)
dV_total = sum(vecnorm(u_cl_crit_penalty .* (tspan(2) - tspan(1)), 2, 2))
m_f = mass * exp(-dV_total * 1000 / (Isp * g_0))


%% Helper Functions

function [f_0] = f0_keplerian(x, mu)
    a = x(1);

    n = sqrt(mu / a ^ 3);

    f_0 = [zeros([5, 1]); ...
           n];
end

function [B] = B_keplerian(x, mu)
    a = x(1);
    e = x(2);
    i = x(3);
    Omega = x(4);
    omega = x(5);
    M = x(6);

    p = a * (1 - e ^ 2);
    b = a * sqrt(1 - e ^ 2);
    E = mean_to_eccentric_anomaly(M, e);
    nu = eccentric_to_true_anomaly(E, e);
    r = a * (1 - e * cos(E));
    h = sqrt(mu * p);

    B = 1 / h * [2 * a ^ 2 * e * sin(nu), 2 * a ^ 2 * p / r, 0;
         p * sin(nu), (p + r) * cos(nu) + r * e, 0;
         0, 0, r * cos(nu + omega);
         0, 0, r * sin(nu + omega) / sin(i);
         -p * cos(nu) / e, (p + r) * sin(nu) / e, -r * sin(nu + omega) / tan(i);
         b * p * cos(nu) / (a * e) - 2 * b * r / a, -b * (p + r) * sin(nu) / (a * e), 0];
    
    % Make B convert disturbance into RTN frame from cartesian frame
    B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)';
end

function [x_dot] = gauss_planetary_eqn(f_0, B, a_d)
    x_dot = f_0 + B * a_d;
end


function [a_kep_lyapunov_bounded] = penalized_keplerian_lyapunov_bounded(x, x_c, K, mu, w, g, g_partial_xslow, k, epsilon, type, u_max)
    delta_x_slow = x(1:5) - x_c(1:5);

    B = B_keplerian(x, mu);
    B_slow = B(1:5, :);

    [P, P_partial_g] = penalty_function(g, k, epsilon, type);
    L = (2 * (1 + w * P) * K + w * K * delta_x_slow * P_partial_g * g_partial_xslow) * B_slow;

    a_kep_lyapunov = -(L' * L) ^ (-1) * L' * delta_x_slow;

    if norm(a_kep_lyapunov) > u_max
        a_kep_lyapunov_bounded = a_kep_lyapunov / norm(a_kep_lyapunov) * u_max;
    else
        a_kep_lyapunov_bounded = a_kep_lyapunov;
    end
end

function [P, P_partial_g] = penalty_function(g, k, epsilon, type)

    if type == "linear"
        P = k * max(0, g + epsilon);
        P_partial_g = k * max(0, g + epsilon) / (g + epsilon);
    elseif type == "quadratic"
        P = k / 2 * max(0, g + epsilon) ^ 2;
        P_partial_g = k * max(0, g + epsilon);
    elseif type == "exponential"
        P = exp(k * (g + epsilon)); 
        P_partial_g = P * k;
    end
end

function [l_star, t_star] = nondimensionalized_quantities(a_c, mu)
    l_star = a_c; % [km]

    t_star = sqrt(l_star ^ 3 / mu); % [s]
end

function [] = plot_orbit_transfer_histories(t_hr, x_c, x_d, u)
    delta_x = x_c(:, 1:5) - x_d(:, 1:5);
        
    nexttile
    plot(t_hr, delta_x(:, 1), t_hr, delta_x(:, 2), t_hr, delta_x(:, 3), t_hr, delta_x(:, 4), t_hr, delta_x(:, 5)); hold off 
    title("Nondimensionalized Orbital Element Error")
    xlabel("Time [hr]")
    ylabel("\delta x")
    legend("a", "e", "i", "\omega", "\Omega", Location="northeast")
    grid on
    
    nexttile
    plot(t_hr, u(:, 1), t_hr, u(:, 2), t_hr, u(:, 3)); hold on 
    plot(t_hr, vecnorm(u, 2, 2), LineStyle= "--"); hold off
    title("Control")
    xlabel("Time [hr]")
    ylabel("u [km / s2]")
    legend("u_1", "u_2", "u_3", "||u||", Location="northeast")
    grid on
end
