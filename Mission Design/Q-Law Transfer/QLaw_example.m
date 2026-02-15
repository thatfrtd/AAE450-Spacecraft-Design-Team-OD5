%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Q-Law Orbit Transfer Example 
% Author: Travis Hastreiter 
% Created On: 8 February, 2026
% Description: Orbit transfer using Q-Law with minimum periapsis constraint
% Most Recent Change: 8 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: - Add eclipse detection and forbide thrusting during them
%       - Add J2 perturbations
%       - Add option for drag for low orbits

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
e_d = 1e-4; % [] eccentricity
i_d = deg2rad(10); % [rad] inclination
Omega_d = deg2rad(0); % [rad] right ascension of ascending node
omega_d = deg2rad(0); % [rad] argument of periapsis
nu_d = deg2rad(0); % [rad] true anomaly at epoch

M_d = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_d, e_d), e_d);
x0_d_keplerian = [a_d; e_d; i_d; Omega_d; omega_d; M_d];
x0_d_cartesian = keplerian_to_cartesian(x0_d_keplerian, nu_d, mu_E);

char_star = load_charecteristic_values_Earth();

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 3000; % [s]
spacecraft_params.m_0 = 800; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 1; % [N]

% Integration error tolerance
default_tolerance = 1e-6;

%% Q Law Transfer
a_d_0 = @(t, x) zeros([3, 1]); % Disturbance function

% Min Periapsis soft constraint
penalty_params = struct();
penalty_params.k = 100; % Smoothing parameter
penalty_params.W_p = 1; % Penalty weight
penalty_params.r_p_min = R_E + 400; % [km] min periapsis

% Define Q-Law feedback controller: W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
Q_params = struct();
Q_params.W_oe = 1 * ones([5, 1]); % Element weights 
Q_params.eta_a_min = 0.1; % Minimum absolute efficiency for thrusting instead of coasting
Q_params.eta_r_min = 0.1; % Minimum relative efficiency for thrusting instead of coasting
Q_params.m = 3;
Q_params.n = 4;
Q_params.r = 2;
Q_params.Theta_rot = 0;

% Parameters for the optimization needed to determine efficiencies
Qdot_opt_params = struct();
Qdot_opt_params.num_start_points = 10;
Qdot_opt_params.strategy = "Best Start Points";
Qdot_opt_params.plot_minQdot_vs_L = false;

[Qtransfer] = QLaw_transfer(x0_d_keplerian, x0_c_keplerian, mu_E, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 50000, angular_step=deg2rad(1));

if Qtransfer.converged
    fprintf("Q-Law Transfer Converged! Took %.3f Days Using %.3f kg Propellant\n", Qtransfer.dt / 60 / 60 / 24, Qtransfer.delta_m)
else
    fprintf("Q-Law Transfer Failed with %s\n", Qtransfer.errors)
end

%% Integrate Target Orbit
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);
[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_keplerian(x, 1), B_keplerian(x, 1), a_d_0(t,x)), Qtransfer.t / char_star.t, x0_c_keplerian .* [1 / char_star.l, ones([1, 5])]', tolerances);
x_keplerian_c = x_keplerian_c';
x_keplerian_c(1, :) = x_keplerian_c(1, :) .* char_star.l; % Redimensionalize

%% Plot Orbit
x_keplerian_cartesian_d = keplerian_to_cartesian_array(Qtransfer.x_keplerian_mass(1:6, :), [], mu_E);

not_coast_colors = interp1(1:numel(Qtransfer.not_coast), double(Qtransfer.not_coast), linspace(1, numel(Qtransfer.not_coast), numel(Qtransfer.t)), 'nearest');

figure
plot_cartesian_orbit_color_varying(x_keplerian_cartesian_d(:, 1:1:end), not_coast_colors, 3); hold on
c = colorbar;
plotOrbit3(Omega_c, i_c, omega_c, a_c * (1 - e_c ^2), e_c, linspace(0, 2 * pi, 1000), "r", 1, 1, [0, 0, 0], 0.1, 2); hold on
plotOrbit3(Omega_d, i_d, omega_d, a_d * (1 - e_d ^2), e_d, linspace(0, 2 * pi, 1000), "g", 1, 1, [0, 0, 0], 0.1, 2)
grid on
earthy(R_E, "Earth", 0.5, [0;0;0]); hold on; 
clim([0, 1])
c.Ticks = [0, 0.5, 1];
c.TickLabels = {'Coast', 'Eclipse', 'Thrust'};
axis equal
title("Q-Law Orbit Transfer")
legend("Spacecraft", "Target", "Initial")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")

%% Plot Orbit Error and Control Histories
figure
plot_orbit_transfer_histories(Qtransfer.t / 60, x_keplerian_c', Qtransfer.x_keplerian_mass(1:6, :)', interp1(1:numel(Qtransfer.not_coast), Qtransfer.u', linspace(1, numel(Qtransfer.not_coast), numel(Qtransfer.t)), 'nearest'));
sgtitle("Q-Law Orbit Transfer with Periapsis Constraint Results")

%% Plot Q Function
figure
plot(interp1(1:numel(Qtransfer.t), Qtransfer.t, linspace(1, numel(Qtransfer.t), numel(Qtransfer.Q))) / 60 / 60 / 24, Qtransfer.Q)
xlabel("Time [days]")
ylabel("Q-Function []")
title("Q-Function vs Time")
yscale("log") % what's the best for plotting this?
grid on


%% Plot Constraint Satisfaction
% figure
% tiledlayout(1, 2)
% 
% nexttile
% plot(Qtransfer.t / 60, x_keplerian_d(:, 1) .* (1 - x_keplerian_d(:, 2)) * l_star - R_E); hold on
% yline(r_p_min - R_E); hold off
% xlabel("Time [hr]")
% ylabel("Periapsis [km]")
% title("Q-Law Orbit Transfer Periapsis Constraint")
% grid on
% 
% nexttile
% plot(Qtransfer.t / 60, Qtransfer.P)
% xlabel("Time [hr]")
% ylabel("Periapsis Constraint []")
% title("Q-Law Orbit Transfer Periapsis Constraint Nondimensionalized")

%% Calculate Spiral Transfer Estimates
dV_spiral = sqrt(mu_E / a_d) - sqrt(mu_E / a_c)
m_f_spiral = spacecraft_params.m_0 * (1 - exp(-abs(dV_spiral) * 1000 / (spacecraft_params.Isp * 9.81)))

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
