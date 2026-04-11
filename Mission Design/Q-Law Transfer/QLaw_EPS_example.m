%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Q-Law Orbit Transfer with Integrated Electrical Power System Example 
% Author: Travis Hastreiter 
% Created On: 8 February, 2026
% Description: Orbit transfer using Q-Law with minimum periapsis constraint
% and state of charge constraint
% Most Recent Change: 14 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: - Add J2 perturbations
%       - Add option for drag for low orbits

R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter

thrust_during_eclipse = true;

% Electrical Power System parameters: .battery_capacity, .panel_area, .panel_efficiency, .min_state_of_charge, .thruster_power
EPS_params = struct();
EPS_params.battery_capacity = 7200; % [Wh]
EPS_params.panel_area = 40; % [m2]
I_d = 0.7; % Degredation factor
EPS_params.panel_efficiency = 0.259 * I_d; % [] max 40% - worth going on lower end (starts at 0.295)
EPS_params.min_state_of_charge = 0.25; % min state of charge [fraction of capacity]
EPS_params.thruster_power = 8000; % [W] thruster power use
EPS_params.nominal_power = 700; % [W] draw of spacecraft during operation
% Check panel 
solar_irradiance = 1361; % [W / m2]
power_gen = solar_irradiance * EPS_params.panel_efficiency * EPS_params.panel_area;

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
r_a_c = R_E + 900; % [km] apoapsis
r_p_c = R_E + 140; % [km] periapsis
e_c = (1 - r_p_c / r_a_c) / (1 + r_p_c / r_a_c); % [] eccentricity
a_c = r_p_c / (1 - e_c); % [km] semi-major axis
i_c = deg2rad(98); % [rad] inclination
Omega_c = deg2rad(10); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu_c = deg2rad(0); % [rad] true anomaly at epoch

% Initial conditions for spacecraft
r_a_d = R_E + 901; % [km] apoapsis
r_p_d = R_E + 900; % [km] periapsis
e_d = (1 - r_p_d / r_a_d) / (1 + r_p_d / r_a_d); % [] eccentricity
a_d = r_p_d / (1 - e_d); % [km] semi-major axis
i_d = deg2rad(98); % [rad] inclination
Omega_d = deg2rad(10); % [rad] right ascension of ascending node
omega_d = deg2rad(0); % [rad] argument of periapsis
nu_d = deg2rad(0); % [rad] true anomaly at epoch


M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);


M_d = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_d, e_d), e_d);
x0_d_keplerian = [a_d; e_d; i_d; Omega_d; omega_d; M_d];
x0_d_cartesian = keplerian_to_cartesian(x0_d_keplerian, nu_d, mu_E);

char_star = load_charecteristic_values_Earth();

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 4100; % [s]
spacecraft_params.m_0 = 6000; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.235; % [N]

% Integration error tolerance
default_tolerance = 1e-12;

%% Q Law Transfer
a_d_0 = @(t, x) zeros([3, 1]); % Disturbance function

% Min Periapsis soft constraint
penalty_params = struct();
penalty_params.k = 100; % Smoothing parameter
penalty_params.W_p = 1; % Penalty weight
penalty_params.r_p_min = R_E + 0; % [km] min periapsis

% Parameters for the optimization needed to determine efficiencies
Qdot_opt_params = struct();
Qdot_opt_params.num_start_points = 10;
Qdot_opt_params.strategy = "Best Start Points";
Qdot_opt_params.plot_minQdot_vs_L = false;

N_i = 1;
eta = linspace(0.3, 0.3, N_i);
clear Qtransfer
for i = 1 : N_i
    % Define Q-Law feedback controller: W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
    Q_params = struct();
    Q_params.W_oe = 1 * ones([5, 1]); % Element weights 
    Q_params.eta_a_min = eta(i); % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = eta(i); % Minimum relative efficiency for thrusting instead of coasting
    Q_params.m = 3;
    Q_params.n = 4;
    Q_params.r = 2;
    Q_params.Theta_rot = 0;

    tic;
    [Qtransfer(i)] = QLaw_transfer(x0_d_keplerian, x0_c_keplerian, mu_E, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, EPS_params, return_dt_dm_only = false, iter_max = 1500000, angular_step=deg2rad(5), thrust_during_eclipse = thrust_during_eclipse, integration_tolerance=1e-10, max_t=356.25*1.5*60*60*24);
    opt_time = toc
end

dVs = zeros([N_i, 1]);
ToFs = zeros([N_i, 1]);
for i = 1 : N_i
    dVs(i) = Qtransfer(i).delta_V;
    ToFs(i) = Qtransfer(i).dt / 60 / 60 / 24;
    if Qtransfer(i).converged
        fprintf("Q-Law Transfer Converged! Took %.3f Days Using %.3f kg Propellant\n", Qtransfer(i).dt / 60 / 60 / 24, Qtransfer(i).delta_m)
    else
        fprintf("Q-Law Transfer Failed with %s\n", Qtransfer(i).errors)
    end
end

Qtransfer = Qtransfer(round(N_i / 2));

%%
figure
plot(ToFs, dVs)
title("Delta V vs ToF Pareto")
xlabel("Tof [days]")
ylabel("Delta V [km / s]")
grid on

%% Integrate Target Orbit
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);
[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_keplerian(x, 1), B_keplerian(x, 1), a_d_0(t,x)), Qtransfer.t / char_star.t, x0_c_keplerian .* [1 / char_star.l, ones([1, 5])]', tolerances);
x_keplerian_c = x_keplerian_c';
x_keplerian_c(1, :) = x_keplerian_c(1, :) .* char_star.l; % Redimensionalize

%% Plot Orbit
x_keplerian_cartesian_d = keplerian_to_cartesian_array(Qtransfer.x_keplerian_mass(1:6, :), [], mu_E);

if ~thrust_during_eclipse
    not_coast_colors = interp1(1:numel(Qtransfer.not_coast), double(Qtransfer.not_coast) + 0.5 * double(Qtransfer.eclipsed), linspace(1, numel(Qtransfer.not_coast), numel(Qtransfer.t)), 'nearest');
else
    not_coast_colors = interp1(1:numel(Qtransfer.not_coast), double(Qtransfer.not_coast), linspace(1, numel(Qtransfer.not_coast), numel(Qtransfer.t)), 'nearest');
end

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

%%
figure
plot(Qtransfer.t / 60 / 60 / 24, vecnorm(x_keplerian_cartesian_d(1:3, :)) - R_E)
grid on

%%
figure
plot(Qtransfer.t_disc / 60 / 60 / 24, Qtransfer.battery_SoC); hold on
stairs(Qtransfer.t_disc(1:end-1) / 60 / 60 / 24, Qtransfer.eclipsed, "--")
stairs(Qtransfer.t_disc(1:end-1) / 60 / 60 / 24, Qtransfer.not_coast, "--")
xlabel("Time [days]")
ylabel("Value")
legend("SoC", "Eclipsed", "Thrusting")
grid on

%% Plot Orbit Error and Control Histories
figure
plot_orbit_transfer_histories(Qtransfer.t / 60 / 60, x_keplerian_c' ./ [R_E, ones([1, 5])], Qtransfer.x_keplerian_mass(1:6, :)' ./ [R_E, ones([1, 5])], Qtransfer.u_cont');
sgtitle("Q-Law Orbit Transfer with Periapsis Constraint Results")

%% Plot Q Function
figure
plot(interp1(1:numel(Qtransfer.t), Qtransfer.t, linspace(1, numel(Qtransfer.t), numel(Qtransfer.Q))) / 60 / 60 / 24, Qtransfer.Q)
xlabel("Time [days]")
ylabel("Q-Function []")
title("Q-Function vs Time")
yscale("log") % what's the best for plotting this?
grid on


%% Validate Solution
% tight_tolerance = 1e-9;
% tight_tolerances = odeset(RelTol=tight_tolerance, AbsTol=tight_tolerance);
% g_0 = 9.81e-3; % [km / s2]
% a_disturbance = @(t,x) [0;0;0];
% x0_kep_mass_d = [x0_d_keplerian; spacecraft_params.m_0];
% u = @(t) interp1(Qtransfer.t, Qtransfer.u_cont', t)';
% %u = @(t) Qtransfer.u_cont(:, min(floor(t ./ Qtransfer.t(end) * numel(Qtransfer.t)) + 1, numel(Qtransfer.t)));
% a_control = @(t, x) u(t) / x(end);
% mdot = @(t) -norm(u(t)) / (spacecraft_params.Isp * g_0);
% [~, x_kep_mass_d_ck] = ode45(@(t, x) [gauss_planetary_eqn(f0_keplerian(x, mu_E), B_keplerian(x, mu_E), a_control(t, x) + a_disturbance(t, x)); mdot(t)], Qtransfer.t, x0_kep_mass_d, tight_tolerances);
% x_kep_mass_d_ck = x_kep_mass_d_ck';
% 
% %%
% x_keplerian_cartesian_ck = [keplerian_to_cartesian_array(x_kep_mass_d_ck(1:6, :), [], mu_E); x_kep_mass_d_ck(7, :)];
% 
% %%
% 
% figure
% plot_cartesian_orbit(x_keplerian_cartesian_d(:, 1:1:end)); hold on
% plot_cartesian_orbit(x_keplerian_cartesian_ck(:, 1:1:end)); hold on
% c = colorbar;
% plotOrbit3(Omega_c, i_c, omega_c, a_c * (1 - e_c ^2), e_c, linspace(0, 2 * pi, 1000), "r", 1, 1, [0, 0, 0], 0.1, 2); hold on
% plotOrbit3(Omega_d, i_d, omega_d, a_d * (1 - e_d ^2), e_d, linspace(0, 2 * pi, 1000), "g", 1, 1, [0, 0, 0], 0.1, 2)
% grid on
% earthy(R_E, "Earth", 0.5, [0;0;0]); hold on; 
% clim([0, 1])
% c.Ticks = [0, 0.5, 1];
% c.TickLabels = {'Coast', 'Eclipse', 'Thrust'};
% axis equal
% 
% title("Q-Law Orbit Transfer")
% legend("Nondim", "Check", "Target", "Initial")
% xlabel("X [km]")
% ylabel("Y [km]")
% zlabel("Z [km]")


% %% Calculate Spiral Transfer Estimates
% dV_spiral = sqrt(mu_E / a_d) - sqrt(mu_E / a_c)
% m_f_spiral = spacecraft_params.m_0 * (1 - exp(-abs(dV_spiral) * 1000 / (spacecraft_params.Isp * 9.81)))

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
    % B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)';
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
