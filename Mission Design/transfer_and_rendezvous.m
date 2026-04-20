%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Combined Transfer from Insertion Orbit and Target Rendezvous
% Author: Travis Hastreiter 
% Created On: 25 February, 2026
% Description: Orbit transfer using Q-Law and rendezvous using SCP in Hill 
% frame.
% Most Recent Change: 25 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_debris = 3000; % [kg]

R_E = 6378.137; % [km] Earth radius
mu_E = 398600.4418; % [km3 / s2] Earth gravitational parameter
J_2_val = 1.08262668e-3*0; % [] Earth J2

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = R_E + 800; % [km] semi-major axis
e_c = 0.003390; % [] eccentricity
i_c = deg2rad(98.3 ); % [rad] inclination
Omega_c = deg2rad(320.5520 ); % [rad] right ascension of ascending node
omega_c = deg2rad(301.2069 ); % [rad] argument of periapsis
nu_c = deg2rad(58.6658 ); % [rad] true anomaly at epoch
n = sqrt(mu_E / a_c ^ 3);

M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);

% Initial conditions for spacecraft
a_d = 7139; % [km] semi-major axis
e_d = 0.00001; % [] eccentricity
i_d = deg2rad(98.2); % [rad] inclination
Omega_d = deg2rad(320.5520 ); % [rad] right ascension of ascending node
omega_d = deg2rad(301.2069); % [rad] argument of periapsis
nu_d = deg2rad(58.6658 ); % [rad] true anomaly at epoch

M_d = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_d, e_d), e_d);
x0_d_keplerian = [a_d; e_d; i_d; Omega_d; omega_d; M_d];
x0_d_cartesian = keplerian_to_cartesian(x0_d_keplerian, nu_d, mu_E);

char_star = load_charecteristic_values_Earth();

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 4100; % [s]
spacecraft_params.m_0 = 1500; % [kg]
spacecraft_params.m_dry = 1000; % [kg]
spacecraft_params.F_max = 0.235; % [N]

% Integration error tolerance
default_tolerance = 1e-10;

%% Q-Law Transfer to Object (First find time then adjust to phase properly)
a_d_0 = @(t, x) zeros([3, 1]); % Disturbance function
a_d_J2 = @(t,x) J2_perturbation(modified_equinoctial_to_cartesian(x(1:6), mu_E), mu_E, R_E, J_2_val);

% Min Periapsis soft constraint
penalty_params = struct();
penalty_params.k = 100; % Smoothing parameter
penalty_params.W_p = 1; % Penalty weight
penalty_params.r_p_min = R_E + 90; % [km] min periapsis

% Parameters for the optimization needed to determine efficiencies
Qdot_opt_params = struct();
Qdot_opt_params.num_start_points = 10;
Qdot_opt_params.strategy = "Best Start Points";
Qdot_opt_params.plot_minQdot_vs_L = false;

eta = 0;
clear Qtransfer
% Define Q-Law feedback controller: W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
Q_params = struct();
Q_params.W_oe = 1 * ones([5, 1]); % Element weights 
Q_params.eta_a_min = eta; % Minimum absolute efficiency for thrusting instead of coasting
Q_params.eta_r_min = eta; % Minimum relative efficiency for thrusting instead of coasting
Q_params.m = 3;
Q_params.n = 4;
Q_params.r = 2;
Q_params.Theta_rot = 0;

[Qtransfer] = QLaw_transfer_fast(x0_d_keplerian, x0_c_keplerian, mu_E, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 500000, angular_step=deg2rad(2), R_c = 0.01, a_disturbance = a_d_0, integration_tolerance=1e-10);

dVs = Qtransfer.delta_V;
ToFs = Qtransfer.dt / 60 / 60 / 24;
if Qtransfer.converged
    fprintf("Q-Law Transfer Converged! Took %.3f Days Using %.3f kg Propellant\n", Qtransfer.dt / 60 / 60 / 24, Qtransfer.delta_m)
else
    fprintf("Q-Law Transfer Failed with %s\n", Qtransfer.errors)
end

plot_Q_transfer(Qtransfer, x0_c_keplerian, x0_d_keplerian, char_star)

%% Make sure transfer correct
a_d_J2_ECI = @(t,x) J2_perturbation_ECI(x, mu_E, R_E, J_2_val);
a_d_RTNcontrol = @(t,x) a_d_J2_ECI(t,x) + RTN_to_ECI(x(1:3), x(4:6)) * interp1(Qtransfer.t, Qtransfer.u_cont' ./ Qtransfer.x_keplerian_mass(7, :)', t, "previous")';

% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

% Run first to see how far off orbit gets from desired because of J2
[~,x_keplerian_cartesian_preck] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_RTNcontrol(t,x)), Qtransfer.t, x0_d_cartesian, tolerances);
x_keplerian_preck = cartesian_to_keplerian_array(x_keplerian_cartesian_preck', [0; 0; 1], [1; 0; 0], mu_E);

% Run with adjusted initial orbit attempting to account for J2 drifting
% RAAN and arg periapsis
x0_d_keplerian_adj = [a_d; e_d; i_d; Omega_d - wrapToPi(x_keplerian_preck(4, end) - x0_c_keplerian(4)); omega_d - wrapToPi(x_keplerian_preck(5, end) - x0_c_keplerian(5)); M_d];
x0_d_cartesian_adj = keplerian_to_cartesian(x0_d_keplerian_adj, nu_d, mu_E);
[~,x_keplerian_cartesian_ck] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d_RTNcontrol(t,x)), Qtransfer.t, x0_d_cartesian_adj, tolerances);
x_keplerian_ck = cartesian_to_keplerian_array(x_keplerian_cartesian_ck', [0; 0; 1], [1; 0; 0], mu_E);

%%

x_Qtrans_cart = keplerian_to_cartesian_array(Qtransfer.x_keplerian_mass(1:6, :), [], mu_E); % QLaw output (not corrected for J2!)
figure
earthy(R_E, "Earth", 1, [0;0;0]); hold on;
plot3(x_keplerian_cartesian_preck(:, 1), x_keplerian_cartesian_preck(:, 2), x_keplerian_cartesian_preck(:, 3))
plot3(x_Qtrans_cart(1,:), x_Qtrans_cart(2,:), x_Qtrans_cart(3,:))

legend("Earth", "Check", "Sol")
axis equal
grid on

%% Find Transition from Transfer to Rendezvous
% Find where QLaw gets 1 km away and use those as initial conditions for
% SCP rendezvous
max_engagement_dist = 6; % [km] Distance to switch to rendezvous from transfer

x_d_cartesian = x_keplerian_cartesian_ck';%keplerian_to_cartesian_array(Qtransfer.x_keplerian_mass(1:6, :), [], mu_E);

delta_M_transfer = Qtransfer.x_keplerian_mass(6, end) - M_d;
delta_M_chief = sqrt(char_star.mu / a_c ^ 3) * Qtransfer.t(end);
M0_c_refined = wrapToPi(delta_M_transfer - delta_M_chief); % Assuming M_c is 0 DOES NOT WORK - JUST USE OPTIMIZATION FUNCTION

x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M0_c_refined + deg2rad(55.9)];
%x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M0_c_refined + deg2rad(77.845)];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, [], mu_E);
x_c_cartesian = propagate_conic_array(x0_c_cartesian, Qtransfer.t, char_star.mu);
%[~, x_c_cartesian] = ode45(@(t, x) gauss_planetary_eqn(f0_cartesian(x, char_star.mu), B_cartesian(x, char_star.mu), a_d_J2_ECI(t, x)), -flip(Qtransfer.t), x0_c_cartesian, tolerances); % Integrate back in time including J2...
%x_c_cartesian = flip(x_c_cartesian, 2);

relative_distances = vecnorm(x_c_cartesian(1:3, :) - x_d_cartesian(1:3, :));
i_engage = find(relative_distances < max_engagement_dist, 1);

if isempty(i_engage)
    % figure
    % plot(Qtransfer.t / 3600, relative_distances, HandleVisibility="off"); hold on
    % yscale("log")
    % xlabel("Time [hr]")
    % ylabel("Relative Distance [km]")
    % title("Distance from Target vs Time for Transfer")
    % grid on

    error("Spacecraft minimum distance to target, %.1f km, is not within engagement distance of %.1f km", min(relative_distances), max_engagement_dist)
end

x_c_engage_ECI = x_c_cartesian(:, i_engage);
x_d_engage_ECI = [x_d_cartesian(:, i_engage); Qtransfer.x_keplerian_mass(7, i_engage)];
x_d_engage_hill = ECI_to_Hill(x_c_engage_ECI, x_d_engage_ECI);
x_d_engage_ECI_hill = Hill_to_ECI(x_c_engage_ECI, x_d_engage_hill);

figure
plot(Qtransfer.t / 3600, relative_distances, HandleVisibility="off"); hold on
xline(Qtransfer.t(i_engage) / 3600, Label="Engagement Start")
yscale("log")
xlabel("Time [hr]")
ylabel("Relative Distance [km]")
title("Distance from Target vs Time for Transfer")
grid on

%%
x_Qtransfer_cart = [x_d_cartesian(:, 1:i_engage); Qtransfer.x_keplerian_mass(7, 1:i_engage)];
u_trans_cont = Qtransfer.u_cont(:, 1:i_engage);

t_trans = Qtransfer.t(1:i_engage);
x_trans = (x_Qtransfer_cart )';
u_accel_trans = reshape(pagemtimes(RTN_to_ECI_array(x_Qtransfer_cart(1:3, :), x_Qtransfer_cart(4:6, :)), reshape(u_trans_cont(:, 1:i_engage), 3, 1, [])), 3, [])' ./ x_trans(:, 7);

trans_table = table(t_trans, x_trans, u_accel_trans);
writetable(trans_table, "Wes_QLaw_TTC_test.csv")

%% Relative Orbit Transfer using SCP 
x_0_hill = x_d_engage_hill;
b = sqrt(40^2/2);
c = sqrt(10^2/2);
aROE = [0; % [km] delta semimajor axis
        0; % [km] delta lambda
        b*1e-3; % [km] delta e_x
        b*1e-3; % [km] delta e_y
        c*1e-3; % [km] delta i_x 
        c*1e-3]; % [km] delta i_y

cart_ROE = ROE_to_cart_matrix(n, 0) * aROE;
x_f_hill = cart_ROE; % [km / s]
ToF_rendezvous = 3600 * 5; % [s]
spacecraft_params_mod = spacecraft_params;
spacecraft_params_mod.Isp = [4155; 300];
spacecraft_params_mod.F_max = [0.235; 88];
% Actually use QLaw control for "nocont" part to check better
%[optimal_rendezvous, rendezvous_SCP_info, nocont_rendezvous] = nonlinear_rendezvous_func(x_0_hill, x_f_hill, ToF_rendezvous, x0_c_keplerian, spacecraft_params_mod, N = 150, dynamics = "Nonlinear", u_hold = "ZOH", max_iters = 15, integration_tolerance = 5e-12);
[optimal_rendezvous, rendezvous_SCP_info] = nonlinear_rendezvous_multithruster_func(x_0_hill, x_f_hill, ToF_rendezvous, x0_c_keplerian, spacecraft_params_mod, N = 250, u_hold = "ZOH", max_iters = 15, integration_tolerance = 5e-12);

%% Save Rendezvous
output_array = [optimal_rendezvous.t(2:end), optimal_rendezvous.x(1:3, 2:end)', optimal_rendezvous.u'];
state_names = ["r_1", "r_2", "r_3"];
control_names = ["u_main_1", "u_main_2", "u_main_3", "u_RCS_1", "u_RCS_2", "u_RCS_3"];
output_names = ["Time", state_names, control_names];

output_table = array2table(output_array, VariableNames = output_names);
writetable(output_table,"./Animation/close_rendezvous.csv")


%% Package Output
pickle = 800;
x_Qtransfer_cart = [x_d_cartesian(:, 1:i_engage); Qtransfer.x_keplerian_mass(7, 1:i_engage)];
x_Qtransfer_cart_plus = [x_d_cartesian(:, 1:i_engage + pickle); Qtransfer.x_keplerian_mass(7, 1:i_engage + pickle)];
x_Qtransfer_hill = ECI_to_Hill(x_c_cartesian(:, 1:i_engage), x_Qtransfer_cart);
x_Qtransfer_hill_plus = ECI_to_Hill(x_c_cartesian(:, 1:i_engage + pickle), x_Qtransfer_cart_plus);
x_Qtransfer_hill_ck = Hill_to_ECI(x_c_cartesian(:, 1:i_engage), x_Qtransfer_hill); % DOES NOT MATCH
t_cont_traj = [Qtransfer.t; optimal_rendezvous.t];
x_cont_traj = [x_Qtransfer_hill, optimal_rendezvous.x];

trajectory_to_target.transfer = Qtransfer;
trajectory_to_target.rendezvous = rendezvous_SCP_info;

%% CHECK ECI<->HILL CONVERSION -- DOESN'T MATCH :(
% ckdiff = x_Qtransfer_hill_ck - x_Qtransfer_cart;
% figure
% plot(x_Qtransfer_hill_ck(4, :)'); hold on
% plot(x_Qtransfer_cart(4, :)');

%%


%%
figure
plt_start_i = 23000; 
i_end = size(x_cont_traj, 2);
plot3(x_cont_traj(1, i_engage:i_end), x_cont_traj(2, i_engage:i_end), x_cont_traj(3, i_engage:i_end)); hold on
plot3(x_Qtransfer_hill_plus(1, plt_start_i:i_engage + pickle), x_Qtransfer_hill_plus(2, plt_start_i:i_engage + pickle), x_Qtransfer_hill_plus(3, plt_start_i:i_engage + pickle))
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")
title("Hill View of Transfer and Rendezvous")
legend("Rendezvous", "Transfer")
grid on
axis equal

%% Deorbit Trajectory
x_c_deorbit_ECI = propagate_conic(x_c_engage_ECI, ToF_rendezvous, mu_E);
x_d_deorbit = Hill_to_ECI(x_c_deorbit_ECI, optimal_rendezvous.x(1:7, end));
x0_d_keplerian_deorbit = [cartesian_to_keplerian(x_d_deorbit, [0; 0; 1], [1; 0; 0], mu_E); x_d_deorbit(7)];

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
r_a_c = R_E + 660; % [km] periapsis
r_p_c = R_E + 95; % [km] periapsis
e_c = (1 - r_p_c / r_a_c) / (1 + r_p_c / r_a_c); % [] eccentricity
a_c = r_p_c / (1 - e_c); % [km] semi-major axis
i_c = deg2rad(71); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu_c = deg2rad(0); % [rad] true anomaly at epoch

M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);

%%
eta = 0.6;
clear Qtransfer
% Define Q-Law feedback controller: W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
Q_params = struct();
Q_params.W_oe = 1 * ones([5, 1]); % Element weights 
Q_params.eta_a_min = eta; % Minimum absolute efficiency for thrusting instead of coasting
Q_params.eta_r_min = eta; % Minimum relative efficiency for thrusting instead of coasting
Q_params.m = 3;
Q_params.n = 4;
Q_params.r = 2;
Q_params.Theta_rot = 0;

spacecraft_params.m_0 = optimal_rendezvous.x(7, end) + m_debris;
[Qtransfer] = QLaw_transfer(x0_d_keplerian_deorbit(1:6), x0_c_keplerian, mu_E, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 150000, angular_step=deg2rad(20), R_c = 1);

dVs = Qtransfer.delta_V;
ToFs = Qtransfer.dt / 60 / 60 / 24;
if Qtransfer.converged
    fprintf("Q-Law Transfer Converged! Took %.3f Days Using %.3f kg Propellant\n", Qtransfer.dt / 60 / 60 / 24, Qtransfer.delta_m)
else
    fprintf("Q-Law Transfer Failed with %s\n", Qtransfer.errors)
end

plot_Q_transfer(Qtransfer, x0_c_keplerian, x0_d_keplerian, char_star)

%% Plot Eclipsing Function
eclipses = zeros(size(t_keplerian));
for i = 1 : numel(t_keplerian)
    eclipses(i) = check_eclipse(t_keplerian(i), x_keplerian_cartesian(1:3, i), R_E, AU, n_E);
end

figure
plot(t_keplerian / 60 / 60 / 24, eclipses); hold on
grid on
xlabel("Time [days]")
ylabel("Value")
title("Eclipsing vs Time")
subtitle("Positive means eclipsed")


%% Helper Functions
function [] = plot_Q_transfer(Qtransfer, x0_c_keplerian, x0_d_keplerian, char_star, options)
    arguments
        Qtransfer
        x0_c_keplerian
        x0_d_keplerian
        char_star
        options.tolerance = 1e-12
    end
    a_d_0 = @(t, x) zeros([3, 1]); % Disturbance function

    %% Integrate Target Orbit
    tolerances = odeset(RelTol=options.tolerance, AbsTol=options.tolerance);
    [t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_keplerian(x, 1), B_keplerian(x, 1), a_d_0(t,x)), Qtransfer.t / char_star.t, x0_c_keplerian .* [1 / char_star.l, ones([1, 5])]', tolerances);
    x_keplerian_c = x_keplerian_c';
    x_keplerian_c(1, :) = x_keplerian_c(1, :) .* char_star.l; % Redimensionalize
    
    %% Plot Orbit
    x_keplerian_cartesian_d = keplerian_to_cartesian_array(Qtransfer.x_keplerian_mass(1:6, :), [], char_star.mu);
    
    not_coast_colors = interp1(1:numel(Qtransfer.not_coast), double(Qtransfer.not_coast), linspace(1, numel(Qtransfer.not_coast), numel(Qtransfer.t)), 'nearest');
    
    a_c = x0_c_keplerian(1);
    e_c = x0_c_keplerian(2);
    i_c = x0_c_keplerian(3);
    Omega_c = x0_c_keplerian(4);
    omega_c = x0_c_keplerian(5);

    a_d = x0_d_keplerian(1);
    e_d = x0_d_keplerian(2);
    i_d = x0_d_keplerian(3);
    Omega_d = x0_d_keplerian(4);
    omega_d = x0_d_keplerian(5);

    figure
    plot_cartesian_orbit_color_varying(x_keplerian_cartesian_d(:, 1:1:end), not_coast_colors, 3); hold on
    c = colorbar;
    plotOrbit3(Omega_c, i_c, omega_c, a_c * (1 - e_c ^2), e_c, linspace(0, 2 * pi, 1000), "r", 1, 1, [0, 0, 0], 0.1, 2); hold on
    plotOrbit3(Omega_d, i_d, omega_d, a_d * (1 - e_d ^2), e_d, linspace(0, 2 * pi, 1000), "g", 1, 1, [0, 0, 0], 0.1, 2)
    grid on
    earthy(char_star.l, "Earth", 0.5, [0;0;0]); hold on; 
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
    plot_orbit_transfer_histories(Qtransfer.t / 60 / 60, x_keplerian_c' ./ [char_star.l, ones([1, 5])], Qtransfer.x_keplerian_mass(1:6, :)' ./ [char_star.l, ones([1, 5])], Qtransfer.u_cont');
    sgtitle("Q-Law Orbit Transfer with Periapsis Constraint Results")
    
    %% Plot Q Function
    figure
    plot(interp1(1:numel(Qtransfer.t), Qtransfer.t, linspace(1, numel(Qtransfer.t), numel(Qtransfer.Q))) / 60 / 60 / 24, Qtransfer.Q)
    xlabel("Time [days]")
    ylabel("Q-Function []")
    title("Q-Function vs Time")
    yscale("log") % what's the best for plotting this?
    grid on
end


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
    plot(t_hr, u(:, 1) * 1e3, t_hr, u(:, 2) * 1e3, t_hr, u(:, 3) * 1e3); hold on 
    plot(t_hr, vecnorm(u, 2, 2) * 1e3, LineStyle= "--"); hold off
    title("Control")
    xlabel("Time [hr]")
    ylabel("u [N]")
    legend("u_1", "u_2", "u_3", "||u||", Location="northeast")
    grid on
end

function [a_J2_RTN] = J2_perturbation(x_cartesian, mu, r_o, J_2)
    rvec = x_cartesian(1:3);
    r = norm(rvec);
    x = rvec(1);
    y = rvec(2);
    z = rvec(3);

    cons = 1 - 5 * z .^ 2 ./ r .^ 2;

    a_J2 = -3 * mu * J_2 * r_o ^ 2 ./ (2 * r .^ 5) ...
        .* ([cons .* x; cons .* y; (2 + cons) .* z]);

    x_keplerian = cartesian_to_keplerian(x_cartesian, [0; 0; 1], [1; 0; 0], mu);
    e = x_keplerian(2);
    i = x_keplerian(3);
    Omega = x_keplerian(4);
    omega = x_keplerian(5);
    M = x_keplerian(6);

    nu = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(M, e), e);

    a_J2_RTN = cartesian_to_RTN_DCM(i, Omega, omega, nu)' * a_J2;
end

function [a_J2] = J2_perturbation_ECI(x_cartesian, mu, r_o, J_2)
    rvec = x_cartesian(1:3);
    r = norm(rvec);
    x = rvec(1);
    y = rvec(2);
    z = rvec(3);

    cons = 1 - 5 * z .^ 2 ./ r .^ 2;

    a_J2 = -3 * mu * J_2 * r_o ^ 2 ./ (2 * r .^ 5) ...
        .* ([cons .* x; cons .* y; (2 + cons) .* z]);
end



function [f_0] = f0_cartesian(x, mu)
    rvec = x(1:3);
    vvec = x(4:6);
    r = norm(rvec);

    a_cartesian = -mu ./ r .^ 3 * rvec;

    f_0 = [vvec; a_cartesian];
end


function [B] = B_cartesian(x, mu)
    B = [zeros(3); eye(3)];
end


function [x_dot] = gauss_planetary_eqn(f_0, B, a_d)
    x_dot = f_0 + B * a_d;
end