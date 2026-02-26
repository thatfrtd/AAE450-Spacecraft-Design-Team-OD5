%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Combined Transfer from Insertion Orbit and Target Rendezvous
% Author: Travis Hastreiter 
% Created On: 25 February, 2026
% Description: Orbit transfer using Q-Law and rendezvous using SCP in Hill 
% frame.
% Most Recent Change: 25 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = R_E + 660; % [km] semi-major axis
e_c = 1e-3; % [] eccentricity
i_c = deg2rad(71); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu_c = deg2rad(0); % [rad] true anomaly at epoch

M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);

% Initial conditions for spacecraft
r_a_d = R_E + 600; % [km] periapsis
r_p_d = R_E + 96; % [km] periapsis
e_d = (1 - r_p_d / r_a_d) / (1 + r_p_d / r_a_d); % [] eccentricity
a_d = r_p_d / (1 - e_d); % [km] semi-major axis
i_d = deg2rad(69); % [rad] inclination
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

%% Q-Law Transfer to Object (First find time then adjust to phase properly)


%% Relative Orbit Transfer using SCP (start at 1 km away?)


%% Package Output
