%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% First Order Low Thrust Calculations
% Author: Travis Hastreiter 
% Created On: 13 February, 2026
% Description: Simple calculation of mission delta V using low thrust
% Last Modified On: 13 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumed dynamical parameter values
R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter
J_2_val = 1.0826e-3; % [] Earth J2

% Orbit
a = R_E + 1000;
e = 0.01;
p = a * (1 - e ^ 2);
i = deg2rad(71);

% Calculate RAAN change from J2
n_bar = (1 + 3 / 2 * J_2_val * (R_E / p) ^ 2 * sqrt(1 - e ^ 2) * (1 - 3 / 2 * sin(i) ^ 2)) * sqrt(mu_E / a ^ 3);
Omegadot = rad2deg((-3 / 2 * J_2_val * (R_E / p) ^ 2 * cos(i)) * n_bar);