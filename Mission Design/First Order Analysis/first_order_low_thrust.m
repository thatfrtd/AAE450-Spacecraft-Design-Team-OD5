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
a = R_E + 3000;
e = 0.01;
i = deg2rad(71);

% Calculate RAAN (Omega) drift
Omega_drift_per_day = rad2deg(J2_RAAN_drift(a, e, i, mu_E, R_E, J_2_val)) * 60 * 60 * 24;