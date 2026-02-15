function [Omegadot] = J2_RAAN_drift(a, e, i, mu, R, J_2_val)
%J2_RAAN_DRIFT Summary of this function goes here
%   Detailed explanation goes here
arguments
    a % Semimajor axis [km]
    e % Eccentricity []
    i % Inclination [rad]
    mu % Gravitational parameter
    R % Body radius [km]
    J_2_val = 1.0826e-3 % [] Earth J2
end
% Orbit parameters
p = a * (1 - e ^ 2);

% Calculate RAAN change from J2
n = sqrt(mu / a ^ 3); % Mean motion
n_bar = (1 + 3 / 2 * J_2_val * (R / p) ^ 2 * sqrt(1 - e ^ 2) * (1 - 3 / 2 * sin(i) ^ 2)) * n;
Omegadot = (-3 / 2 * J_2_val * (R / p) ^ 2 * cos(i)) * n_bar; % [rad / s]
end