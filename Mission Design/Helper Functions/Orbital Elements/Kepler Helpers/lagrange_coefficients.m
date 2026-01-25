function [f, g, fdot, gdot, M2, conic_anomalies] = lagrange_coefficients(x1_kep, ToF, mu)
%LAGRANGE_COEFFICIENTS Calculates lagrange coefficients from keplerian state
%   Detailed explanation goes here
arguments
    x1_kep 
    ToF 
    mu 
    %options.conic_anomaly % Maybe already know eccentric or hyperbolic
    %anomaly
end

% Extract Keplerian elements
a = x1_kep(1);
e = x1_kep(2); 
i = x1_kep(3);
Omega =  x1_kep(4);
omega = x1_kep(5);
M1 = x1_kep(6);

% Get final mean motion
mean_motion = sqrt(mu / abs(a) ^ 3);
M2 = M1 + mean_motion * ToF;

% Wrangle Lagrange
if x1_kep(2) < 1 % e < 1 - eccentric
    E1 = mean_to_eccentric_anomaly(M1, e);
    E2 = mean_to_eccentric_anomaly(M2, e);
    delta_E = E2 - E1;

    r1 = eccentric_radius(a, e, E1);
    r2 = eccentric_radius(a, e, E2);

    f = 1 - a / r1 * (1 - cos(delta_E));
    g = ToF + (sin(delta_E) - delta_E) / mean_motion;
    fdot = -mean_motion * a ^ 2 / (r2 * r1) * sin(delta_E);
    gdot = 1 - a / r2 * (1 - cos(delta_E));
        
    conic_anomaly_1 = E1;    
    conic_anomaly_2 = E2;
elseif x1_kep(2) > 1 % e > 1 - hyperbolic
    H1 = mean_to_hyperbolic_anomaly(M1, e);
    H2 = mean_to_hyperbolic_anomaly(M2, e);
    delta_H = H2 - H1;

    r1 = hyperbolic_radius(a, e, H1);
    r2 = hyperbolic_radius(a, e, H2);

    f = 1 - a / r1 * (1 - cosh(delta_H));
    g = ToF - (sinh(delta_H) - delta_H) / mean_motion;
    fdot = -mean_motion * a ^ 2 / (r2 * r1) * sinh(delta_H);
    gdot = 1 - a / r2 * (1 - cosh(delta_H));
        
    conic_anomaly_1 = H1;    
    conic_anomaly_2 = H2;
end
conic_anomalies = [conic_anomaly_1, conic_anomaly_2];

end