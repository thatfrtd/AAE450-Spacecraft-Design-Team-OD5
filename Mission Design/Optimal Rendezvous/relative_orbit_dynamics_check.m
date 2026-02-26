%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Check of Relative Orbit Dynamics (CWH, Linearized, Nonlinear) 
% Author: Travis Hastreiter 
% Created On: 25 February, 2026
% Description: Check of relative orbit dynamics (CWH, Linearized, Nonlinear) 
% by propagating with and without control and looking at divergence rates.
% Also, compare against directly propagating orbits in inertial frame
% Last Modified On: 25 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

char_star = load_charecteristic_values_Earth();
nd_scalar = [char_star.l * ones([3, 1]); char_star.v * ones([3, 1]); char_star.m];

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 3000; % [s]
spacecraft_params.m_0 = 800; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.1; % [N]
F_max_nd = spacecraft_params.F_max / 1000 / char_star.F; % F_max in N, char_star.F in kN

alpha = 1 / (spacecraft_params.Isp * 9.81e-3);

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = 6728; % [km] semi-major axis
e_c = 0.01; % [] eccentricity
i_c = deg2rad(10); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu0_c = deg2rad(0); % [rad] true anomaly at epoch
M0_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu0_c, e_c), e_c);
x_keplerian_c = [a_c; e_c; i_c; Omega_c; omega_c; M0_c];

n = sqrt(char_star.mu / a_c ^ 3); % [rad / s]

% Rendezvous time
tf = 3600; % / char_star.t; % [s] (nondimensionalized)
tspan = linspace(0, tf, 1000);

% Initial conditions for spacecraft - specify orbit instead?
r_0 = [-0.5; -0.5; 0.2]; % [km]
v_0 = [0e-3; 1e-3; 0]; % [km / s]
x_0 = [r_0; v_0; spacecraft_params.m_0];

% Algorithm parameters
default_tolerance = 1e-13;
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

%% Define Dynamics
f_nonlinear = @(t, x, u, p) nonlinear_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
f_linearized = @(t, x, u, p) linearized_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
f_CWH = @(t, x, u, p) CWH_relative_orbit_EoM(t, x, u, p, [a_c; spacecraft_params.Isp]);

%% Propagate
p = [];
u = [0; 1e-4; 1e-4];
a_d = @(t, x) [0; 0; 0];

mu = char_star.mu;

p = a_c * (1 + e_c) * (1 - e_c);
h = sqrt(p * mu);
theta_0 = 0;
M = @(t) sqrt(mu / a_c ^ 3) * tf;
theta = @(t) theta_0 + eccentric_to_true_anomaly(mean_to_eccentric_anomaly(M(t), e_c), e_c);
E = @(t) M(t) + 2 * besselj(1, e_c) * sin(M(t)) + besselj(2, 2 * e_c) * sin(2 * M(t));
r = @(t) a_c * (1 - e_c * cos(E(t)));
v = @(t) sqrt(mu * (2 / r(t) - 1 / a_c));
r_cdot = @(t) sin(acos(h / (r(t) * v(t)))) * v(t);
thetadot = @(t) h / r(t) ^ 2;
thetaddot = @(t) -2 * h  / r(t) ^ 3 * r_cdot(t);

% Propagate relative orbit dynamics
[~, x_CWH] = ode45(@(t, x) CWH_EoM_manual(t, x, u, char_star.mu, a_c, alpha), tspan, x_0, tolerances);
%[~, x_CWH] = ode45(@(t, x) f_CWH(t, x, u * 1000 / char_star.F, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
[~, x_linearized] = ode45(@(t, x) linearized_relative_orbit_EoM_manual(t, x, u, char_star.mu, a_c, thetadot(t), thetaddot(t), alpha), tspan, x_0, tolerances);
%[~, x_linearized] = ode45(@(t, x) f_linearized(t, x, u, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
[~, x_nonlinear] = ode45(@(t, x) nonlinear_relative_orbit_EoM_manual(t, x, u, char_star.mu, a_c, r_cdot(t), thetadot(t), alpha), tspan, x_0, tolerances);
%[~, x_nonlinear] = ode45(@(t, x) f_nonlinear(t, x, u, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);

% Propagate 2-body orbit dynamics and convert to relative
%[~, x_cartesian_c] = ode45(@(t, x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d(t,x)), tspan, x_0_cartesian_c, tolerances = tolerances);
%[~, x_cartesian_d] = ode45(@(t, x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d(t,x)), tspan, x_0_cartesian_d, tolerances = tolerances);
%x_cartesian = ECI_to_Hill(x_cartesian_c, x_cartesian_d);

% Redimensionalize
x_CWH = x_CWH' .* 1;%nd_scalar;
x_linearized = x_linearized' .* 1;%nd_scalar;
x_nonlinear = x_nonlinear' .* 1;%nd_scalar;

%% Compare
figure
plot3(x_CWH(1, :), x_CWH(2, :), x_CWH(3, :)); hold on
plot3(x_linearized(1, :), x_linearized(2, :), x_linearized(3, :)); hold on
plot3(x_nonlinear(1, :), x_nonlinear(2, :), x_nonlinear(3, :)); hold off
grid on
legend("CWH", "Linearized", "Nonlinear")
axis equal

%% Helper Functions
function [x_hill] = ECI_to_Hill(x_c, x_d)
    % Convert from chief and deputy cartesian orbits to relative orbit of
    % deputy w.r.t. chief in Hill frame (chief orbital frame)
    RTN_to_ECI_DCM = RTN_to_ECI_array(x_c(1:3, :), x_c(4:6));
    ECI_to_RTN_DCM = pagetranspose(RTN_to_ECI_DCM);

    x_relative = reshape(x_d(1:6, :) - x_c(1:6, :), 6, 1, []);
    x_hill = [pagemtimes(ECI_to_RTN_DCM, x_relative); x_d(7, :)];
end

function [xdot] = CWH_EoM_manual(t, x, u, mu, r_c, alpha)
    r = x(1:3);
    v = x(4:6);
    m = x(7);

    n = sqrt(mu / r_c ^ 3);

    rdot = v;
    vdot = [3 * n ^ 2, 0, 0; 0, 0, 0; 0, 0, -n ^ 2]  * r ...
        +  [0, 2 * n, 0; -2 * n, 0, 0; 0, 0, 0] * v ...
        + u / m;
    mdot = -alpha * sqrt(u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2);

    xdot = [rdot; vdot; mdot];
end

function [xdot] = linearized_relative_orbit_EoM_manual(t, x, u, mu, r_c, thetadot, thetaddot, alpha)
    r = x(1:3);
    v = x(4:6);
    m = x(7);

    rdot = v;
    vdot = [thetadot ^ 2 + 2 * mu / r_c ^ 3, thetaddot, 0; 
           -thetaddot, thetadot ^ 2 - mu / r_c ^ 3, 0; 
            0, 0, -mu / r_c ^ 3] * r ...
         + [0, 2 * thetadot, 0; -2 * thetadot, 0, 0; 0, 0, 0] * v ...
         + u / m;
    mdot = -alpha * sqrt(u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2);

    xdot = [rdot; vdot; mdot];
end

function [xdot] = nonlinear_relative_orbit_EoM_manual(t, x, u, mu, r_c, r_cdot, nudot, alpha)
    r = x(1:3);
    v = x(4:6);
    m = x(7);

    r_d = norm([r_c; 0; 0] + r);

    rdot = v;
    vdot = [2 * nudot * (v(2) - r(2) * r_cdot / r_c) + r(1) * nudot ^ 2 + mu / r_c ^ 2 - mu / r_d ^ 3 * (r_c + r(1));
           -2 * nudot * (v(1) - r(1) * r_cdot / r_c) + r(2) * nudot ^ 2 - mu / r_d ^ 3 * r(2);
           -mu / r_d ^ 3 * r(3)] ...
         + u / m;
    mdot = -alpha * sqrt(u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2);

    xdot = [rdot; vdot; mdot];
end