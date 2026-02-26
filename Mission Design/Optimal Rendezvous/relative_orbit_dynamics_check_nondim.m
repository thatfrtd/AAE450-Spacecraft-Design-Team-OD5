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
i_c = deg2rad(80); % [rad] inclination
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

%% Convert IC to Cartesian Elements
x_0_cartesian_c = keplerian_to_cartesian(x_keplerian_c, nu0_c, char_star.mu);
[x_0_cartesian_d] = Hill_to_ECI(x_0_cartesian_c, x_0);
x_0_ck = ECI_to_Hill(x_0_cartesian_c, x_0_cartesian_d); % Slightly off...

%% Define Dynamics
f_nonlinear = @(t, x, u, p) nonlinear_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
f_linearized = @(t, x, u, p) linearized_relative_orbit_EoM(t, x, u, p, [x_keplerian_c; spacecraft_params.Isp]);
f_CWH = @(t, x, u, p) CWH_relative_orbit_EoM(t, x, u, p, [a_c; spacecraft_params.Isp]);

%% Propagate
p = [];
u = [0; 1; 0];
u = u / norm(u) * F_max_nd;
a_d = @(t, x) RTN_to_ECI_array(x(1:3), x(4:6)) * u / x(7) / F_max_nd / 1000 * spacecraft_params.F_max;

mu = 1;

a_c_nd = a_c / char_star.l;
alpha_nd = alpha * char_star.l / char_star.t; % alpha has [s / km]

p = a_c_nd * (1 + e_c) * (1 - e_c);
h = sqrt(p * mu);
M = @(t) sqrt(mu / a_c_nd ^ 3) * t;
E = @(t) M(t) + 2 * besselj(1, e_c) * sin(M(t)) + besselj(2, 2 * e_c) * sin(2 * M(t));
r = @(t) a_c_nd * (1 - e_c * cos(E(t)));
v = @(t) sqrt(mu * (2 / r(t) - 1 / a_c_nd));
r_cdot = @(t) sin(acos(h / (r(t) * v(t)))) * v(t);
thetadot = @(t) h / r(t) ^ 2;
thetaddot = @(t) -2 * h  / r(t) ^ 3 * r_cdot(t);

% Propagate relative orbit dynamics
%[~, x_CWH] = ode45(@(t, x) CWH_EoM_manual(t, x, u, mu, a_c_nd, alpha_nd), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
[~, x_CWH] = ode45(@(t, x) f_CWH(t, x, u, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
%[~, x_linearized] = ode45(@(t, x) linearized_relative_orbit_EoM_manual(t, x, u, mu, a_c_nd, thetadot(t), thetaddot(t), alpha_nd), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
[~, x_linearized] = ode45(@(t, x) f_linearized(t, x, u, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
%[~, x_linearized] = ode45(@(t, x) linearized_relative_orbit_EoM_sym(t, x, u, x_keplerian_c ./ [char_star.l; ones([5, 1])], alpha), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
%[~, x_nonlinear] = ode45(@(t, x) nonlinear_relative_orbit_EoM_manual(t, x, u, mu, a_c_nd, r_cdot(t), thetadot(t), alpha_nd), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
[~, x_nonlinear] = ode45(@(t, x) f_nonlinear(t, x, u, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
%[~, x_nonlinear] = ode45(@(t, x) nonlinear_relative_orbit_EoM_sym(t, x, u, x_keplerian_c ./ [char_star.l; ones([5, 1])], alpha), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);

% Propagate 2-body orbit dynamics and convert to relative
[~, x_0_cartesian_c] = ode45(@(t, x) gauss_planetary_eqn(f0_cartesian(x, char_star.mu), B_cartesian(x, char_star.mu), zeros([3, 1])), tspan, x_0_cartesian_c, tolerances);
[~, x_cartesian_d] = ode45(@(t, x) [gauss_planetary_eqn(f0_cartesian(x, char_star.mu), B_cartesian(x, char_star.mu), a_d(t,x)); -alpha * sqrt(u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2)], tspan, x_0_cartesian_d, tolerances);
x_cartesian_hill = ECI_to_Hill(x_0_cartesian_c', x_cartesian_d');

% Redimensionalize
x_CWH = x_CWH' .* nd_scalar;
x_linearized = x_linearized' .* nd_scalar;
x_nonlinear = x_nonlinear' .* nd_scalar;

%% Compare
figure
plot3(x_CWH(1, :), x_CWH(2, :), x_CWH(3, :)); hold on
plot3(x_linearized(1, :), x_linearized(2, :), x_linearized(3, :));
plot3(x_nonlinear(1, :), x_nonlinear(2, :), x_nonlinear(3, :));
plot3(x_cartesian_hill(1, :), x_cartesian_hill(2, :), x_cartesian_hill(3, :)); hold off
grid on
legend("CWH", "Linearized", "Nonlinear", "Cartesian")
axis equal
title("Relative Orbit Position Propagation Comparison")

figure
plot3(x_CWH(4, :), x_CWH(5, :), x_CWH(6, :)); hold on
plot3(x_linearized(4, :), x_linearized(5, :), x_linearized(6, :));
plot3(x_nonlinear(4, :), x_nonlinear(5, :), x_nonlinear(6, :));
plot3(x_cartesian_hill(4, :), x_cartesian_hill(5, :), x_cartesian_hill(6, :)); hold off
grid on
legend("CWH", "Linearized", "Nonlinear", "Cartesian")
axis equal
title("Relative Orbit Velocity Propagation Comparison")

%% Helper Functions
function [x_hill] = ECI_to_Hill(x_c, x_d)
    % Convert from chief and deputy cartesian orbits to relative orbit of
    % deputy w.r.t. chief in Hill frame (chief orbital frame)
    RTN_to_ECI_DCM = RTN_to_ECI_array(x_c(1:3, :), x_c(4:6, :));
    ECI_to_RTN_DCM = pagetranspose(RTN_to_ECI_DCM);

    h_mag = vecnorm(cross(x_c(1:3, :), x_c(4:6, :)));
    r_mag = vecnorm(x_c(1:3, :));
    omega_c = h_mag ./ r_mag .^ 2; % Orbital angular velocity

    x_relative = reshape(x_d(1:6, :) - x_c(1:6, :), 6, 1, []);

    r_hill = pagemtimes(ECI_to_RTN_DCM, x_relative(1:3, :, :));
    v_hill = pagemtimes(ECI_to_RTN_DCM, (x_relative(4:6, :, :)) - cross(reshape([zeros([2, numel(omega_c)]); omega_c], 3, 1, []), r_hill(1:3, :, :)));
    x_hill = [reshape(r_hill, 3, []); 
              reshape(v_hill, 3, []); 
              x_d(7, :)];
end

function [x_d] = Hill_to_ECI(x_c, x_hill)
    % Convert to chief and deputy cartesian orbits from relative orbit of
    % deputy w.r.t. chief in Hill frame (chief orbital frame)
    RTN_to_ECI_DCM = RTN_to_ECI_array(x_c(1:3, :), x_c(4:6, :));

    h_mag = vecnorm(cross(x_c(1:3, :), x_c(4:6, :)));
    r_mag = vecnorm(x_c(1:3, :));
    omega_c = h_mag ./ r_mag .^ 2; % Orbital angular velocity

    r_relative = pagemtimes(RTN_to_ECI_DCM, x_hill(1:3, :, :));
    v_relative = pagemtimes(RTN_to_ECI_DCM, (x_hill(4:6, :, :) + cross([0; 0; omega_c], x_hill(1:3, :, :))));
    x_relative = [reshape(r_relative, 3, []); 
                  reshape(v_relative, 3, [])];

    x_d = [reshape(x_relative + x_c(1:6, :), 6, 1, []);
           x_hill(7, :)]; % Mass
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

function [xdot] = linearized_relative_orbit_EoM_sym(t, x, u, x_keplerian_c, alpha)
    r = x(1:3);
    v = x(4:6);
    m = x(7);

    mu = 1; % Should be nondimensionalized so this is the case

    a_c = x_keplerian_c(1);
    e_c = x_keplerian_c(2);
    M0_c = x_keplerian_c(6);

    p_c = a_c * (1 + e_c) * (1 - e_c);
    h = sqrt(p_c * mu);
    M = sqrt(mu / a_c ^ 3) * t + M0_c;
    % Three term bessel expansion of Kepler's Equation so we can take its derivative
    E = M + 2 * besselj(1, e_c) * sin(M) + besselj(2, 2 * e_c) * sin(2 * M);
    r_c = a_c * (1 - e_c * cos(E));
    v_c = sqrt(mu * (2 / r_c - 1 / a_c));
    r_cdot = sin(acos(h / (r_c * v_c))) * v_c;
    thetastardot = h / r_c ^ 2;
    thetastarddot = -2 * h  / r_c ^ 3 * r_cdot;

    rdot = v;
    vdot = [thetastardot ^ 2 + 2 * mu / r_c ^ 3, thetastarddot, 0; 
           -thetastarddot, thetastardot ^ 2 - mu / r_c ^ 3, 0; 
            0, 0, -mu / r_c ^ 3] * r ...
         + [0, 2 * thetastardot, 0; -2 * thetastardot, 0, 0; 0, 0, 0] * v ...
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

function [xdot] = nonlinear_relative_orbit_EoM_sym(t, x, u, x_keplerian_c, alpha)
    r = x(1:3);
    v = x(4:6);
    m = x(7);

    mu = 1; % Should be nondimensionalized so this is the case

    a_c = x_keplerian_c(1);
    e_c = x_keplerian_c(2);
    M0_c = x_keplerian_c(6);

    p_c = a_c * (1 + e_c) * (1 - e_c);
    h = sqrt(p_c * mu);
    M = sqrt(mu / a_c ^ 3) * t + M0_c;
    % Three term bessel expansion of Kepler's Equation so we can take its derivative
    E = M + 2 * besselj(1, e_c) * sin(M) + besselj(2, 2 * e_c) * sin(2 * M);
    r_c = a_c * (1 - e_c * cos(E));
    v_c = sqrt(mu * (2 / r_c - 1 / a_c));
    r_cdot = sin(acos(h / (r_c * v_c))) * v_c;
    thetastardot = h / r_c ^ 2;

    r_d = norm([r_c; 0; 0] + r);

    rdot = v;
    vdot = [2 * thetastardot * (v(2) - r(2) * r_cdot / r_c) + r(1) * thetastardot ^ 2 + mu / r_c ^ 2 - mu / r_d ^ 3 * (r_c + r(1));
           -2 * thetastardot * (v(1) - r(1) * r_cdot / r_c) + r(2) * thetastardot ^ 2 - mu / r_d ^ 3 * r(2);
           -mu / r_d ^ 3 * r(3)] ...
         + u / m;
    mdot = -alpha * sqrt(u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2);

    xdot = [rdot; vdot; mdot];
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