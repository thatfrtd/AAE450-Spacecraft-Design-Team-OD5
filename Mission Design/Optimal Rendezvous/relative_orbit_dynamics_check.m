%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Check of Relative Orbit Dynamics (CWH, Linearized, Nonlinear) 
% Author: Travis Hastreiter 
% Created On: 25 February, 2026
% Description: Check of relative orbit dynamics (CWH, Linearized, Nonlinear) 
% by propagating with and without control and looking at divergence rates.
% Also, compare against directly propagating orbits in inertial frame
% Last Modified On: 9 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Data
data_table = readtable("Full_Detumble_History.csv");

%%
char_star = load_charecteristic_values_Earth();
nd_scalar = [char_star.l * ones([3, 1]); char_star.v * ones([3, 1]); char_star.m];

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 4100; % [s]
spacecraft_params.m_0 = 1500; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.25; % [N]
F_max_nd = spacecraft_params.F_max / 1000 / char_star.F; % F_max in N, char_star.F in kN

alpha = 1 / (spacecraft_params.Isp * 9.81e-3);

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = char_star.l + 800; % [km] semi-major axis
e_c = 0.000; % [] eccentricity
i_c = deg2rad(98); % [rad] inclination
Omega_c = deg2rad(0); % [rad] right ascension of ascending node
omega_c = deg2rad(0); % [rad] argument of periapsis
nu0_c = deg2rad(0); % [rad] true anomaly at epoch
M0_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu0_c, e_c), e_c);
x_keplerian_c = [a_c; e_c; i_c; Omega_c; omega_c; M0_c];

n = sqrt(char_star.mu / a_c ^ 3); % [rad / s]

% Rendezvous time
tf = 7200; % / char_star.t; % [s] (nondimensionalized)
tspan = linspace(0, tf, 1000);

b = sqrt(40^2/2);
c = sqrt(10^2/2);
aROE = [0; % [km] delta semimajor axis
        0; % [km] delta lambda
        b*1e-3; % [km] delta e_x
        b*1e-3; % [km] delta e_y
        c*1e-3; % [km] delta i_x 
        c*1e-3]; % [km] delta i_y

cart_ROE = ROE_to_cart_matrix(n, 0) * aROE;
%%
ROE_to_cart_matrix(n, 0) * cart_to_ROE_matrix(n, 0)
%%
cart_to_ROE_matrix(n, 0) * [data_table.chaser_x(1); data_table.chaser_y(1); data_table.chaser_z(1); data_table.chaser_vx(1); data_table.chaser_vy(1); data_table.chaser_vz(1)]
%%

% Initial conditions for spacecraft - specify orbit instead?
% r_0 = [0; -1; 0]; % [km]
% v_0 = [0e-3; 4e-3; -0.024]; % [km / s]
% x_0 = [cart_ROE; spacecraft_params.m_0];

r_0 = [data_table.chaser_x(1); data_table.chaser_y(1); data_table.chaser_z(1)]*1e-3; % [km]
v_0 = [data_table.chaser_vx(1); data_table.chaser_vy(1); data_table.chaser_vz(1)]*1e-3; % [km / s]
x_0 = [cart_ROE; spacecraft_params.m_0];

% Algorithm parameters
default_tolerance = 1e-13;
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

%%
aROE_ck = cart_to_ROE_matrix(n, 0) * [r_0; v_0];

%%


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
u = [0; 0; 0];
%u = u / norm(u) * spacecraft_params.F_max / 1000;
a_d_0 = @(t, x) [0; 0; 0];
a_d = @(t, x) RTN_to_ECI_array(x(1:3), x(4:6)) * u / x(7);

mu = char_star.mu;

p = a_c * (1 + e_c) * (1 - e_c);
h = sqrt(p * mu);
M = @(t) sqrt(mu / a_c ^ 3) * t + M0_c;
E = @(t) M(t) + 2 * besselj(1, e_c) * sin(M(t)) + besselj(2, 2 * e_c) * sin(2 * M(t));
r = @(t) a_c * (1 - e_c * cos(E(t)));
v = @(t) sqrt(mu * (2 / r(t) - 1 / a_c));
r_cdot = @(t) sin(acos(h / (r(t) * v(t)))) * v(t);
thetadot = @(t) h / r(t) ^ 2;
thetaddot = @(t) -2 * h  / r(t) ^ 3 * r_cdot(t);

% Propagate relative orbit dynamics
[~, x_CWH] = ode45(@(t, x) CWH_EoM_manual(t, x, u, mu, a_c, alpha), tspan, x_0, tolerances);
%[~, x_CWH] = ode45(@(t, x) f_CWH(t, x, u * 1000 / char_star.F, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);
%[~, x_linearized] = ode45(@(t, x) linearized_relative_orbit_EoM_manual(t, x, u, mu, a_c, thetadot(t), thetaddot(t), alpha), tspan, x_0, tolerances);
[~, x_linearized] = ode45(@(t, x) f_linearized(t / char_star.t, x ./ nd_scalar, u / char_star.F, p) / char_star.t .* nd_scalar, tspan, x_0, tolerances);
%[~, x_nonlinear] = ode45(@(t, x) nonlinear_relative_orbit_EoM_manual(t, x, u, mu, a_c, r_cdot(t), thetadot(t), alpha), tspan, x_0, tolerances);
[~, x_nonlinear] = ode45(@(t, x) f_nonlinear(t / char_star.t, x ./ nd_scalar, u / char_star.F, p) / char_star.t .* nd_scalar, tspan, x_0, tolerances);
%[~, x_nonlinear] = ode45(@(t, x) f_nonlinear(t, x, u, p), tspan / char_star.t, x_0 ./ nd_scalar, tolerances);

% Propagate 2-body orbit dynamics and convert to relative
[~, x_0_cartesian_c] = ode45(@(t, x) gauss_planetary_eqn(f0_cartesian(x, char_star.mu), B_cartesian(x, char_star.mu), a_d_0(t,x)), tspan, x_0_cartesian_c, tolerances);
[~, x_cartesian_d] = ode45(@(t, x) [gauss_planetary_eqn(f0_cartesian(x, char_star.mu), B_cartesian(x, char_star.mu), a_d(t,x)); -alpha * sqrt(u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2)], tspan, x_0_cartesian_d, tolerances);
x_cartesian_hill = ECI_to_Hill(x_0_cartesian_c', x_cartesian_d');

% Redimensionalize
x_CWH = x_CWH' .* 1;%nd_scalar;
x_linearized = x_linearized' .* 1;%nd_scalar;
x_nonlinear = x_nonlinear' .* 1;%nd_scalar;

%% Compare
figure
plot3(x_CWH(1, :), x_CWH(2, :), x_CWH(3, :)); hold on
scatter3(data_table.chaser_x(1:100:end)*1e-3, data_table.chaser_y(1:100:end)*1e-3, data_table.chaser_z(1:100:end)*1e-3); hold on
plot3(x_linearized(1, :), x_linearized(2, :), x_linearized(3, :)); hold on
plot3(x_nonlinear(1, :), x_nonlinear(2, :), x_nonlinear(3, :)); hold on
plot3(x_cartesian_hill(1, :), x_cartesian_hill(2, :), x_cartesian_hill(3, :)); hold off
grid on
legend("CWH", "Atharva", "Linearized", "Nonlinear", "Cart")
axis equal
xlabel("R")
ylabel("T")
zlabel("N")
title("Position Comparison")

figure
plot3(x_CWH(4, :), x_CWH(5, :), x_CWH(6, :)); hold on
plot3(x_linearized(4, :), x_linearized(5, :), x_linearized(6, :)); hold on
plot3(x_nonlinear(4, :), x_nonlinear(5, :), x_nonlinear(6, :)); hold on
plot3(x_cartesian_hill(4, :), x_cartesian_hill(5, :), x_cartesian_hill(6, :)); hold off
grid on
legend("CWH", "Linearized", "Nonlinear", "Cart")
axis equal
xlabel("R")
ylabel("T")
zlabel("N")
title("Velocity Comparison")

%%
r_cont_sol = x_CWH(1:3, :);

output_array = [tspan', r_cont_sol'];
state_names = ["r_1", "r_2", "r_3"];
output_names = ["Time", state_names];

output_table = array2table(output_array, VariableNames = output_names);
writetable(output_table,"./Animation/inspection_orbit.csv")

%% Helper Functions
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


function [x_star, tstar] = nondimensionalize_from_cartesian(x_cartesian, t, p_star, mu_barycenter, d)
    l_star = p_star(1);
    t_star = p_star(3);
    
    tstar = t / t_star;

    omega_star_NR = sqrt(mu_barycenter / d ^ 3) * t_star; % [rad / s] Earth-Moon mean motion about barycenter
    C_RN = make_R(omega_star_NR * tstar, 3);

    rvec_star_N = x_cartesian(1:3)' / l_star;
    vvec_star_N = x_cartesian(4:6)' / l_star * t_star;

    rvec_star_R = C_RN * rvec_star_N;
    % Transport theorem
    vvec_star_R = C_RN * (vvec_star_N + cross([0; 0; omega_star_NR], rvec_star_N));

    x_star = [rvec_star_R; vvec_star_R];
end

function [x_star_array, tstar_array] = nondimensionalize_from_cartesian_array(x_cartesian_array, t_array, p_star, mu_barycenter, d)
    x_star_array = zeros(size(x_cartesian_array, 1), 6);
    tstar_array = zeros(size(x_cartesian_array, 1), 1);

    for ind = 1:size(x_cartesian_array, 1)
        x_cartesian = x_cartesian_array(ind, :);

        [x_star_array(ind, :), tstar_array(ind)] = nondimensionalize_from_cartesian(x_cartesian, t_array(ind), p_star, mu_barycenter, d);
    end
end

function [x_cartesian_N, t] = dimensionalize_to_cartesian(x_star, tstar, p_star, mu_barycenter, d)
    l_star = p_star(1);
    t_star = p_star(3);

    t = tstar * t_star;

    omega_RN = -sqrt(mu_barycenter / d ^ 3); % [rad / s] Earth-Moon mean motion about barycenter
    C_NR = make_R(omega_RN * t, 3);

    rvec_cartesian_R = x_star(1:3)' * l_star;
    vvec_cartesian_R = x_star(4:6)' * l_star / t_star;

    rvec_cartesian_N = C_NR * rvec_cartesian_R;
    % Transport theorem
    vvec_cartesian_N = C_NR * (vvec_cartesian_R + cross([0; 0; -omega_RN], rvec_cartesian_R));

    x_cartesian_N = [rvec_cartesian_N; vvec_cartesian_N];
end

function [x_cartesian_N, t_array] = dimensionalize_to_cartesian_array(x_star_array, tstar_array, p_star, mu_barycenter, d)
    x_cartesian_N = zeros(size(x_star_array, 1), 6);
    t_array = zeros(size(x_star_array, 1), 1);

    for ind = 1:size(x_star_array, 1)
        x_star = x_star_array(ind, :);

        [x_cartesian_N(ind, :), t_array(ind)] = dimensionalize_to_cartesian(x_star, tstar_array(ind), p_star, mu_barycenter, d);
    end
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