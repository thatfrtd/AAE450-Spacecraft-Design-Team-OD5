%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% HW8 Q1b
% Author: Travis Hastreiter 
% Created On: 11 April, 2025
% Description: Lyapunov orbit rendezvous
% Most Recent Change: 11 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter

% Target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = 10000; % [km] semi-major axis
e_c = 0.3; % [] eccentricity
i_c = deg2rad(50); % [rad] inclination
Omega_c = deg2rad(50); % [rad] right ascension of ascending node
omega_c = deg2rad(40); % [rad] argument of periapsis
nu_c = deg2rad(0); % [rad] true anomaly at epoch

M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);

% Initial conditions for s/c in Earth orbit (in Earth Centered Inertial (ECI) frame)
a_d = 11500; % [km] semi-major axis
e_d = 0.15; % [] eccentricity
i_d = deg2rad(35); % [rad] inclination
Omega_d = deg2rad(50); % [rad] right ascension of ascending node
omega_d = deg2rad(40); % [rad] argument of periapsis
nu_d = deg2rad(0); % [rad] true anomaly at epoch

M_d = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_d, e_d), e_d);
x0_d_keplerian = [a_d; e_d; i_d; Omega_d; omega_d; M_d];
x0_d_cartesian = keplerian_to_cartesian(x0_d_keplerian, nu_d, mu_E);

[l_star, t_star] = nondimensionalized_quantities(a_c, mu_E);

% NEED TO ADD MASS IN STATE
Isp = 3000; % [s]
mass = 800; % [kg]
u_max = 1e-3; % [km / s2]
u_max_star = u_max / l_star * t_star ^ 2;

% Propagation Time 
orbits = 5;
tspan = linspace(0, orbits * period(a_c, mu_E), 1e4);
t_orbits = linspace(0, orbits, numel(tspan));
t_hr = tspan / 60 / 60;

% Integration error tolerance
default_tolerance = 1e-10;

%% a) Simulate orbit - critically damped

a_d_0 = @(t, x) zeros([3, 1]);

% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

% Integrate target orbit
[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_0(t,x)), tspan / t_star, x0_c_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_c = x_keplerian_c ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

% Define Lyapunov feedback controller
K_r = eye(3);
P = 2 * sqrt(K_r);
a_d_cl = @(t, x) cartesian_lyapunov(t, x, interp1(t_keplerian_c, x_keplerian_c, t, "linear","extrap")', K_r, P, 1);

% Integrate spacecraft orbit under Lyapunov feedback controller
[t_keplerian_d, x_keplerian_d] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_cl(t,x)), tspan / t_star, x0_d_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_d = x_keplerian_d ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

%%
u_cl_crit = zeros([numel(tspan), 3]);

for i = 1:numel(tspan)
    u_cl_crit(i, :) = a_d_cl(t_keplerian_d(i), x_keplerian_d(i, :)') * l_star / t_star ^ 2;
end

%% a) Plot 
figure
earthy(R_E, "Earth", 0.5,[0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d, 'b', 0, 1); hold off

title("Cartesian Lyapunov Rendezvous for Critically Damped")
legend("", "Chief", "", "Deputy")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")

%%
plot_relative_orbit_histories(t_hr, x_keplerian_cartesian_c, x_keplerian_cartesian_d, u_cl_crit, "for Critically Damped", t_star, l_star)

plot_relative_orbit(x_keplerian_cartesian_c, x_keplerian_cartesian_d, "for Critically Damped")

%% Bounded Control
%% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_0(t,x)), tspan / t_star, x0_c_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_c = x_keplerian_c ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

K_r = eye(3);
P = 2 * sqrt(K_r);
a_d_cl = @(t, x) cartesian_lyapunov_bounded(t, x, interp1(t_keplerian_c, x_keplerian_c, t, "linear","extrap")', K_r, P, 1, u_max_star);

[t_keplerian_d_bnd, x_keplerian_d] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_cl(t,x)), tspan / t_star, x0_d_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_d_bnd = x_keplerian_d ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

u_cl_crit_bounded = zeros([numel(tspan), 3]);

for i = 1:numel(tspan)
    u_cl_crit_bounded(i, :) = a_d_cl(t_keplerian_d(i), x_keplerian_d(i, :)') * l_star / t_star ^ 2;
end
%%
figure
earthy(R_E, "Earth", 0.5,[0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d_bnd, 'b', 0, 1); hold off

title("Cartesian Lyapunov Rendezvous for Critically Damped and Bounded")
legend("", "Chief", "", "Deputy")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")


%%
plot_relative_orbit_histories(t_hr, x_keplerian_cartesian_c, x_keplerian_cartesian_d_bnd, u_cl_crit_bounded, "for Critically Damped and Bounded", t_star, l_star)

plot_relative_orbit(x_keplerian_cartesian_c, x_keplerian_cartesian_d_bnd, "for Critically Damped and Bounded")

%% a) Simulate orbit - underdamped

a_d_0 = @(t, x) zeros([3, 1]);

% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_0(t,x)), tspan / t_star, x0_c_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_c = x_keplerian_c ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

K_r = eye(3);
P = 2 * sqrt(K_r) * 0.2;
a_d_cl = @(t, x) cartesian_lyapunov(t, x, interp1(t_keplerian_c, x_keplerian_c, t, "linear","extrap")', K_r, P, 1);

[t_keplerian_d, x_keplerian_d] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_cl(t,x)), tspan / t_star, x0_d_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_d = x_keplerian_d ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];
%%
u_cl_crit = zeros([numel(tspan), 3]);

for i = 1:numel(tspan)
    u_cl_crit(i, :) = a_d_cl(t_keplerian_d(i), x_keplerian_d(i, :)') * l_star / t_star ^ 2;
end

%% a) Plot 
figure
earthy(R_E, "Earth", 0.5,[0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d, 'b', 0, 1); hold off

title("Cartesian Lyapunov Rendezvous for Underdamped")
legend("", "Chief", "", "Deputy")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")

%%
plot_relative_orbit_histories(t_hr, x_keplerian_cartesian_c, x_keplerian_cartesian_d, u_cl_crit, "for Underdamped", t_star, l_star)

plot_relative_orbit(x_keplerian_cartesian_c, x_keplerian_cartesian_d, "for Underdamped")

%% Bounded Control
%% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_0(t,x)), tspan / t_star, x0_c_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_c = x_keplerian_c ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

K_r = eye(3);
P = 2 * sqrt(K_r) * 0.2;
a_d_cl = @(t, x) cartesian_lyapunov_bounded(t, x, interp1(t_keplerian_c, x_keplerian_c, t, "linear","extrap")', K_r, P, 1, u_max_star);

[t_keplerian_d_bnd, x_keplerian_d] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_cl(t,x)), tspan / t_star, x0_d_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_d_bnd = x_keplerian_d ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

u_cl_crit_bounded = zeros([numel(tspan), 3]);

for i = 1:numel(tspan)
    u_cl_crit_bounded(i, :) = a_d_cl(t_keplerian_d(i), x_keplerian_d(i, :)') * l_star / t_star ^ 2;
end
%%
figure
earthy(R_E, "Earth", 0.5,[0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d_bnd, 'b', 0, 1); hold off

title("Cartesian Lyapunov Rendezvous for Underdamped and Bounded")
legend("", "Chief", "", "Deputy")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")


%%
plot_relative_orbit_histories(t_hr, x_keplerian_cartesian_c, x_keplerian_cartesian_d_bnd, u_cl_crit_bounded, "for Underdamped and Bounded", t_star, l_star)

plot_relative_orbit(x_keplerian_cartesian_c, x_keplerian_cartesian_d_bnd, "for Underdamped and Bounded")

%% a) Simulate orbit - Overdamped

a_d_0 = @(t, x) zeros([3, 1]);

% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_0(t,x)), tspan / t_star, x0_c_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_c = x_keplerian_c ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

K_r = eye(3);
P = 2 * sqrt(K_r) * 2;
a_d_cl = @(t, x) cartesian_lyapunov(t, x, interp1(t_keplerian_c, x_keplerian_c, t, "linear","extrap")', K_r, P, 1);

[t_keplerian_d, x_keplerian_d] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_cl(t,x)), tspan / t_star, x0_d_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_d = x_keplerian_d ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];
%%
u_cl_crit = zeros([numel(tspan), 3]);

for i = 1:numel(tspan)
    u_cl_crit(i, :) = a_d_cl(t_keplerian_d(i), x_keplerian_d(i, :)') * l_star / t_star ^ 2;
end

%% a) Plot 
figure
earthy(R_E, "Earth", 0.5,[0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d, 'b', 0, 1); hold off

title("Cartesian Lyapunov Rendezvous for Overdamped")
legend("", "Chief", "", "Deputy")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")

%%
plot_relative_orbit_histories(t_hr, x_keplerian_cartesian_c, x_keplerian_cartesian_d, u_cl_crit, "for Overdamped", t_star, l_star)

plot_relative_orbit(x_keplerian_cartesian_c, x_keplerian_cartesian_d, "for Overdamped")

%% Bounded Control
%% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

[t_keplerian_c, x_keplerian_c] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_0(t,x)), tspan / t_star, x0_c_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_c = x_keplerian_c ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

K_r = eye(3);
P = 2 * sqrt(K_r) * 1.5;
a_d_cl = @(t, x) cartesian_lyapunov_bounded(t, x, interp1(t_keplerian_c, x_keplerian_c, t, "linear","extrap")', K_r, P, 1, u_max_star);

[t_keplerian_d_bnd, x_keplerian_d] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, 1), B_cartesian(x, 1), a_d_cl(t,x)), tspan / t_star, x0_d_cartesian .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])]', tolerances);
x_keplerian_cartesian_d_bnd = x_keplerian_d ./ [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

u_cl_crit_bounded = zeros([numel(tspan), 3]);

for i = 1:numel(tspan)
    u_cl_crit_bounded(i, :) = a_d_cl(t_keplerian_d(i), x_keplerian_d(i, :)') * l_star / t_star ^ 2;
end
%%
figure
earthy(R_E, "Earth", 0.5,[0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian_c, 'r', 0, 1); hold on
plot_cartesian_orbit(x_keplerian_cartesian_d_bnd, 'b', 0, 1); hold off

title("Cartesian Lyapunov Rendezvous for Overdamped and Bounded")
legend("", "Chief", "", "Deputy")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")


%%
plot_relative_orbit_histories(t_hr, x_keplerian_cartesian_c, x_keplerian_cartesian_d_bnd, u_cl_crit_bounded, "for Overdamped Bounded", t_star, l_star)

plot_relative_orbit(x_keplerian_cartesian_c, x_keplerian_cartesian_d_bnd, "for Overdamped and Bounded")


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
    B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)';
end

function [x_dot] = gauss_planetary_eqn(f_0, B, a_d)
    x_dot = f_0 + B * a_d;
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

function [a_cart_lyapunov] = cartesian_lyapunov(t, x, x_c, K_r, P, mu)
    delta_x = x - x_c;
    delta_a = f0_cartesian(x, mu) - f0_cartesian(x_c, mu);

    a_cart_lyapunov = (-K_r * delta_x(1:3) ...
        - P * delta_x(4:6) ...
        - delta_a(4:6));
end


function [a_cart_lyapunov_bounded] = cartesian_lyapunov_bounded(t, x, x_c, K_r, P, mu, u_max)
    a_cart_lyapunov = cartesian_lyapunov(t, x, x_c, K_r, P, mu);

    if norm(a_cart_lyapunov) > u_max
        a_cart_lyapunov_bounded = a_cart_lyapunov / norm(a_cart_lyapunov) * u_max;
    else
        a_cart_lyapunov_bounded = a_cart_lyapunov;
    end
end

function [f_0] = f0_cartesian(x, mu)
    rvec = x(1:3);
    vvec = x(4:6);
    r = norm(rvec);

    a_cartesian = -mu ./ r .^ 3 * rvec;

    f_0 = [vvec; a_cartesian];
end

function [l_star, t_star] = nondimensionalized_quantities(a_c, mu)
    l_star = a_c; % [km]

    t_star = sqrt(l_star ^ 3 / mu); % [s]
end

function [] = plot_relative_orbit(x_c, x_d, title_tag)
    delta_x = x_c - x_d;

    figure

    plot3(delta_x(:, 1), delta_x(:, 2), delta_x(:, 3)); hold on
    scatter3(delta_x(1, 1), delta_x(1, 2), delta_x(1, 3), Marker="o"); hold on
    scatter3(delta_x(end, 1), delta_x(end, 2), delta_x(end, 3), Marker="diamond"); hold on
    scatter3(0, 0, 0, Marker="diamond"); hold off
    xlabel("\delta r_1 [km]")
    ylabel("\delta r_2 [km]")
    zlabel("\delta r_3 [km]")
    legend("", "initial", "final", "target")
    title("Lyapunov Controller Rendezvous Deputy Relative Position " + title_tag)
    grid on
end

function [] = plot_relative_orbit_histories(t_hr, x_c, x_d, u, title_tag, t_star, l_star)
    delta_x = x_c - x_d;
    
    figure
    tiledlayout(2, 2)
    
    nexttile
    plot(t_hr, delta_x(:, 1), t_hr, delta_x(:, 2), t_hr, delta_x(:, 3)); hold on 
    plot(t_hr, vecnorm(delta_x(:, 1:3), 2, 2), LineStyle="--"); hold off
    title("Position Error")
    xlabel("Time [hr]")
    ylabel("\delta r [km]")
    legend("\delta r_1", "\delta r_2", "\delta r_3", "||\delta r||", Location="northeast")
    grid on
    
    nexttile
    plot(t_hr, delta_x(:, 4), t_hr, delta_x(:, 5), t_hr, delta_x(:, 6)); hold on 
    plot(t_hr, vecnorm(delta_x(:, 4:6), 2, 2), LineStyle="--"); hold off
    title("Velocity Error")
    xlabel("Time [hr]")
    ylabel("\delta v [km / s]")
    legend("\delta v_1", "\delta v_2", "\delta v_3", "||\delta v||", Location="northeast")
    grid on
    
    nexttile
    plot(t_hr, u(:, 1), t_hr, u(:, 2), t_hr, u(:, 3)); hold on 
    plot(t_hr, vecnorm(u, 2, 2), LineStyle= "--"); hold off
    title("Control")
    xlabel("Time [hr]")
    ylabel("u [km / s2]")
    legend("u_1", "u_2", "u_3", "||u||", Location="northeast")
    grid on

    delta_x = delta_x .* [1 / l_star * ones([1, 3]), t_star / l_star * ones([1, 3])];

    nexttile
    plot(t_hr, vecnorm(delta_x(:, [1, 4]), 2, 2), t_hr, vecnorm(delta_x(:, [2, 5]), 2, 2), t_hr, vecnorm(delta_x(:, [3, 6]), 2, 2)); hold on 
    plot(t_hr, vecnorm(delta_x, 2, 2), LineStyle="--"); hold off
    title("Delta")
    xlabel("Time [hr]")
    ylabel("\Delta")
    legend("\Delta_1", "\Delta_2", "\Delta_3", "||\Delta||", Location="northeast")
    yscale("log")
    grid on
    
    sgtitle("Cartesian Lyapunov Controller Rendezvous Time Histories " + title_tag)
end

function [B] = B_cartesian(x, mu)
    B = [zeros(3); eye(3)];
end