%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Author: Travis Hastreiter 
% Created On: 13 February, 2026
% Description: Simulate orbits with J2 perturbations, thrust, and drag
% Most Recent Change: 13 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assumed dynamical parameter values
R_E = 6378.1; % [km] Earth radius
mu_E = 398600; % [km3 / s2] Earth gravitational parameter
J_2_val = 1.0826e-3; % [] Earth J2

% Spacecraft parameters
m = 800; % [kg]

% Initial conditions for s/c in Earth orbit (in Earth Centered Inertial (ECI) frame)
r_a_0 = R_E + 600; % [km] periapsis
r_p_0 = R_E + 96; % [km] periapsis
e_0 = (1 - r_p_0 / r_a_0) / (1 + r_p_0 / r_a_0); % [] eccentricity
a_0 = r_p_0 / (1 - e_0); % [km] semi-major axis
i_0 = deg2rad(80); % [rad] inclination
Omega_0 = deg2rad(30); % [rad] right ascension of ascending node
omega_0 = deg2rad(20); % [rad] argument of periapsis
nu_0 = deg2rad(0); % [rad] true anomaly at epoch

M_0 = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_0, e_0), e_0);
x0_keplerian = [a_0; e_0; i_0; Omega_0; omega_0; M_0];
x0_cartesian = keplerian_to_cartesian(x0_keplerian, nu_0, mu_E);

% Propagation Time 
orbits = 40;
tspan = linspace(0, orbits * period(a_0, mu_E), 1e5);
t_orbits = linspace(0, orbits, numel(tspan));
t_hr = tspan / 60 / 60;

% Integration error tolerance
default_tolerance = 1e-6;

% Control
F_mag = 1; % [N]
accel_thrust = F_mag / m / 1000; % [km / s2]

% Drag
A_over_m = pi * 1 ^ 2 / m * 1e-6; % [km2 / kg] spacecraft area to mass ratio
C_D = 2.31; % [] drag coefficient

stop_altitude = 1; % [km]

%%
t0_datetime = datetime('2030-01-01','InputFormat','yyyy-MM-dd');

%% a) Simulate orbit under Earth J2 perturbation and thrust
a_d_J2 = @(t,x) J2_perturbation(x, mu_E, R_E, J_2_val);

thrust_DCM = angle2dcm(deg2rad(0), pi / 2 - deg2rad(0), deg2rad(0),"ZYX");

a_d_thrust = @(t,x) [accel_thrust * sind(0); ...
                     accel_thrust * cosd(0); ... 
                     accel_thrust * sind(0)]; % [km / s2] Orbital RTN frame (Radial-Theta-Normal frame)

a_d = @(t,x) a_d_J2(t,x) + cartesian_to_RTN_DCM_from_cart(x, mu_E) * a_d_thrust(t,x) ...
                         + drag_perturbation_nrlmsise(t, x, mu_E, C_D, A_over_m, t0_datetime, false);
                         %+ drag_perturbation_simple(t, x, mu_E, R_E, C_D, A_over_m);

% Simulate 
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance, Events = @(t, x) altitude_event_function(t, x, R_E, stop_altitude, t0_datetime));
    
[t_keplerian,x_keplerian_cartesian] = ode45(@(t,x) gauss_planetary_eqn(f0_cartesian(x, mu_E), B_cartesian(x, mu_E), a_d(t,x)), tspan, x0_cartesian, tolerances);
x_keplerian_cartesian = x_keplerian_cartesian';
x_keplerian = cartesian_to_keplerian_array(x_keplerian_cartesian, [0; 0; 1], [1; 0; 0], mu_E);

%%
% clear as
% ind = 1;
% for i = 1 :10: numel(t_keplerian)
%     [a_drag_RTN] = drag_perturbation_nrlmsise(t_keplerian(i), x_keplerian_cartesian(:, i), mu_E, C_D, A_over_m, t0_datetime, false);
%     as(ind) = norm(a_drag_RTN);
%     ind = ind + 1;
% end
%%
% plot(t_keplerian(1:10:end) / 60 / 60, as)
% yscale("log")
% %%
% plot(vecnorm(x_keplerian_cartesian(4:6, :)))

%% Delta V
Isp = 2000;
g_0 = 9.81;
delta_m = 1 / (Isp * g_0) * F_mag * tspan(end)
dV_spiral = sqrt(mu_E / a_0) - sqrt(mu_E / x_keplerian(1, end))
dV_sim = Isp * g_0 * log(m / (m - delta_m)) / 1000

dt_est_hr = dV_spiral * 1000 / (F_mag / m) / 60 / 60

%% a) Plot 
figure
earthy(R_E, "Earth", 1, [0;0;0]); hold on;
axis equal

plot_cartesian_orbit(x_keplerian_cartesian', 'r', 0, 1);

title("Orbit Propagated with J2 with Constant Thrust")
xlabel("X [km]")
ylabel("Y [km]")
zlabel("Z [km]")

%%
figure
tiledlayout(2, 3)

nexttile
plot(t_keplerian, x_keplerian(1, :) .* (1 - x_keplerian(2, :) .* cos(mean_to_eccentric_anomaly(x_keplerian(6, :), x_keplerian(2, :)))) - R_E)
title("Radius vs Time")
grid on

nexttile
plot(t_keplerian, x_keplerian(2, :))
title("Eccentricity vs Time")
grid on

nexttile
plot(t_keplerian, rad2deg(x_keplerian(3, :)))
title("Inclination vs Time")
grid on

nexttile
plot(t_keplerian, rad2deg(x_keplerian(4, :)))
title("RAAN vs Time")
grid on

nexttile
plot(t_keplerian, rad2deg(x_keplerian(5, :)))
title("Argument of Periapsis Axis vs Time")
grid on

nexttile
plot(t_keplerian, wrapTo360(rad2deg(x_keplerian(6, :))))
title("Mean Anomaly vs Time")
grid on

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
    %B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)';
end

function [x_dot] = gauss_planetary_eqn(f_0, B, a_d)
    x_dot = f_0 + B * a_d;
end

function [a_J2] = J2_perturbation(x_cartesian, mu, r_o, J_2)
    rvec = x_cartesian(1:3);
    r = norm(rvec);
    x = rvec(1);
    y = rvec(2);
    z = rvec(3);

    cons = 1 - 5 * z .^ 2 ./ r .^ 2;

    a_J2 = -3 * mu * J_2 * r_o ^ 2 ./ (2 * r .^ 5) ...
        .* ([cons .* x; cons .* y; (2 + cons) .* z]);
    % 
    % x_keplerian = cartesian_to_keplerian(x_cartesian, [0; 0; 1], [1; 0; 0], mu);
    % e = x_keplerian(2);
    % i = x_keplerian(3);
    % Omega = x_keplerian(4);
    % omega = x_keplerian(5);
    % M = x_keplerian(6);
    % 
    % nu = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(M, e), e);
    % 
    % a_J2_RTN = cartesian_to_RTN_DCM(i, Omega, omega, nu)' * a_J2;
end

function [a_drag] = drag_perturbation_simple(t, x_cartesian, mu, r_o, C_D, A_over_m)
    rvec = x_cartesian(1:3);
    vvec = x_cartesian(4:6);
    r = norm(rvec);
    v = norm(vvec);
    alt = r - r_o;

    rho = 1.45*1.02e7 * alt ^ -7.172 * 1e9; % Formula (up to 1000 km ish) from 
     
    a_drag = -1 / 2 * rho * v * vvec * A_over_m * C_D;
    % 
    % % Convert to RTN
    % x_keplerian = cartesian_to_keplerian(x_cartesian, [0; 0; 1], [1; 0; 0], mu);
    % e = x_keplerian(2);
    % i = x_keplerian(3);
    % Omega = x_keplerian(4);
    % omega = x_keplerian(5);
    % M = x_keplerian(6);
    % 
    % nu = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(M, e), e);
    % 
    % a_drag_RTN = cartesian_to_RTN_DCM(i, Omega, omega, nu)' * a_drag;
end

function [DCM_ECI_RTN] = cartesian_to_RTN_DCM_from_cart(x_cartesian, mu)
    x_keplerian = cartesian_to_keplerian(x_cartesian, [0; 0; 1], [1; 0; 0], mu);
    e = x_keplerian(2);
    i = x_keplerian(3);
    Omega = x_keplerian(4);
    omega = x_keplerian(5);
    M = x_keplerian(6);

    nu = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(M, e), e);

    DCM_ECI_RTN = cartesian_to_RTN_DCM(i, Omega, omega, nu);
end

function [a_drag] = drag_perturbation_nrlmsise(t, x_cartesian, mu, C_D, A_over_m, t0_datetime, include_Oxygen)
    arguments
        t
        x_cartesian
        mu
        C_D
        A_over_m
        t0_datetime
        include_Oxygen = false
    end

    rvec = x_cartesian(1:3);
    vvec = x_cartesian(4:6);
    r = norm(rvec);
    v = norm(vvec);

    t_datetime = t0_datetime + seconds(t);

    lla = eci2lla(rvec' * 1e3, datevec(t_datetime));

    if include_Oxygen
        otype = "Oxygen";
    else
        otype = "NoOxygen";
    end

    if lla(3) > 80000
        [~, rho_outputs] = atmosnrlmsise00(lla(3), lla(1), lla(2), t_datetime.Year, day(t_datetime, 'dayofyear'), second(t_datetime, 'secondofday'), otype, 'None'); % (up to 1000 km ish)
        rho_kg_m3 = rho_outputs(6);
    else
        [~, a, ~, rho_kg_m3] = atmosisa(lla(3),extended=true);
        M = v * 1000 / a;
        fprintf("Mach: %.3f\n",M);
    end

    rho_kg_km3 = rho_kg_m3 * (1e3)^3;
    a_drag = -1 / 2 * rho_kg_km3 * v * vvec * A_over_m * C_D;
    % 
    % % Convert to RTN
    % x_keplerian = cartesian_to_keplerian(x_cartesian, [0; 0; 1], [1; 0; 0], mu);
    % e = x_keplerian(2);
    % i = x_keplerian(3);
    % Omega = x_keplerian(4);
    % omega = x_keplerian(5);
    % M = x_keplerian(6);
    % 
    % nu = eccentric_to_true_anomaly(mean_to_eccentric_anomaly(M, e), e);
    % 
    % a_drag_RTN = cartesian_to_RTN_DCM(i, Omega, omega, nu)' * a_drag;
end

function [altitude, isterminal, direction] = altitude_event_function(t, x_cartesian, R_E, stop_altitude, t0_datetime)
    % Stops simulation when specific altitude is hit

    t_datetime = t0_datetime + seconds(t);

    lla = eci2lla(x_cartesian(1:3)' * 1e3, datevec(t_datetime));
    
    fprintf("Time (hr): %.3f, Alt (km): %.3f\n", t / 60 / 60, lla(3) / 1000)
    altitude = lla(3) / 1000 - stop_altitude;
    isterminal = 1; % Stop simulation
    direction = -1; % Trigger when zero is approached from positive side
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
