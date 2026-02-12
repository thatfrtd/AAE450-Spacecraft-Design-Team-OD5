%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Symbolic Nonlinear Relative Orbit EoM
% Author: Travis Hastreiter 
% Created On: 11 February, 2026
% Description: 3DoF rocket landing dynamics with changing mass
% Most Recent Change: 11 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g_0 = 9.81e-3; % [km / s2]

char_star = load_charecteristic_values_Earth();

mu = 1; % Nondimensionalized

% Constants
x_keplerian_c = sym("x_keplerian_c", [6, 1]);
a_c = x_keplerian_c(1) / char_star.l;
e_c = x_keplerian_c(2);
i_c = x_keplerian_c(3);
Omega_c = x_keplerian_c(4);
omega_c = x_keplerian_c(5);
M0_c = x_keplerian_c(6);
Isp = sym("Isp", [1, 1]);
c = [x_keplerian_c; Isp];

alpha = char_star.v / (Isp * g_0);

% State
t = sym("t");
r = sym("r", [3, 1]);
v = sym("v", [3, 1]);
m = sym("m", [1, 1]);
x = [r; v; m];

thrust = sym("thrust", [3, 1]); % Thrust (nondimensionalized by kN)
u = [thrust];
p = sym("p", [0, 1]);

%% Calculate Chief Orbit Properties
p_c = a_c * (1 + e_c) * (1 - e_c);
h = sqrt(p_c * mu);
M = @(t) sqrt(mu / a_c ^ 3) * (t * char_star.t) + M0_c;
% Three term bessel expansion of Kepler's Equation so we can take its derivative
E = @(t) M(t) + 2 * besselj(1, e_c) * sin(M(t)) + besselj(2, 2 * e_c) * sin(2 * M(t));
r_c = @(t) a_c * (1 - e_c * cos(E(t)));
v_c = @(t) sqrt(mu * (2 / r_c(t) - 1 / a));
r_cdot = @(t) sin(acos(h / (r_c(t) * v_c(t)))) * v_c(t);
thetastardot = @(t) h / r_c(t) ^ 2;

%% Nonlinear Relative Dynamics
r_d = norm([r_c(t); 0; 0] + r);

rdot = v;
vdot = [2 * thetastardot(t) * (v(2) - r(2) * r_cdot / r_c(t)) + r(1) * thetastardot(t) ^ 2 + mu / r_c(t) ^ 2 - mu / r_d ^ 3 * (r_c(t) + r(1));
       -2 * thetastardot(t) * (v(1) - r(1) * r_cdot / r_c(t)) + r(2) * thetastardot(t) ^ 2 - mu / r_d ^ 3 * r(2);
       -mu / r_d ^ 3 * r(3)] ...
     + thrust / m;
mdot = -alpha * sqrt(thrust(1) ^ 2 + thrust(2) ^ 2 + thrust(3) ^ 2);

xdot = [rdot; vdot; mdot];

%% Create Equations of Motion Function for Optimizer
matlabFunction(xdot,"File","nonlinear_relative_orbit_EoM","Vars", [{t}; {x}; {u}; {p}; {c}]);