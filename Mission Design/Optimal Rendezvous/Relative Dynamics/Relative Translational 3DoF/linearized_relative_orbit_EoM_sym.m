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
M = sqrt(mu / a_c ^ 3) * t + M0_c;
% Three term bessel expansion of Kepler's Equation so we can take its derivative
E = M + 2 * besselj(1, e_c) * sin(M) + besselj(2, 2 * e_c) * sin(2 * M);
r_c = a_c * (1 - e_c * cos(E));
v_c = sqrt(mu * (2 / r_c - 1 / a_c));
r_cdot = sin(acos(h / (r_c * v_c))) * v_c;
thetastardot = h / r_c ^ 2;
thetastarddot = -2 * h  / r_c ^ 3 * r_cdot;

%% Linearized Relative Dynamics
rdot = v;
vdot = [thetastardot ^ 2 + 2 * mu / r_c ^ 3, thetastarddot, 0; 
       -thetastarddot, thetastardot ^ 2 - mu / r_c ^ 3, 0; 
        0, 0, -mu / r_c ^ 3] * r ...
     + [0, 2 * thetastardot, 0; -2 * thetastardot, 0, 0; 0, 0, 0] * v ...
     + u / m;
mdot = -alpha * sqrt(u(1) ^ 2 + u(2) ^ 2 + u(3) ^ 2);

xdot = [rdot; vdot; mdot];

%% Create Equations of Motion Function for Optimizer
matlabFunction(xdot,"File","linearized_relative_orbit_EoM","Vars", [{t}; {x}; {u}; {p}; {c}]);