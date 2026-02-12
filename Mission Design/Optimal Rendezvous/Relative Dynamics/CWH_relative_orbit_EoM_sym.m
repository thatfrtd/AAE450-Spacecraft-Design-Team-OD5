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
r_c = sym("r_c", [1, 1]);
Isp = sym("Isp", [1, 1]);
c = [r_c; Isp];

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

%% CWH Relative Dynamics
n = sqrt(mu / (r_c / char_star.l) ^ 3);

rdot = v;
vdot = [3 * n ^ 2, 0, 0; 0, 0, 0; 0, 0, -n ^ 2]  * r ...
    +  [0, 2 * n, 0; -2 * n, 0, 0; 0, 0, 0] * v ...
    + u / m;

mdot = -alpha * sqrt(thrust(1) ^ 2 + thrust(2) ^ 2 + thrust(3) ^ 2);

xdot = [rdot; vdot; mdot];

%% Create Equations of Motion Function for Optimizer
matlabFunction(xdot,"File","CWH_relative_orbit_EoM","Vars", [{t}; {x}; {u}; {p}; {c}]);