%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Symbolic Nonlinear Relative Orbit EoM Dimensionalized
% Author: Travis Hastreiter 
% Created On: 23 March, 2026
% Description: Nonlinear relative orbital motion dynamics with two
% thrusters with different Isp (and max thrust)
% Most Recent Change: 23 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

char_star = load_charecteristic_values_Earth();

% Constants
mu = sym("mu", [1, 1]);
g_0 = sym("g_0", [1, 1]);
x_keplerian_c = sym("x_keplerian_c", [6, 1]);
a_c = x_keplerian_c(1);
e_c = x_keplerian_c(2);
i_c = x_keplerian_c(3);
Omega_c = x_keplerian_c(4);
omega_c = x_keplerian_c(5);
M0_c = x_keplerian_c(6);
Isp = sym("Isp", [2, 1]);
c = [x_keplerian_c; Isp; mu; g_0];

alpha_1 = 1 / (Isp(1) * g_0);
alpha_2 = 1 / (Isp(2) * g_0);

% State
t = sym("t");
r = sym("r", [3, 1]);
v = sym("v", [3, 1]);
m_nd = sym("m", [1, 1]);
x = [r; v; m_nd];

u = sym("thrust", [6, 1]); % Thrust
u_1 = u(1:3);
u_2 = u(4:6);
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

%% Nonlinear Relative Dynamics
r_d = norm([r_c; 0; 0] + r);

m = m_nd * char_star.m;

rdot = v;
vdot = [2 * thetastardot * (v(2) - r(2) * r_cdot / r_c) + r(1) * thetastardot ^ 2 + mu / r_c ^ 2 - mu / r_d ^ 3 * (r_c + r(1));
       -2 * thetastardot * (v(1) - r(1) * r_cdot / r_c) + r(2) * thetastardot ^ 2 - mu / r_d ^ 3 * r(2);
       -mu / r_d ^ 3 * r(3)] ...
     + (u_1 + u_2) / 1000 / m;
mdot = (-(alpha_1 * sqrt(u_1(1) ^ 2 + u_1(2) ^ 2 + u_1(3) ^ 2) ...
        + alpha_2 * sqrt(u_2(1) ^ 2 + u_2(2) ^ 2 + u_2(3) ^ 2))) / char_star.m;

xdot = [rdot; vdot; mdot];

%% Create Equations of Motion Function for Optimizer
matlabFunction(xdot,"File","nonlinear_relative_orbit_EoM_twothruster_dim","Vars", [{t}; {x}; {u}; {p}; {c}]);