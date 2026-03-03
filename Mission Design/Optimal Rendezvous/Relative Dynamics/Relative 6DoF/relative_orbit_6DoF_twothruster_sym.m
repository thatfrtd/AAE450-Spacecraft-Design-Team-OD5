%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Symbolic Nonlinear 6DoF Relative Orbit EoM
% Author: Travis Hastreiter 
% Created On: 1 March, 2026
% Description: 6DoF relative orbital dynamics with changing mass and a main
% thruster with a specific Isp and RCS thrusters with different Isp along
% with reaction wheel (or other propellant less torquer) torque
% Most Recent Change: 1 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = sym("t");

g_0 = 9.81e-3; % [km / s2]

char_star = load_charecteristic_values_Earth();

mu = 1; % Nondimensionalized

% Constants
x_keplerian_c = sym("x_keplerian_c", [6, 1]);
I = sym("I", [3; 1]);
T_1_direction = sym("T_1_direction", [3, 1]);
a_c = x_keplerian_c(1) / char_star.l;
e_c = x_keplerian_c(2);
i_c = x_keplerian_c(3);
Omega_c = x_keplerian_c(4);
omega_c = x_keplerian_c(5);
M0_c = x_keplerian_c(6);
Isp = sym("Isp", [2, 1]);
c = [x_keplerian_c; Isp; I; T_1_direction];

alpha_1 = char_star.v / (Isp(1) * g_0);
alpha_2 = char_star.v / (Isp(2) * g_0);

% states

r = sym("r", [3;1]);
v = sym("v", [3,1]);
q = sym("theta", [4;1]);
w = sym("w", [3;1]);
m = sym("m", [1;1]);

x = [r; v; q; w; m];

% controls

T_1_mag = sym("T_1_mag", [1; 1]); % Main thruster
T_2 = sym("T_2", [3; 1]); % Should also provide torque... need allocation matrix
tau = sym("tau", [3; 1]); % Reaction wheels

u = [T_1_mag; T_2; tau];

thrust_mag_2 = sqrt(T_2(1) ^ 2 + T_2(2) ^ 2 + T_2(3) ^ 2);

% Parameters
p = sym("p", [0, 1]);

% calculate qdot
qdot = 1 / 2 * q_mul(q, [w; 0]);

% r dot
r_dot = v;

% v dot (nonlinear relative orbit dynamics)
T_1_I = quat_rot(q, T_1_direction * T_1_mag);
T_2_I = quat_rot(q, T_2);
rel_orbit_xdot = nonlinear_relative_orbit_EoM_twothruster(t, x([1:6, 14]'), [T_1_I; T_2_I], p, c(1:8));
v_dot = rel_orbit_xdot(4:6);
%
M = tau; % should have RCS thrusters too
% w dot
w_dot = (M([1; 2; 3]) + (I([2; 3; 1]) - I([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ I([1; 2; 3;]);

% mdot
m_dot = -alpha_1 * T_1_mag - alpha_2 * thrust_mag_2;

x_dot = [r_dot; v_dot; qdot; w_dot; m_dot];

vars = [{t}; {x}; {u}; {p}; {c}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","relative_orbit_6DoF_twothruster_EoM","Vars", [{t}; {x}; {u}; {p}; {c}]);