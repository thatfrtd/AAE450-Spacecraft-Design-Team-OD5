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

char_star = load_charecteristic_values_Earth();

% Constants
x_keplerian_c = sym("x_keplerian_c", [6, 1]);
I = sym("I", [3; 1]);
T_1_direction = sym("T_1_direction", [3, 1]);
a_c = x_keplerian_c(1);
e_c = x_keplerian_c(2);
i_c = x_keplerian_c(3);
Omega_c = x_keplerian_c(4);
omega_c = x_keplerian_c(5);
M0_c = x_keplerian_c(6);
Isp = sym("Isp", [2, 1]);
g_0 = sym("g_0", [1, 1]);
mu = sym("mu", [1, 1]);
R_C_S = sym("R_C_S", [6, 8]); % RCS thrust and torque allocation matrix
c = [x_keplerian_c; Isp; mu; g_0; I; T_1_direction; R_C_S(:)];

alpha_1 = 1 / (Isp(1) * g_0);
alpha_2 = 1 / (Isp(2) * g_0);

% states

r = sym("r", [3;1]);
v = sym("v", [3,1]);
q = sym("theta", [4;1]);
w = sym("w", [3;1]);
m_nd = sym("m", [1;1]);

x = [r; v; q; w; m_nd];

% controls

T_1_mag = sym("T_1_mag", [1; 1]); % Main thruster
T_RCS = sym("T_RCS", [8, 1]); % RCS thruster forces
tau = sym("tau", [3; 1]); % Reaction wheels

u = [T_1_mag; T_RCS; tau];

thrust_mag_2 = norms(T_RCS, 1);

RCS_net = R_C_S * T_RCS; % [thrust; torque] vectors

% Parameters
p = sym("p", [0, 1]);

% calculate qdot
qdot = 1 / 2 * q_mul(q, [w; 0]);

% r dot
r_dot = v;

% v dot (nonlinear relative orbit dynamics)
T_1_I = quat_rot(q, T_1_direction * T_1_mag);
T_2_I = quat_rot(q, RCS_net(1:3));
rel_orbit_xdot = nonlinear_relative_orbit_EoM_twothruster_dim(t, x([1:6, 14].'), [T_1_I; T_2_I], p, c(1:10));
v_dot = rel_orbit_xdot(4:6);
%
M = RCS_net(4:6) + tau;
% w dot
w_dot = (M([1; 2; 3]) + (I([2; 3; 1]) - I([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ I([1; 2; 3]);

% mdot
m_dot = (-alpha_1 * T_1_mag - alpha_2 * thrust_mag_2) / char_star.m;

x_dot = [r_dot; v_dot; qdot; w_dot; m_dot];

vars = [{t}; {x}; {u}; {p}; {c}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","relative_orbit_6DoF_twothruster_EoM_noxRCS","Vars", [{t}; {x}; {u}; {p}; {c}]);