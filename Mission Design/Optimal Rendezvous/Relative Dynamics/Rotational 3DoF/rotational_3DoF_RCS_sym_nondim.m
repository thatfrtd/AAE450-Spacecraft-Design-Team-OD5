%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Symbolic Rotational 3DoF EoM with RCS Thrusters Nondimensionalized
% Author: Travis Hastreiter 
% Created On: 21 March, 2026
% Description: 3DoF rotational 3DoF dynamics with reaction wheel (or other 
% propellant less torquer) torque and RCS thrusters where quantities are
% nondimensionalized
% Most Recent Change: 21 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = sym("t");
char_star = load_charecteristic_values_Earth();

% Constants
I = sym("I", [3; 1]);
I_nd = I / char_star.m / char_star.l ^ 2;
R_C_S = sym("R_C_S", [6, 12]); % RCS thrust and torque allocation matrix
c = [I; R_C_S(:)];

% states
q = sym("theta", [4;1]);
w = sym("w", [3;1]);

x = [q; w];

% controls
T_RCS = sym("T_RCS", [12, 1]); % RCS thruster forces
tau = sym("tau", [3; 1]); % Reaction wheels

u = [T_RCS; tau]; 
T_RCS_nd = T_RCS / 1000 / char_star.F / char_star.l;
tau_nd = tau / 1000 / char_star.F / char_star.l;

% Parameters
p = sym("p", [0, 1]);

RCS_net = R_C_S * T_RCS_nd; % [thrust; torque] vectors

% calculate qdot
qdot = 1 / 2 * q_mul(q, [w; 0]);

% w dot
M = (RCS_net(4:6) + tau_nd);

w_dot = (M([1; 2; 3]) + (I_nd([2; 3; 1]) - I_nd([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ I_nd([1; 2; 3]);

x_dot = [qdot; w_dot];

vars = [{t}; {x}; {u}; {p}; {c}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","rotational_3DoF_RCS_EoM_nondim","Vars", [{t}; {x}; {u}; {p}; {c}]);