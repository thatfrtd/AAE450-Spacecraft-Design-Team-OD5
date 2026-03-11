%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Symbolic Rotational 3DoF EoM
% Author: Travis Hastreiter 
% Created On: 1 March, 2026
% Description: 3DoF rotational 3DoF dynamics with reaction wheel (or other 
% propellant less torquer) torque
% Most Recent Change: 1 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = sym("t");

% Constants
I = sym("I", [3; 1]);
c = [I];

% states

q = sym("theta", [4;1]);
w = sym("w", [3;1]);

x = [q; w];

% controls

tau = sym("tau", [3; 1]); % Reaction wheels

u = [tau];

% Parameters
p = sym("p", [0, 1]);

% calculate qdot
qdot = 1 / 2 * q_mul(q, [w; 0]);

% w dot
M = tau; % should have RCS thrusters too

w_dot = (M([1; 2; 3]) + (I([2; 3; 1]) - I([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ I([1; 2; 3;]);

x_dot = [qdot; w_dot];

vars = [{t}; {x}; {u}; {p}; {c}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","rotational_3DoF_EoM","Vars", [{t}; {x}; {u}; {p}; {c}]);