%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Q-Law Orbit Transfer Example 
% Author: Travis Hastreiter 
% Created On: 9 February, 2026
% Description: Calculate partial derivatives needed to calculate D1, D2, D3,
% for Q-Law's coasting heuristic
% Most Recent Change: 9 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constants
mu = 1; % Assume everything is nondimensionalized

%% Inputs
% Target orbit
oe_t = sym("oe_t", [5, 1]);

% Q-Law parameters
W_oe = sym("W_oe", [5, 1]);
m = sym("m", [1, 1]);
n = sym("n", [1, 1]);
r = sym("r", [1, 1]);

% Spacecraft parameters
F_max = sym("F_max", [1, 1]);

% Penalty
W_p = sym("W_p", [1, 1]);
r_p_min = sym("r_p_min", [1, 1]); % Nondimensionalized
k_p = sym("k", [1, 1]);

%% Q Function Calculation
% Symbolic
oe = sym("oe", [5, 1]);
a = oe(1);
f = oe(2);
g = oe(3);
h = oe(4);
k = oe(5);
L = sym("L", [1, 1]); % True longitude
x_me = [oe; L];

% Intermediate
e = sqrt(f ^ 2 + g ^ 2);
r_p = a * (1 - e);

% Penalty
P = periapsis_penalty(r_p, r_p_min, k_p);

matlabFunction(P, "File","periapsis_penalty_sym", "Vars", [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");

% Q function
S_oe = QLaw_scaling(a, oe_t(1), m, n, r);

matlabFunction(S_oe, "File","QLaw_scaling_sym", "Vars", [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");

[oedot_xx] = QLaw_oe_max_rate(oe, mu, F_max);

matlabFunction(oedot_xx, "File","QLaw_oe_max_rate_sym", "Vars", [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");

Q = Q_function(oe, oe_t, W_p, P, S_oe, W_oe, oedot_xx);

matlabFunction(Q, "File","Q_func_sym", "Vars", [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");

%% Partial Derivatives
partial_Q_partial_oe = jacobian(Q, oe);
matlabFunction(partial_Q_partial_oe, "File","partial_Q_partial_oe_func", "Vars", [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");

% Gauss Planetary Equations
B = B_modified_equinoctial_with_a(x_me, mu);
B_oe = B(1:5, :); % Partial oedot w.r.t. partial F

matlabFunction(B_oe, "File","B_oe_func", "Vars", [{oe}; {L}],"Optimize",true, "Comments","Inputs: [{oe}; {L}]");

%% DRIZZLING DRACONIAN DAREDEVILS
D = partial_Q_partial_oe * B_oe;

D_norm = sqrt(D(1) ^ 2 + D(2) ^ 2 + D(3) ^ 2);

matlabFunction(D,"File","D_func", "Vars", [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");
matlabFunction(D_norm,"File","D_norm_func", "Vars", [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");

%% Create Function for Qdot Reusing D Computations
[a_control, alpha, beta] = QLaw_thrust_mapping(D, 1);
xdot_me = gauss_planetary_eqn(f0_modified_equinoctial_with_a(x_me, mu), B, a_control);
odot = xdot_me(1:5);

Qdot = partial_Q_partial_oe * odot;

matlabFunction(Qdot,"File","Qdot_func", "Vars", [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");

% Combined D, Q, Qdot, P, and partial_Q_partial_oe function
matlabFunction(D, Q, Qdot, P, partial_Q_partial_oe, "File", "D_Q_Qdot_P_partial_Q_partial_oe_func", "Vars", [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}],"Optimize",true, "Comments","Inputs: [{oe}; {L}; {oe_t}; {W_oe}; {m}; {n}; {r}; {F_max}; {W_p}; {r_p_min}; {k_p}]");