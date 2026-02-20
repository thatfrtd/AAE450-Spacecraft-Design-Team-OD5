function [eta_a, eta_r] = QLaw_efficiencies(Qdot, Qdot_min, Qdot_max)
%QLAW_EFFICIENCIES Summary of this function goes here
%   Used as a heuristic to determine if thrusting should happen. If above a
%   threshold then thrusting happens
arguments
    Qdot
    Qdot_min
    Qdot_max
end

% Absolute efficiency
eta_a = Qdot / Qdot_min;

% Relative efficiency
eta_r = (Qdot - Qdot_max) / (Qdot_min - Qdot_max);

end