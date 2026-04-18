function [Qdot_min, Qdot_max] = Qdot_extremize_fast(oe, partial_Q_partial_oe, num_start_points)
%QDOT_EXTREMIZE Summary of this function goes here
%   The minimum and maximum w.r.t. L of the minimum Qdot w.r.t. thrust 
% direction needs to be found so that the current L can be compared to the
% best and worst Ls to compute effeciences used to determine if coasting
% should be done or not. There is usually two minimums and maximums so the
% right ones should try to be chosen. Care must be taken to make this
% process fast.
arguments
    oe
    partial_Q_partial_oe
    num_start_points
end

% Know min w.r.t. alpha, beta of Qdot in terms of L
min_ab_Qdot = @(L) -norm(partial_Q_partial_oe * B_oe_func(oe, L));

% Initialize search for L
L_start = linspace(0, 2 * pi, num_start_points + 1);
Qdot_start = zeros([num_start_points, 1]);
for s = 1 : num_start_points
    Qdot_start(s) = min_ab_Qdot(L_start(s));
end

% Just look at the min and max of start point array, no iteration

% Min w.r.t. L (min w.r.t. alpha, beta (Qdot))
Qdot_min = min(Qdot_start);
    
% Max w.r.t. L (min w.r.t. alpha, beta (Qdot))
Qdot_max = max(Qdot_start);

end