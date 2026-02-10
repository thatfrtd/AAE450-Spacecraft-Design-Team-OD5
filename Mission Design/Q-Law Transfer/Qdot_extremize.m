function [Qdot_min, Qdot_max] = Qdot_extremize(oe, L, partial_Q_partial_oe, Qdot_opt_params)
%QDOT_EXTREMIZE Summary of this function goes here
%   Detailed explanation goes here
arguments
    oe
    L
    partial_Q_partial_oe
    Qdot_opt_params
end

% Know min w.r.t. alpha, beta of Qdot in terms of L
min_ab_Qdot = @(L) -norm(partial_Q_partial_oe * B_oe_func(oe, L));

% Initialize search for L


% Min w.r.t. L (min w.r.t. alpha, beta (Qdot))
Qdot_min = ;


% Max w.r.t. L (min w.r.t. alpha, beta (Qdot))
Qdot_max = ;

end