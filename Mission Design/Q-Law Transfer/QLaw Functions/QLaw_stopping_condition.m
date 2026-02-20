function [Q_stop] = QLaw_stopping_condition(R_c, W_oe)
%QLAW_STOPPING_CONDITION Summary of this function goes here
%   If the Q-function is less than the stopping condition then the transfer
%   is within the tolerances and is deemed to be converged
arguments
    R_c % Parameter controlling error of the remaining distance from target orbit
    W_oe % Scaling factors for orbital elements
end

Q_stop = R_c * sqrt(sum(W_oe));

end