function [P] = periapsis_penalty(r_p, r_p_min, k)
%PERIAPSIS_PENALTY Exponential soft penalty function for min periapsis
%   Detailed explanation goes here
arguments
    r_p
    r_p_min
    k
end

P = exp(k * (1 - r_p / r_p_min));

end