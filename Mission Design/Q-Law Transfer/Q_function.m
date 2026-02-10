function [Q] = Q_function(oe, oe_t, W_p, P, S_oe, W_oe, oedot_xx)
%Q_FUNCTION Compute Q-function to determine distance to desired orbit
%   The Q-function uses the slow elements of modified equinoctial elements
%   - HOWEVER - it uses semimajor axis, a, instead of semilatus rectum, p,
%   because the authors found it to yield better control
arguments
    oe % Current orbital elements (a, f, g, h, k)
    oe_t % Desired orbital elements (a, f, g, h, k)
    W_p % Penalty weight term
    P % Penalty violation (positive is violation)
    S_oe % Scaling factors
    W_oe % Weighting factors for orbital elements
    oedot_xx % Max rate of change of elements over the thrust direction and 
    % true anomaly on the osculating orbit (current orbit)
end

penalty_multiplier = (1 + W_p * P);

% Q describes the distance between current and desired orbits
Q = penalty_multiplier * sum(S_oe .* W_oe .* ((oe - oe_t) ./ oedot_xx) .^ 2);

end