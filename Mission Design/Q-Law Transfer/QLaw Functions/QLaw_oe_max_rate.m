function [oedot_xx] = QLaw_oe_max_rate(oe, mu, F)
%QLAW_OE_MAX_RATE Calculate max rate of change of orbital elements over current orbit
%   Calculates the maximum rate of change of the slow modified equinoctial 
% elements over the current orbit. The f and g element formulas are
% approximate.
arguments
    oe % Current orbital elements (a, f, g, h, k)
    mu % Gravitational parameter
    F % Perturbing force magnitude
end

% Read in orbital elements
a = oe(1);
f = oe(2);
g = oe(3);
h = oe(4);
k = oe(5);

p = a * (1 - f ^ 2 - g ^ 2);
e = sqrt(f ^ 2 + g ^ 2);
s_sqr = 1 + h ^ 2 + k ^ 2;

% Calculate 
adot_xx = 2 * F * a * sqrt(a / mu) * sqrt((1 + e) / (1 - e));
fdot_xx = 2 * F * sqrt(p / mu); % APPROXIMATE
gdot_xx = fdot_xx; % APPROXIMATE
hdot_xx = 1 / 2 * F * sqrt(p / mu) * s_sqr / (sqrt(1 - g ^ 2) + f);
kdot_xx = 1 / 2 * F * sqrt(p / mu) * s_sqr / (sqrt(1 - f ^ 2) + g);

% Package outputs
oedot_xx = [adot_xx; fdot_xx; gdot_xx; hdot_xx; kdot_xx];

end