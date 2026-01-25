function [r] = hyperbolic_radius(a, e, H)
%HYPERBOLIC_RADIUS Summary of this function goes here
%   Detailed explanation goes here

r = a * (1 - e * cosh(H));

end