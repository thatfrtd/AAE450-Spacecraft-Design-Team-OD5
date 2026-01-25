function [r] = eccentric_radius(a, e, E)
%ECCENTRIC_RADIUS Summary of this function goes here
%   Detailed explanation goes here

r = a * (1 - e * cos(E));

end