function [x_norm] = dnorm(x)
%DNORM Summary of this function goes here
%   Detailed explanation goes here

x_norm = sqrt(sum(x .^ 2));

end