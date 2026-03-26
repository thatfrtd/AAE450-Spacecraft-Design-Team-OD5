function [x_norm] = dnorm(x)
%DNORM Summary of this function goes here
%   Detailed explanation goes here

x_norm = sqrt(sum(abs(x)));

end