function [DCM] = RTN_to_ECI(r, v)
%RTN_TO_ECI Summary of this function goes here
%   Detailed explanation goes here

R_hat = r / sqrt(r(1) ^ 2 + r(2) ^ 2 + r(3) ^ 2);
N = cross(r, v);
N_hat = N / sqrt(N(1) ^ 2 + N(2) ^ 2 + N(3) ^ 2);
T_hat = cross(N_hat, R_hat);

DCM = [R_hat, T_hat, N_hat];

end