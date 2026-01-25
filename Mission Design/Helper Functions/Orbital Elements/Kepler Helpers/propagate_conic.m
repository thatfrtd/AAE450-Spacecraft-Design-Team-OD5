function [x2_cart_ck, x1_kep, M2_ck] = propagate_conic(x1_cart, ToF, mu, options)
%PROPAGATE_CONIC Summary of this function goes here
%   Input state x_1 has to be in an inertial reference frame
arguments
    x1_cart 
    ToF 
    mu 
    options.lagrange_tolerance = 1e-8;
end

% Extract initial position and velocity
rvec_1 = x1_cart(1:3);
vvec_1 = x1_cart(4:6);

% Convert state to Keplerian elements
x1_kep = cartesian_to_keplerian(x1_cart, [0; 0; 1], [1; 0; 0], mu);

% Calculate Lagrange coefficients f, g, fdot, gdot
[f, g, fdot, gdot, M2_ck] = lagrange_coefficients(x1_kep, ToF, mu);

% Check coefficients
lagrange_ck = abs(f * gdot - fdot * g - 1); % Should be zero
if lagrange_ck > options.lagrange_tolerance
    warning("Lagrange coefficient accuracy %.3g is more than %.3g", lagrange_ck, options.lagrange_tolerance);
end

% Find final state
rvec_2 = f * rvec_1 + g * vvec_1;
vvec_2 = fdot * rvec_1 + gdot * vvec_1;

% Return as a 6x1 column vector [r; v]
x2_cart_ck = [rvec_2; vvec_2];

end