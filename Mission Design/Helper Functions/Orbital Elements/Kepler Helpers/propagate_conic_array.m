function [x2_cart_ck, x1_kep, M2_ck] = propagate_conic_array(x1_cart, ToF, mu, options)
%PROPAGATE_CONIC Summary of this function goes here
%   Input state x_1 has to be in an inertial reference frame
arguments
    x1_cart 
    ToF 
    mu 
    options.lagrange_tolerance = 1e-8
    options.one_ToF_to_one_state = false
end

x2_cart_ck = zeros([6, numel(ToF)]);
x1_kep = zeros([6, numel(ToF)]);
M2_ck = zeros([1, numel(ToF)]);

if ~options.one_ToF_to_one_state
    for i = 1 : numel(ToF)
        [x2_cart_ck(:, i), x1_kep(:, i), M2_ck(i)] = propagate_conic(x1_cart, ToF(i), mu, lagrange_tolerance = options.lagrange_tolerance);
    end
else
    for i = 1 : numel(ToF)
        [x2_cart_ck(:, i), x1_kep(:, i), M2_ck(i)] = propagate_conic(x1_cart(:, i), ToF(i), mu, lagrange_tolerance = options.lagrange_tolerance);
    end
end

end