function [x2_cart_ck, x1_kep, M2_ck] = propagate_conic_array_R(x1_cart, ToF, mu, options)
arguments
    x1_cart
    ToF
    mu
    options.n_pts = 200;
end

nT = numel(ToF);
x2_cart_ck = zeros(6, nT);
x1_kep = zeros(6, nT);
M2_ck = zeros(1, nT);

for i = 1:nT
    [x2_cart_ck(:,i), x1_kep(:,i), M2_ck(i)] = propagate_conic_R( ...
        x1_cart, ToF(i), mu, n_pts = options.n_pts);
end
end
