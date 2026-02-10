function [f_0] = f0_modified_equinoctial_with_a(x, mu)
    a = x(1, :);
    f = x(2, :);
    g = x(3, :);
    L = x(6, :);

    q = 1 + f .* cos(L) + g .* sin(L);
    p = a * (1 - f ^ 2 - g ^ 2);

    f_0 = [zeros(5, 1); ...
           sqrt(mu * p) .* (q ./ p) .^ 2];
end