function [f_0] = f0_modified_equinoctial(x, mu)
    p = x(1, :);
    f = x(2, :);
    g = x(3, :);
    L = x(6, :);

    q = 1 + f .* cos(L) + g .* sin(L);

    f_0 = [zeros(5, 1); ...
           sqrt(mu * p) .* (q ./ p) .^ 2];
end