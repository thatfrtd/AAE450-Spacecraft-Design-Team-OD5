function [B] = B_modified_equinoctial_with_a(x, mu)
    a = x(1);
    f = x(2);
    g = x(3);
    h = x(4);
    k = x(5);
    L = x(6);

    q = 1 + f * cos(L) + g * sin(L);
    s_sqr = 1 + h ^ 2 + k ^ 2;

    p = (1 - f ^ 2 - g ^ 2) * a;
    r = p / q;
    e = sqrt(f ^ 2 + g ^ 2);
    nu = L - atan(g / f);

    cons_1 = (h * sin(L) - k * cos(L)) / q;
    B = sqrt(p / mu) * ...
        [2 * a ^ 2 * e * sin(nu) / p, 2 * a ^ 2 / r, 0;
        sin(L), ((q + 1) * cos(L) + f) / q, -g * cons_1;
        -cos(L), ((q + 1) * sin(L) + g) / q, f * cons_1;
        0, 0, s_sqr / (2 * q) * cos(L);
        0, 0, s_sqr / (2 * q) * sin(L);
        0, 0, cons_1];
end