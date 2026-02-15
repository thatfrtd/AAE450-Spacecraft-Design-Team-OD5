function [B] = B_modified_equinoctial(x, mu)
    p = x(1);
    f = x(2);
    g = x(3);
    h = x(4);
    k = x(5);
    L = x(6);

    q = 1 + f * cos(L) + g * sin(L);
    s_sqr = 1 + h ^ 2 + k ^ 2;

    cons_1 = (h * sin(L) - k * cos(L)) / q;
    B = sqrt(p / mu) * [0, 2 * p / q, 0;
        sin(L), ((q + 1) * cos(L) + f) / q, -g * cons_1;
        -cos(L), ((q + 1) * sin(L) + g) / q, f * cons_1;
        0, 0, s_sqr / (2 * q) * cos(L);
        0, 0, s_sqr / (2 * q) * sin(L);
        0, 0, cons_1];
    % 
    % % Make B convert disturbance into RTN frame from cartesian frame
    % i = atan2(2 * sqrt(h ^ 2 + k ^ 2), 1 - h ^ 2 - k ^ 2);
    % Omega = atan2(k, h);
    % omega = atan2(g * h - f * k, f * h + g * k);
    % nu = L - (Omega + omega);
    % B = B * cartesian_to_RTN_DCM(i, Omega, omega, nu)';
end