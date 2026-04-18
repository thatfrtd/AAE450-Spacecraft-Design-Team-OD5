function [T_inv] = cart_to_ROE_matrix(n, t)
    nt = n * t;
    T_inv = [4, 0, 0, 0, 2 / n, 0;
             0, 1, 0, -2 / n, 0, 0;
             3 * cos(nt), 0, 0, sin(nt) / n, 2 * cos(nt) / n, 0;
             3 * sin(nt), 0, 0, -cos(nt) / n, 2 * sin(nt) / n, 0;
             0, 0, sin(nt), 0, 0, cos(nt) / n;
             0, 0, -cos(nt), 0, 0, sin(nt) / n];
end