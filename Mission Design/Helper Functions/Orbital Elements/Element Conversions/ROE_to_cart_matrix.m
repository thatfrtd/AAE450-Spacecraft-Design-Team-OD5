function [T] = ROE_to_cart_matrix(n, t)
    nt = n * t;
    T = [1, 0, -cos(nt), -sin(n * t), 0, 0;
         0, 1, 2 * sin(nt), -2 * cos(nt), 0, 0;
         0, 0, 0, 0, sin(nt), -cos(nt);
         0, 0, n * sin(nt), -n * cos(nt), 0, 0;
         -3 * n / 2, 0, 2 * n * cos(nt), 2 * n * sin(nt), 0, 0;
         0, 0, 0, 0, n * cos(nt), n * sin(nt)];
end
