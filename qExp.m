function [q] = qExp(tau)
    theta = norm(tau);
    u = tau / theta;
    q = [u * sin(theta / 2); cos(theta / 2)];
end