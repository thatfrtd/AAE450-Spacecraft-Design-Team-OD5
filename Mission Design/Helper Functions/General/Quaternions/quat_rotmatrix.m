function [R] = quat_rotmatrix(q)
    w = q(4);
    v = q(1:3);

    R = (w ^ 2 - v' * v) * eye(3) + 2 * v * v' + 2 * w * skew(v);
end