function [q_mm] = q_mul_matrix(q)
%Q_MUL_MATRIX Summary of this function goes here
%   Detailed explanation goes here

q_mm = [q(4) * eye(3) + skew(q(1:3)), q(1:3); -q(1:3)', q(4)];

end

