function [q_mm_star] = q_conj_mul_matrix(q)
%Q_CONJ_MUL_MATRIX Summary of this function goes here
%   Detailed explanation goes here

q_mm_star = [q(4) * eye(3) + skew(q(1:3)), q(1:3); -q(1:3)', q(4)];

end

