function [q_prod] = q_mul(q_A, q_B)
%Q_MUL Summary of this function goes here
%   Detailed explanation goes here

q_mm_A =  q_mul_matrix(q_A);
q_prod = q_mm_A * q_B;

end

