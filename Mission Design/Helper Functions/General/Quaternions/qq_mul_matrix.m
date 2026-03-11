function [qq_mm] = qq_mul_matrix(qq)
%QQ_MUL_MATRIX Summary of this function goes here
%   Detailed explanation goes here

q1_mm = q_mul_matrix(qq(1:4));
q2_mm = q_mul_matrix(qq(5:8));

qq_mm = [q1_mm, zeros([4, 4]); q2_mm, q1_mm];

end

