function [qq_prod] = qq_mul(qq_A, qq_B)
%QQ_MUL Summary of this function goes here
%   Detailed explanation goes here

qq_mm_A = qq_mul_matrix(qq_A);
qq_prod = qq_mm_A * qq_B;

end

