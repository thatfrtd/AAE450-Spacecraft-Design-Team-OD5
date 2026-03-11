function [q_prod] = q_mul_array(q_A, q_B)
%Q_MUL Summary of this function goes here
%   Detailed explanation goes here

N = size(q_A, 2);

q_prod = zeros(4, N);

for i = 1 : N
    q_prod(:, i) = q_mul(q_A(:, i), q_B(:, i));
end

end

