function [qq] = construct_dual_quat_array(q, r_I)
%CONSTRUCT_DUAL_QUAT Summary of this function goes here
%   Detailed explanation goes here

N = size(q, 2);

qq = zeros(8, N);

for i = 1 : N
    qq(:, i) = construct_dual_quat(q(:, i), r_I(:, i));
end

end

