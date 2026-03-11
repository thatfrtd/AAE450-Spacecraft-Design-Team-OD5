function [qq] = construct_dual_quat(q, r_I)
%CONSTRUCT_DUAL_QUAT Summary of this function goes here
%   Detailed explanation goes here

qq = [q; 1/2 * q_mul([r_I; 0], q)];

end

