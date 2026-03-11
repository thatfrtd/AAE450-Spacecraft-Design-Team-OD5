function [v_rot] = quat_rot(q, v)
%QUAT_ROT Summary of this function goes here
%   Detailed explanation goes here

v_rot = eye([3, 4]) * q_mul(q, q_mul([v; 0], q_conj(q)));

end

