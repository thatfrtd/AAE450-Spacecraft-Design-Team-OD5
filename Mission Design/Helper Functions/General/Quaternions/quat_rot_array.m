function [v_rot] = quat_rot_array(q, v)
%QUAT_ROT Summary of this function goes here
%   Detailed explanation goes here

N = size(q, 2);

v_rot = zeros(3, N);

for i = 1 : N
    v_rot(:, i) = quat_rot(q(:, i), v(:, i));
end

end

