function [q_cross_mat] = q_cross(q)
%Q_CROSS Summary of this function goes here
%   Detailed explanation goes here

q_cross_mat = [q(4) * eye(3) + skew(q(1:3)), q(1:3); zeros([1, 3]), 0];

end

