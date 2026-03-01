function [qq_cross_mat] = qq_cross(qq)
%QQ_CROSS Summary of this function goes here
%   Detailed explanation goes here

q_A = qq(1:4);
q_B = qq(5:8);

q_A_cross = q_cross(q_A);
q_B_cross = q_cross(q_B);

qq_cross_mat = [q_A_cross, zeros(4, 4); q_B_cross, q_A_cross];

end

