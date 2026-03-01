function [q_star] = q_conj(q)
%Q_CONJ Quaternion conjugate
%   Detailed explanation goes here

q_star = [-q(1:3, :); q(4, :)];

end

