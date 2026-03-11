function [qq_star] = qq_conj(qq)
%QQ_CONJ Summary of this function goes here
%   Detailed explanation goes here

qq_star = [q_conj(qq(1:4)); q_conj(qq(5:8))];

end

