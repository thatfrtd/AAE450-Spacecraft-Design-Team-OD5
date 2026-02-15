function [u_star, alpha_star, beta_star] = QLaw_thrust_mapping(D, F)
%QLAW_THRUST_MAPPING Calculate optimal thrusting angles to minimize Q-function
%   Detailed explanation goes here
arguments
    D % Q derivative w.r.t. RTN directions (D1, D2, D3)
    F % Force applied
end

alpha_star = atan2(-D(1), -D(2));
beta_star = atan(-D(3) / sqrt(D(1) ^ 2 + D(2) ^ 2));

u_star = F * [cos(beta_star) * sin(alpha_star); % radial
              cos(beta_star) * cos(alpha_star); % theta
              sin(beta_star)]; % normal
% 
% Qdot_ab = D * [cos(beta_star) * sin(alpha_star); ...
%                cos(beta_star) * cos(alpha_star);
%                sin(beta_star)];
end