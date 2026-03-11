function [R_C_S] = calculate_RCS_allocation_matrix(thrust_vector, thrust_position)
%CALCULATE_RCS_ALLOCATION_MATRIX Calculates RCS thruster force and torque allocation matrix
%   [F_net; tau_net] = R_C_S * T_RCS where R_C_S is the allocation matrix
%   and T_RCS is a vector of directional thrusts.
arguments
    thrust_vector
    thrust_position
end

R_C_S = [thrust_vector; % Thrust
         cross(thrust_position, thrust_vector)]; % Torque

end