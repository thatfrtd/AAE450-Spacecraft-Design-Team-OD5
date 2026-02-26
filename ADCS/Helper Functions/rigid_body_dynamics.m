function dwdt = rigid_body_dynamics(omega, J, u)
% Rigid Body Dynamics
% J: 3x3 Inertia Tensor
% u: 3x1 Net Torque (Control + Disturbances)

    % dw/dt = J^-1 * (u - omega x (J * omega))
    dwdt = J \ (u - cross(omega, J * omega));
end