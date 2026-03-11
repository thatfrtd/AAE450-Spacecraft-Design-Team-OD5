function [x_d] = Hill_to_ECI(x_c, x_hill)
% Convert to chief and deputy cartesian orbits from relative orbit of
% deputy w.r.t. chief in Hill frame (chief orbital frame)
RTN_to_ECI_DCM = RTN_to_ECI_array(x_c(1:3, :), x_c(4:6, :));

h = cross(x_c(1:3, :), x_c(4:6, :));
r_mag = vecnorm(x_c(1:3, :));
omega_c = h ./ r_mag .^ 2; % Orbital angular velocity

r_relative = pagemtimes(RTN_to_ECI_DCM, x_hill(1:3, :, :));
r_ECI = reshape(r_relative, 3, 1, []) + x_c(1:3, :);
v_relative = reshape(pagemtimes(RTN_to_ECI_DCM, x_hill(4:6, :, :)), 3, 1, []) + cross(omega_c, r_ECI);
v_ECI = v_relative + x_c(4:6, :);

x_d = [r_ECI;
       v_ECI;
       x_hill(7, :)]; % Mass
end