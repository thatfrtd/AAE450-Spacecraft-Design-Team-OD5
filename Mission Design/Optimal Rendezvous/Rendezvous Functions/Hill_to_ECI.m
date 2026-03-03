function [x_d] = Hill_to_ECI(x_c, x_hill)
% Convert to chief and deputy cartesian orbits from relative orbit of
% deputy w.r.t. chief in Hill frame (chief orbital frame)
RTN_to_ECI_DCM = RTN_to_ECI_array(x_c(1:3, :), x_c(4:6, :));

h_mag = vecnorm(cross(x_c(1:3, :), x_c(4:6, :)));
r_mag = vecnorm(x_c(1:3, :));
omega_c = h_mag ./ r_mag .^ 2; % Orbital angular velocity

r_relative = pagemtimes(RTN_to_ECI_DCM, x_hill(1:3, :, :));
v_relative = pagemtimes(RTN_to_ECI_DCM, (x_hill(4:6, :, :) + cross([0; 0; omega_c], x_hill(1:3, :, :))));
x_relative = [reshape(r_relative, 3, []); 
              reshape(v_relative, 3, [])];

x_d = [reshape(x_relative + x_c(1:6, :), 6, 1, []);
       x_hill(7, :)]; % Mass
end