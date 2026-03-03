function [x_hill] = ECI_to_Hill(x_c, x_d)
    % Convert from chief and deputy cartesian orbits to relative orbit of
    % deputy w.r.t. chief in Hill frame (chief orbital frame)
    RTN_to_ECI_DCM = RTN_to_ECI_array(x_c(1:3, :), x_c(4:6, :));
    ECI_to_RTN_DCM = pagetranspose(RTN_to_ECI_DCM);

    h_mag = vecnorm(cross(x_c(1:3, :), x_c(4:6, :)));
    r_mag = vecnorm(x_c(1:3, :));
    omega_c = h_mag ./ r_mag .^ 2; % Orbital angular velocity

    x_relative = reshape(x_d(1:6, :) - x_c(1:6, :), 6, 1, []);

    r_hill = pagemtimes(ECI_to_RTN_DCM, x_relative(1:3, :, :));
    v_hill = pagemtimes(ECI_to_RTN_DCM, (x_relative(4:6, :, :)) - cross(reshape([zeros([2, numel(omega_c)]); omega_c], 3, 1, []), r_hill(1:3, :, :)));
    x_hill = [reshape(r_hill, 3, []); 
              reshape(v_hill, 3, []); 
              x_d(7, :)];
end