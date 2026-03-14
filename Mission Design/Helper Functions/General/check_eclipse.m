function [eclipse_val, eclipsed, on_close_side] = check_eclipse(t, rvec, R_E, a_E, n_E)
    % Check if spacecraft is eclipsed by the Earth
    % Should get solar/Earth position from ephemeris...
    % How to properly account for Earth's tilt?
    axial_tilt = deg2rad(23.44);

    nu = n_E * t;
    sun_vector = [cos(nu); sin(nu); 0]; 

    eclipse_angle = atan2(R_E, a_E);

    rvec_solarframe = eul2rotm([0, -axial_tilt, 0]) * rvec; % Yaw probably shouldn't be 0 but some true anomaly for equinox or something
    rvec_E = sun_vector * a_E;
    rvec_sc = rvec_E + rvec_solarframe;
    r = norm(rvec_sc);
    rhat = rvec_sc / r;
    on_close_side = r < a_E;

    eclipse_val = (dot(rhat, sun_vector) - cos(eclipse_angle)) * 1e10; % > 0
    if on_close_side
        eclipse_val = min(eclipse_val, -1e-1);
    end

    eclipsed = eclipse_val > 0;
end
