function [x2_cart_ck, x1_kep, M2_ck] = propagate_conic_R(x1_cart, ToF, mu, options)
    arguments
        x1_cart
        ToF
        mu
        options.n_pts = 200;
    end
    
    % Extract elements
    x1_kep = cartesian_to_keplerian(x1_cart, [0;0;1], [1;0;0], mu);
    a = x1_kep(1); e = x1_kep(2); i = x1_kep(3);
    RAAN = x1_kep(4); omega = x1_kep(5); M1 = x1_kep(6);
    
    % Mean motion and final true anomaly
    n = sqrt(mu / abs(a)^3);
    %M1 = keplerian_to_mean_anomaly(theta1, e);
    M2_ck = M1 + n * ToF;
    E2_ck = mean_to_eccentric_anomaly(M2_ck,e);
    theta2 = eccentric_to_true_anomaly(E2_ck,e);
    
    % Orbital parameter
    p = a * (1 - e^2);
    
    % Build rotation matrix
    theta = theta2 + omega;
    ICR = [cos(RAAN)*cos(theta)-sin(RAAN)*cos(i)*sin(theta), -cos(RAAN)*sin(theta)-sin(RAAN)*cos(i)*cos(theta), sin(RAAN)*sin(i);
           sin(RAAN)*cos(theta)+cos(RAAN)*cos(i)*sin(theta), -sin(RAAN)*sin(theta)+cos(RAAN)*cos(i)*cos(theta), -cos(RAAN)*sin(i);
           sin(i)*sin(theta), sin(i)*cos(theta), cos(i)];
    
    % Radius and position
    r = p / (1 + e*cos(theta2 - omega));
    r_vec = ICR * [r; 0; 0];
    
    % Velocity in perifocal frame
    v_pf = sqrt(mu/p) * [-sin(theta2); e + cos(theta2); 0];
    v_vec = ICR * v_pf;
    
    x2_cart_ck = [r_vec; v_vec];
end
