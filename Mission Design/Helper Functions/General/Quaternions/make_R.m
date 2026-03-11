function R = make_R(theta, axis)
    R = zeros([3, 3, numel(theta)]);

    for i = 1:numel(theta)
        c_theta = cos(theta(i));
        s_theta = sin(theta(i));
        if axis == 3
            R_i = [c_theta, s_theta, 0; 
                   -s_theta,  c_theta, 0; 
                   0,        0,       1];
        elseif axis == 2
            R_i = [c_theta, 0, -s_theta; 
                   0,       1,       0;
                   s_theta,0, c_theta];
        elseif axis == 1
            R_i = [1,        0,        0;
                   0,  c_theta, s_theta; 
                   0,  -s_theta,  c_theta];
        end

        if numel(theta) > 1
            R(:, :, i) = R_i;
        else
            R = R_i;
        end
    end
end