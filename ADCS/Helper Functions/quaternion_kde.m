function deps = quaternion_kde(t, epsilon, omega)
% Quaternion KDE

    % Extract angular velocity components for readability
    w1 = omega(1);
    w2 = omega(2);
    w3 = omega(3);

    Omega = [  0,   w3, -w2,  w1;
             -w3,   0,   w1,  w2;
              w2, -w1,   0,   w3;
             -w1, -w2, -w3,   0];

    % Compute the derivative
    deps = 0.5 * Omega * epsilon;
end