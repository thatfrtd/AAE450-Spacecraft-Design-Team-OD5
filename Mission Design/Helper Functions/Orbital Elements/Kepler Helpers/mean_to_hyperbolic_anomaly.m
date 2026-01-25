function H = mean_to_hyperbolic_anomaly(M,e)
    tol = 1e-12;

    % Solve Kepler's equation M = e sinh(H) - H
    H(1:numel(M)) = M;
    for index = 1:numel(M)
        % Initial guess
        if e < 1.6
            if (-pi < M(index) && M(index) < 0) || M(index) > e
                H(index) = M(index) - e;
            else
                H(index) = M(index) + e;
            end
        else
            if e < 3.6 && abs(M(index)) > pi
                H(index) = M(index) - sign(M(index)) * e;
            else
                H(index) = M(index) / (e - 1);
            end
        end

        % Iterate
        while abs(M(index) - (e * sinh(H(index)) - H(index))) > tol
            H(index) = H(index) + (M(index) - e * sinh(H(index)) + H(index)) / (e * cosh(H(index)) - 1);
        end
    end
end