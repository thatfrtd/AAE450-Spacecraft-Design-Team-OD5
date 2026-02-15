function E = mean_to_eccentric_anomaly(M,e)
    tol = 1e-10;

    % Solve Kepler's equation M = E - e sin(E)
    E(1:numel(M)) = M;
    for index = 1:numel(M)
        while abs(M(index) - (E(index) - e(index) * sin(E(index)))) > tol
            E(index) = E(index) - (E(index) - e(index) * sin(E(index)) - M(index)) / (1 - e(index) * cos(E(index)));
        end
    end
end