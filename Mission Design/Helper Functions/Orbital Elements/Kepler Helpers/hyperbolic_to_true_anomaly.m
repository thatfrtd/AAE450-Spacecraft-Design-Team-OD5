function nu = hyperbolic_to_true_anomaly(H,e)
    nu = atanh(-sinh(H) .* sqrt(e .^ 2 - 1)./ (1-e*cosh(H)));
end