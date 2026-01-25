function H = true_to_hyperbolic_anomaly(nu,e)
    H = atanh(sin(nu) .* sqrt(e .^ 2 - 1) / (cos(nu) + e));
end