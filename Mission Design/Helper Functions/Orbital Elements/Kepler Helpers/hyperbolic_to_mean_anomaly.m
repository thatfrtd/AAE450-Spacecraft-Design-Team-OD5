function M = hyperbolic_to_mean_anomaly(H,e)
    M = e .* sinh(H) - H;
end