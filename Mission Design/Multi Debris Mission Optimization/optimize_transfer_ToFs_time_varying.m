function [ToF_best, dVs_best] = optimize_transfer_ToFs_time_varying(paretos_ToF, paretos_dV, ToF_bounds, max_t, options)
%OPTIMIZE_TRANSFER_TOFS Summary of this function goes here
%   Detailed explanation goes here
arguments
    paretos_ToF
    paretos_dV
    ToF_bounds
    max_t
    options.ToF_0 = ToF_bounds(2, :)'
end

opts = optimoptions("fmincon","Display","iter-detailed");
ToF_best = fmincon(@(x) transfer_ToF_objective(x, paretos_ToF, paretos_dV), options.ToF_0, ...
            ones(size(ToF_0))', max_t, [], [], ToF_bounds(1, :), ToF_bounds(2, :), [], opts);
dVs_best = interp_paretos_time_varying(ToF_best, paretos_ToF, paretos_dV);
end

function [J] = transfer_ToF_objective(ToFs, pareto_ToF, pareto_dV)
    % Min dV objective
    J = sum(interp_paretos_time_varying(ToFs, pareto_ToF, pareto_dV));
end

function [dVs] = interp_paretos_time_varying(ToFs, paretos_ToF, paretos_t, paretos_dV)
    dVs = zeros(size(ToFs));
    for t = 1 : numel(ToFs)
        % Interplolate points
        interpolated_paretos_ToF = interp1(paretos_t, paretos_ToF', t0)';
        interpolated_paretos_dV = interp1(paretos_t, paretos_dV', t0)';
    
        % Bound samples
        ToF_bounded = max(min(ToFs, interpolated_paretos_ToF(end)), interpolated_paretos_ToF(1));
    
        % Interpolate
        dVs(t) = interp1(interpolated_paretos_ToF, interpolated_paretos_dV, ToF_bounded, "linear");
    end
end