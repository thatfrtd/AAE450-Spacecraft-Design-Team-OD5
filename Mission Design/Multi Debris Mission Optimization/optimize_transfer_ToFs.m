function [ToF_best, dVs_best] = optimize_transfer_ToFs(paretos_ToF, paretos_dV, ToF_bounds, max_t)
%OPTIMIZE_TRANSFER_TOFS Summary of this function goes here
%   Detailed explanation goes here
arguments
    paretos_ToF
    paretos_dV
    ToF_bounds
    max_t
end

opts = optimoptions("fmincon","Display","none");
ToF_0 = ToF_bounds(2, :)';
ToF_best = fmincon(@(x) transfer_ToF_objective(x, paretos_ToF, paretos_dV), ToF_0, ...
            ones(size(ToF_0))', max_t, [], [], ToF_bounds(1, :), ToF_bounds(2, :), [], opts);
dVs_best = interp_paretos(ToF_best, paretos_ToF, paretos_dV);
end

function [J] = transfer_ToF_objective(ToFs, pareto_ToF, pareto_dV)
    % Min dV objective
    J = sum(interp_paretos(ToFs, pareto_ToF, pareto_dV));
end

function [dVs] = interp_paretos(ToFs, paretos_ToF, paretos_dV)
    dVs = zeros(size(ToFs));
    for t = 1 : numel(ToFs)
        dVs(t) = interp1(paretos_ToF(:, t), paretos_dV(:, t), ToFs(t), "linear");
    end
end