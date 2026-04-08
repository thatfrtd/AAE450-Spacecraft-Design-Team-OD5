function [ToF_best, dVs_best] = optimize_transfer_ToFs_time_varying(paretos, IDs, max_t, ToF_0)
%OPTIMIZE_TRANSFER_TOFS Summary of this function goes here
%   Detailed explanation goes here
arguments
    paretos
    IDs
    max_t
    ToF_0
end

opts = optimoptions("fmincon","Display","iter-detailed");
ToF_best = fmincon(@(x) transfer_ToF_objective(IDs, x, paretos), ToF_0, ...
            ones(size(ToF_0))', max_t, [], [], [], [], [], opts);
dVs_best = interp_paretos_time_varying(IDs, ToF_best, paretos);
end

function [J] = transfer_ToF_objective(IDs, ToFs, paretos)
    % Min dV objective
    J = sum(interp_paretos_time_varying(IDs, ToFs, paretos));
end

function [dVs] = interp_paretos_time_varying(IDs, ToFs, paretos)
    % Need to flatten vars to make interp1 happy....
    % Use IDs into paretos...
    % IDs has shape [N_ships, N_debris_max]
    % ToFs has shape [N_ships, N_debris_max - 1]
    % dVs has shape [N_ships, N_debris_max - 1]
    % Pareto has .ToF and .dV like {N_debris, N_debris}[N_pareto, N_times] 
    % and .t like {N_debris, N_debris}[N_times]

    dVs = zeros(size(ToFs));
    t0 = 0;
    for t = 1 : size(ToFs, 2)
        if IDs(t + 1) ~= 0 && IDs(t) ~= IDs(t + 1)
            ID1 = IDs(t); % From 
            ID2 = IDs(t + 1); % To
            t0 = t0 + paretos.t{ID1, ID2}(1);

            % Calculate starting time of transfer
            if t ~= 1
                ToF_total_t = sum(ToFs(1 : (t - 1)));
            else
                ToF_total_t = 0;
            end
            t0_transfer = min(ToF_total_t + t0, paretos.t{ID1, ID2}(end));

            % Interplolate points
            interpolated_paretos_ToF = interp1(paretos.t{ID1, ID2}, paretos.ToF{ID1, ID2}', t0_transfer)';
            interpolated_paretos_dV = interp1(paretos.t{ID1, ID2}, paretos.dV{ID1, ID2}', t0_transfer)';

            % Bound ToFs
            ToF_bounded = max(min(ToFs(t), interpolated_paretos_ToF(end)), interpolated_paretos_ToF(1));
            
            % Interpolate Pareto
            dVs(t) = interp1(interpolated_paretos_ToF, interpolated_paretos_dV, ToF_bounded, "linear");
        else
            break;
        end
    end
end