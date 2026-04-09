function [ToF_best, dVs_best, dV_total, t_total] = optimize_transfer_ToFs_time_varying(paretos, IDs, max_t, ToF_0, dV_deorbit, spacecraft_mass, debris_mass)
%OPTIMIZE_TRANSFER_TOFS Summary of this function goes here
%   Detailed explanation goes here
arguments
    paretos
    IDs
    max_t
    ToF_0
    dV_deorbit
    spacecraft_mass
    debris_mass
end

opts = optimoptions("fmincon","Display","iter-detailed");
ToF_best = fmincon(@(x) transfer_ToF_objective(IDs, x, paretos, dV_deorbit, spacecraft_mass, debris_mass), ToF_0, ...
            ones(size(ToF_0))', max_t, [], [], [], [], @(x) spacecraft_routing_nonlconstraints(IDs, x, paretos, max_t, size(IDs, 2)), opts);
dVs_best = interp_paretos_time_varying(IDs, ToF_best', paretos);
[~, ~, t_total] = spacecraft_routing_nonlconstraints(IDs, ToF_best', paretos, max_t, numel(IDs));
[~, dV_total] = transfer_ToF_objective(IDs, ToF_best', paretos, dV_deorbit, spacecraft_mass, debris_mass);
end

function [c, ceq, t_per_sc] = spacecraft_routing_nonlconstraints(IDs, ToFs, paretos, max_t, N_debris_max)
    % Constraints:
    %   * One ship per debris (ineq) - not equality because sum debris
    %     might not be visited so number is 0 or 1
    %   * Maximum total time (ineq)
    %   * Maximum total Delta V (ineq)
    %   * ToF within the Pareto bounds (mostly care about mininimum) (ineq)
    % All inequality constraints are written so satisfied when <= 0
    N_ships = 1;
    N_debris = numel(IDs);

    debris_matrix = zeros([N_ships, N_debris]);
    for s = 1 : N_ships
        for i = 1 : N_debris
            debris_matrix(s, i) = any(IDs(s, :) == i);
        end
    end

    ToF_deorbit = zeros([N_ships, 1]);
    for s = 1 : N_ships
        employed_i = find(IDs(s, :) ~= 0);  % Should always be in order
        for t = 1 : numel(employed_i)
            % First element of .t is the deorbit time which doesn't depend
            % on destination
            if IDs(s, t) ~= 1
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 1}(1);
            else
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 2}(1);
            end
        end
    end
    t_per_sc = sum(ToFs, 2) + ToF_deorbit;
    c_max_t = t_per_sc - max_t;
    
    % ToF within Pareto bounds - 2 * N_ships * (N_debris_max - 1)
    % NEED TO INTERPOLATE
    c_pareto_min_ToF = zeros([N_ships, N_debris_max - 1]);
    c_pareto_max_ToF = zeros([N_ships, N_debris_max - 1]);
    for s = 1 : N_ships
        t0 = 0;
        employed_i = find(IDs(s, :) ~= 0);  % Should always be in order
        for t = 1 : numel(employed_i) - 1
            if IDs(s, t) ~= IDs(s, t + 1)
                ID1 = IDs(s, t); % From 
                ID2 = IDs(s, t + 1); % To
                t0 = t0 + paretos.t{ID1, ID2}(1);

                % Calculate starting time of transfer
                if t ~= 1
                    ToF_total_t = sum(ToFs(s, 1 : (t - 1)));
                else
                    ToF_total_t = 0;
                end
                t0_transfer = min(ToF_total_t + t0, paretos.t{ID1, ID2}(end));

                c_pareto_min_ToF(s, t) = interp1(paretos.t{ID1, IDs(s, t + 1)}, paretos.ToF{ID1, ID2}(1, :), t0_transfer) - ToFs(t);
                c_pareto_max_ToF(s, t) = ToFs(t) - interp1(paretos.t{ID1, ID2}, paretos.ToF{ID1, ID2}(end, :), t0_transfer);
            end
        end
    end
    c_pareto_ToF = [c_pareto_min_ToF(:); c_pareto_max_ToF(:)];
    
    % Package constraints
    c = [c_max_t;
         c_pareto_ToF];
    ceq = []; % No equality constraints
end

function [J, dV_per_sc] = transfer_ToF_objective(IDs, ToFs, paretos, dV_deorbit, spacecraft_mass, debris_mass)
    % Min dV objective
    N_ships = 1;
    N_debris = numel(ToFs) + 1;

    dVs = interp_paretos_time_varying(IDs, numel(ToFs) + 1, paretos);

    debris_matrix = zeros([N_ships, N_debris]);
    for s = 1 : N_ships
        for i = 1 : N_debris
            debris_matrix(s, i) = any(IDs(s, :) == i);
        end
    end

    % Max total Delta V - N_ships (Get approximate adjusted delta V as if
    % the spacecraft's mass never changed from the debris, should be within a few percent)
    dV_per_sc = sum(dVs, 2) + debris_matrix * (dV_deorbit .* (spacecraft_mass + debris_mass) ./ spacecraft_mass);

    J = sum(dV_per_sc);
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

