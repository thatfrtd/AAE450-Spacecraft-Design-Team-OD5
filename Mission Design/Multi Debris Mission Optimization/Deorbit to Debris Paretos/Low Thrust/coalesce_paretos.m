%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Transfer from a deorbit orbit (after deorbiting debris) to new debris
% Author: Travis Hastreiter 
% Created On: 13 March, 2026
% Description: Orbit transfer using Q-Law from deorbit orbit (after drop 
% off) to new debris not accounting for rendezvous (assuming not much extra 
% delta V and time).
% Most Recent Change: 15 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transfer_dataset_inputs = load("Multi Debris Mission Optimization\transfer_dataset_inputs_fixedreorbit.mat").transfer_dataset_inputs;

folder = "Mission Design\Multi Debris Mission Optimization\Deorbit to Debris Paretos\Low Thrust\";

indices = 1 : 310;

N_debris = numel(transfer_dataset_inputs.debris_ID);

paretos = struct();
paretos.ToF = cell(N_debris, N_debris);
paretos.dV = cell(N_debris, N_debris);
paretos.t = cell(N_debris, N_debris);
for i = indices
    pareto_i = load(folder + "transfer_" + i).pareto;
    IDs_i = transfer_dataset_inputs.IDs(:, i);
    
    [ToF_sorted, ToF_sorted_i] = sort(pareto_i.fval(:, 2));
    for k = 1 : (numel(ToF_sorted) - 1) % Fix duplicate ToFs...
        if any(ToF_sorted(k) == ToF_sorted(k + 1:end))
            ToF_sorted(k + find(ToF_sorted(k) == ToF_sorted(k + 1:end))) = ToF_sorted(k + find(ToF_sorted(k) == ToF_sorted(k + 1:end))) + 1e-10;
        end
    end
    dV_sorted = pareto_i.fval(ToF_sorted_i, 1);

    if numel(ToF_sorted) ~= 60
        % Resample to be 60 points, this should only happen for very small
        % transfers anyways
        ToF_sorted_new = linspace(min(ToF_sorted), max(ToF_sorted), 60);
        dV_sorted = interp1(ToF_sorted, dV_sorted, ToF_sorted_new);
        ToF_sorted = ToF_sorted_new;
    end
  
    paretos.ToF{IDs_i(1), IDs_i(2)}(:, end + 1) = ToF_sorted / 365.25;
    paretos.dV{IDs_i(1), IDs_i(2)}(:, end + 1) = dV_sorted;
    paretos.t{IDs_i(1), IDs_i(2)}(end + 1) = transfer_dataset_inputs.t0(i);
end
%%
save("Mission Design\Multi Debris Mission Optimization\Deorbit to Debris Paretos\Low Thrust\low_thrust_paretos.mat", "paretos")

%%
ind = [3, 8];

figure
for i = 1 : size(paretos.ToF{ind(1), ind(2)}, 2)
    plot3(paretos.ToF{ind(1), ind(2)}(:, i), paretos.dV{ind(1), ind(2)}(:, i), paretos.t{ind(1), ind(2)}(i) * ones(size(paretos.ToF{ind(1), ind(2)}, 1), 1)); hold on
end

t0_sample = linspace(0.6, 7, 20);
for i = 1 : numel(t0_sample)
    ToF_sample = linspace(0, 2, 100);
    
    dV_sample = interp_paretos_time_varying(paretos.ToF{ind(1), ind(2)}, paretos.t{ind(1), ind(2)}, paretos.dV{ind(1), ind(2)}, ToF_sample, t0_sample(i));
    scatter3(ToF_sample, dV_sample, t0_sample(i) * ones([1, numel(ToF_sample)]))
end

grid on
xlabel("ToF [days]")
ylabel("dV [km / s]")
ylim([0, inf])
zlabel("Initial Time [yr]")
title(sprintf("Debris %g Deorbit to Debris %g Pareto vs Time", transfer_dataset_inputs.debris_ID(ind(1)), transfer_dataset_inputs.debris_ID(ind(2))))

%% Helper Functions
function [dVs] = interp_paretos_time_varying(paretos_ToF, paretos_t, paretos_dV, ToF, t0)
    % Interplolate points
    interpolated_paretos_ToF = interp1(paretos_t, paretos_ToF', t0)';
    interpolated_paretos_dV = interp1(paretos_t, paretos_dV', t0)';

    % Bound samples
    ToF_bounded = max(min(ToF, interpolated_paretos_ToF(end)), interpolated_paretos_ToF(1));

    % Interpolate
    dVs = interp1(interpolated_paretos_ToF, interpolated_paretos_dV, ToF_bounded, "linear");
end