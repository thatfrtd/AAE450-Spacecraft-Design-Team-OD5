function heater = heater_sizing(results_cold, components)
% HEATER_SIZING  Size survival and operational heaters for cold case
%
%  DESCRIPTION:
%    Determines which components require heaters (active thermal control)
%    in the worst-case cold scenario, and sizes the heater power.
%
%    Heater sizing logic:
%      1. If T_min_simulated < T_min_operating → operational heater needed
%      2. If T_min_simulated < T_min_survival  → CRITICAL — survival heater needed
%
%    Heater power sizing (steady-state):
%      Q_heater = C_thermal * (T_setpoint - T_cold) / t_eclipse
%
%    where:
%      t_eclipse  = eclipse duration [s]
%      T_setpoint = desired temperature [K]
%      T_cold     = predicted minimum temp without heater [K]
%
%    Heaters are typically:
%      - Kapton film heaters (flexible, bonded to component surface)
%      - Setpoint controlled via thermostats or thermistors + relay
%      - Redundant A/B heater circuits
%
%  OUTPUTS:
%    heater.total_power_W          — total heater power budget [W]
%    heater.components_needing_heat — cell array of component names
%    heater.power_per_component_W  — power per component [W]
%    heater.setpoints_C            — heater setpoint temperatures [°C]
%    heater.type                   — 'survival' or 'operational' per component
%
%  NOTE:
%    Heater power is a PARASITIC load — must be included in power budget.
%    In eclipse, heaters are the primary power consumer in safe mode.
%    TODO: Confirm available power in eclipse from batteries with EPS team.
% =========================================================================

N_comp  = length(components);
nodes   = results_cold.nodes;

% Heater setpoint: 5°C above survival min (hysteresis band lower bound)
% Thermostat typically: ON at T_setpoint, OFF at T_setpoint + 5°C  % TODO: confirm hysteresis
setpoint_margin_above_surv = 5;   % °C   % TODO: adjust per component

heater_list   = {};
power_list    = [];
setpoint_list = [];
type_list     = {};

fprintf('  [Heater] Checking cold case temperatures:\n');

for i = 1:N_comp
    c    = components(i);
    T_min_simulated = results_cold.T_min_C(i);
    T_op_min  = c.T_op_min_C;
    T_surv_min = c.T_surv_min_C;

    fprintf('    %-25s  T_min=%.1f°C  (op_min=%.0f°C, surv_min=%.0f°C)\n', ...
        c.name, T_min_simulated, T_op_min, T_surv_min);

    need_heater = false;
    heater_type = '';

    if T_min_simulated < T_surv_min
        fprintf('      *** CRITICAL: Below survival temp! Heater REQUIRED. ***\n');
        need_heater = true;
        heater_type = 'survival';
        T_setpoint  = T_surv_min + setpoint_margin_above_surv;

    elseif T_min_simulated < T_op_min
        fprintf('      NOTE: Below operating temp. Operational heater recommended.\n');
        need_heater = true;
        heater_type = 'operational';
        T_setpoint  = T_op_min + setpoint_margin_above_surv;
    end

    if need_heater
        % Heater power: must supply enough energy during eclipse to prevent
        % component from dropping below setpoint.
        %
        % Simple sizing: P_heater = C_thermal * (T_setpoint - T_min_simulated) / t_eclipse
        %
        % This assumes heater runs continuously in eclipse to make up the shortfall.
        % A more accurate model integrates the ODE with heater ON/OFF logic.
        % TODO: Refine with actual duty cycle once heater control law is defined.

        t_eclipse_s = 0.40 * 5765;   % 40% of 1100 km orbit period (s)  % TODO: use orbit.eclipse_dur
        delta_T     = T_setpoint - T_min_simulated;   % °C
        C_th        = c.C_thermal_J_K;                % J/K

        P_heater = C_th * delta_T / t_eclipse_s;

        % Minimum heater power floor (practical minimum: 2W for Kapton heater)
        P_heater = max(P_heater, 2.0);   % W  % TODO: adjust per component size

        % Add 20% design margin
        P_heater = P_heater * 1.20;

        fprintf('      Heater power: %.1f W (setpoint = %.0f°C, type = %s)\n', ...
            P_heater, T_setpoint, heater_type);

        heater_list{end+1}   = c.name;
        power_list(end+1)    = P_heater;
        setpoint_list(end+1) = T_setpoint;
        type_list{end+1}     = heater_type;
    end
end

%% ---- Battery special case ----------------------------------------------
%  Battery is the most temperature-sensitive component (10–25°C narrow range).
%  Even if simulation doesn't show it going below 10°C, recommend a heater
%  because: (1) uncertainty in model, (2) battery degradation starts at ~5°C.
%
%  TODO: Verify battery thermal isolation strategy with EPS team.
%  Strategy options: (a) insulate battery + small heater, (b) mount near other
%  heat-generating components for passive warming.

idx_battery = 3;   % index in components struct
T_bat_min = results_cold.T_min_C(idx_battery);
if T_bat_min > components(idx_battery).T_op_min_C
    fprintf('    NOTE: Battery is OK in simulation, but a small heater is\n');
    fprintf('          STRONGLY RECOMMENDED as a precaution (Li-ion sensitivity).\n');
    fprintf('          Suggest 5W battery heater as insurance.\n');
    % Add precautionary heater
    heater_list{end+1}   = 'Battery_precautionary';
    power_list(end+1)    = 5.0;   % W — conservative precaution
    setpoint_list(end+1) = 15.0;  % °C
    type_list{end+1}     = 'precautionary';
end

%% ---- Assemble output ---------------------------------------------------
heater.components_needing_heat  = heater_list;
heater.power_per_component_W    = power_list;
heater.setpoints_C              = setpoint_list;
heater.type                     = type_list;
heater.total_power_W            = sum(power_list);
heater.n_heaters                = length(heater_list);

fprintf('\n  [Heater] Total heater power budget: %.1f W\n', heater.total_power_W);
fprintf('  [Heater] NOTE: This power must come from batteries during eclipse.\n');
fprintf('  [Heater] TODO: Confirm battery capacity can support %.1f W for %.0f min.\n', ...
    heater.total_power_W, 0.40*5765/60);

end
