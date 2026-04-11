function results_ctrl = run_transient_controlled(nodes, C_vec, G_matrix, ...
                                                   geom, surf, orbit, env, ...
                                                   scenario, components, power, ...
                                                   heater, radiator, case_type)
% RUN_TRANSIENT_CONTROLLED  Re-run transient with heaters and/or radiators active
%
%  DESCRIPTION:
%    Identical ODE structure to run_transient, but with two additions:
%
%    FOR COLD CASE ('cold'):
%      Each component that was assigned a heater in heater_sizing.m has a
%      thermostat-controlled heat source added to its node. The heater fires
%      when the node temperature drops below the setpoint and switches off
%      when it recovers to setpoint + hysteresis_band.
%
%    FOR HOT CASE ('hot'):
%      Each component is given additional radiative coupling to space,
%      representing the sized radiator area apportioned per node based on
%      its fraction of total heat dissipation. This models the effect of
%      body-mounted radiators or thermal straps routing heat to the shell.
%
%    The result is compared against the uncontrolled case to show that
%    all margins are satisfied.
%
%  INPUTS:
%    nodes, C_vec, G_matrix, geom, surf, orbit, env, scenario,
%    components, power   — same as run_transient
%    heater              — struct from heater_sizing.m
%    radiator            — struct from radiator_sizing.m
%    case_type           — string: 'hot' or 'cold'
%
%  OUTPUTS:
%    results_ctrl  — same format as run_transient results, plus:
%      .heater_Q_hist   — heater power history per node [W x timesteps]
%      .ctrl_specs      — struct describing what was applied to each node
%
%  NOTE:
%    The thermostat model here is a smooth sigmoid approximation rather
%    than a hard ON/OFF switch. This avoids ODE stiffness issues from
%    discontinuous forcing while still accurately representing average
%    heater power. The sigmoid width is 0.5°C — effectively a bang-bang
%    controller for all practical purposes.
%
%  TODO: Once heater control law is finalized (PWM vs thermostat), update
%        the control logic here to match the flight software.
% =========================================================================

sigma     = 5.6704e-8;
T_space_K = 2.7;

N_nodes = length(nodes);
N_comp  = length(components);
i_shell = N_nodes;

%% ---- Build per-node control specification ------------------------------
ctrl = struct();
ctrl.case_type = case_type;

% Initialize arrays
ctrl.has_heater      = false(N_nodes, 1);
ctrl.heater_power_W  = zeros(N_nodes, 1);
ctrl.heater_setpt_K  = zeros(N_nodes, 1);
ctrl.heater_hyst_K   = 2.0;   % °C hysteresis band  % TODO: match actual thermostat spec
ctrl.extra_A_rad     = zeros(N_nodes, 1);   % additional radiator area per node [m^2]
ctrl.extra_eps_rad   = zeros(N_nodes, 1);

%% ---- COLD CASE: assign heaters to nodes --------------------------------
if strcmp(case_type, 'cold')
    for h = 1:heater.n_heaters
        % Match heater name to component index
        comp_name = heater.components_needing_heat{h};
        node_idx  = find_node_index(nodes, comp_name, N_comp);
        if node_idx > 0
            ctrl.has_heater(node_idx)     = true;
            ctrl.heater_power_W(node_idx) = heater.power_per_component_W(h);
            ctrl.heater_setpt_K(node_idx) = heater.setpoints_C(h) + 273.15;
        end
    end
end

%% ---- HOT CASE: apportion radiator area to nodes ------------------------
if strcmp(case_type, 'hot')
    % Distribute the required radiator area across nodes proportional to
    % each component's heat dissipation fraction.
    % Shell node gets any leftover area.
    %
    % Physical interpretation: this represents either
    %   (a) body-surface radiator panels painted/coated and thermally
    %       coupled via the shell to each component, OR
    %   (b) dedicated radiator panels with thermal straps per component
    %
    % TODO: Once strap routing is decided, replace proportional split with
    %       actual strap conductance + dedicated radiator area per component.

    total_P = sum([components.P_nom_W]);
    if total_P == 0; total_P = 1; end  % avoid divide-by-zero

    A_total = radiator.area_required_m2;

    for i = 1:N_comp
        frac = components(i).P_nom_W / total_P;
        ctrl.extra_A_rad(i)   = A_total * frac;
        ctrl.extra_eps_rad(i) = 0.85;   % white paint / OSR  % TODO: confirm coating
    end
    % Shell gets its baseline area already set in run_transient;
    % here we add the per-component areas to their own nodes so heat
    % is rejected directly rather than via conduction to shell first.
end

%% ---- Radiation-to-space (baseline, same as run_transient) --------------
A_rad   = zeros(N_nodes, 1);
eps_rad = zeros(N_nodes, 1);

A_rad(i_shell)   = geom.A_radiator_available;
eps_rad(i_shell) = surf.eps_body;

idx_panels = 2;
A_rad(idx_panels)   = geom.A_panel_one_face * 2;
eps_rad(idx_panels) = 0.85;

% Add controlled radiator areas
for i = 1:N_nodes
    A_rad(i)   = A_rad(i)   + ctrl.extra_A_rad(i);
    eps_rad(i) = max(eps_rad(i), ctrl.extra_eps_rad(i));
end

%% ---- Integration parameters --------------------------------------------
N_orbits = 5;
T_end    = N_orbits * orbit.period;
t_span   = [0, T_end];
T0       = ones(N_nodes, 1) * 293.15;   % K — 20°C start

opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-3, ...
              'MaxStep', orbit.period / 200, ...
              'Stats', 'off');

% Storage for heater power history (sampled at ODE output points)
heater_Q_log = [];
t_log        = [];

%% ---- ODE function ------------------------------------------------------
    function dTdt = thermal_ode_ctrl(t, T)
        dTdt = zeros(N_nodes, 1);

        %% External heat fluxes
        Q_ext      = external_heat_fluxes(geom, surf, env, scenario, t, orbit);
        Q_ext_node = zeros(N_nodes, 1);
        Q_ext_node(i_shell) = Q_ext.total_sum_W;

        %% Internal dissipation
        Q_int = zeros(N_nodes, 1);
        for ii = 1:N_comp
            P = components(ii).P_nom_W * scenario.power_fracs(ii);
            Q_int(ii) = P;
        end

        %% Heater control (cold case only)
        %  Smooth sigmoid thermostat: P_heater * sigmoid((Tset - T) / width)
        %  At T = Tset: delivers 50% power; at T = Tset-1K: ~88%; at T = Tset+1K: ~12%
        %  This gives effectively ON below setpoint, OFF above, without discontinuity.
        Q_heater = zeros(N_nodes, 1);
        if strcmp(case_type, 'cold')
            sigmoid_width = 0.3;   % K — sharpness of thermostat transition
            for ii = 1:N_nodes
                if ctrl.has_heater(ii)
                    x = (ctrl.heater_setpt_K(ii) - T(ii)) / sigmoid_width;
                    % Clamp to avoid exp overflow
                    x = max(min(x, 50), -50);
                    frac_on = 1 / (1 + exp(-x));
                    Q_heater(ii) = ctrl.heater_power_W(ii) * frac_on;
                end
            end
        end

        %% Conduction
        for ii = 1:N_nodes
            for jj = 1:N_nodes
                if G_matrix(ii,jj) > 0
                    dTdt(ii) = dTdt(ii) + G_matrix(ii,jj) * (T(jj) - T(ii));
                end
            end
        end

        %% Radiation to space
        for ii = 1:N_nodes
            if A_rad(ii) > 0
                Q_rad_out = sigma * eps_rad(ii) * A_rad(ii) * (T(ii)^4 - T_space_K^4);
                dTdt(ii)  = dTdt(ii) - Q_rad_out;
            end
        end

        %% Sum all loads and divide by capacitance
        for ii = 1:N_nodes
            dTdt(ii) = dTdt(ii) + Q_int(ii) + Q_ext_node(ii) + Q_heater(ii);
            if C_vec(ii) > 0
                dTdt(ii) = dTdt(ii) / C_vec(ii);
            else
                dTdt(ii) = 0;
            end
        end
    end

%% ---- Run solver --------------------------------------------------------
fprintf('  Integrating controlled case (%s) over %.1f orbits...\n', ...
    case_type, N_orbits);
[t_sol, T_sol] = ode15s(@thermal_ode_ctrl, t_span, T0, opts);
fprintf('  Done: %d steps\n', length(t_sol));

%% ---- Extract last orbit ------------------------------------------------
idx_last = t_sol > (N_orbits-1)*orbit.period;
t_last   = t_sol(idx_last) - (N_orbits-1)*orbit.period;
T_last   = T_sol(idx_last, :);

%% ---- Recompute heater power at each output point (for plotting) --------
heater_Q_hist = zeros(N_nodes, sum(idx_last));
if strcmp(case_type, 'cold')
    sigmoid_width = 0.3;
    T_last_K = T_last + 273.15;   % already in K from solver
    % Note: T_sol is in K; T_last = T_sol rows, so T_last is already K
    % Actually T0 was in K so T_sol is in K — correct
    for k = 1:size(T_last, 1)
        for ii = 1:N_nodes
            if ctrl.has_heater(ii)
                x = (ctrl.heater_setpt_K(ii) - T_last(k,ii)) / sigmoid_width;
                x = max(min(x, 50), -50);
                frac_on = 1 / (1 + exp(-x));
                heater_Q_hist(ii, k) = ctrl.heater_power_W(ii) * frac_on;
            end
        end
    end
end

%% ---- Compile results ---------------------------------------------------
results_ctrl.t_s      = t_last;
results_ctrl.T_K      = T_last';
results_ctrl.T_C      = T_last' - 273.15;

results_ctrl.T_peak_C = max(results_ctrl.T_C, [], 2);
results_ctrl.T_min_C  = min(results_ctrl.T_C, [], 2);
results_ctrl.T_mean_C = mean(results_ctrl.T_C, 2);

[~, idx_hot]  = max(results_ctrl.T_peak_C);
[~, idx_cold] = min(results_ctrl.T_min_C);
results_ctrl.hottest_node = nodes(idx_hot).name;
results_ctrl.coldest_node = nodes(idx_cold).name;

results_ctrl.margin_hot  = zeros(N_nodes, 1);
results_ctrl.margin_cold = zeros(N_nodes, 1);
for i = 1:N_nodes
    results_ctrl.margin_hot(i)  = nodes(i).T_max_op_K - 273.15 - results_ctrl.T_peak_C(i);
    results_ctrl.margin_cold(i) = results_ctrl.T_min_C(i) - (nodes(i).T_min_op_K - 273.15);
end

results_ctrl.nodes       = nodes;
results_ctrl.N_nodes     = N_nodes;
results_ctrl.heater_Q_hist = heater_Q_hist;
results_ctrl.ctrl        = ctrl;

end

%% ---- Helper: find node index by name -----------------------------------
function idx = find_node_index(nodes, name, N_comp)
    idx = 0;
    % Strip '_precautionary' suffix added by heater_sizing for battery
    search_name = strrep(name, '_precautionary', '');
    for k = 1:N_comp
        if strcmpi(nodes(k).name, search_name)
            idx = k;
            return;
        end
    end
end
