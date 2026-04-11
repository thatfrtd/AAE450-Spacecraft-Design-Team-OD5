function results = run_transient(nodes, C_vec, G_matrix, ...
                                  geom, surf, orbit, env, ...
                                  scenario, components, power)
% RUN_TRANSIENT  Integrate thermal ODE system over multiple orbits
%
%  DESCRIPTION:
%    Solves the coupled nonlinear ODE system:
%
%      C_i * dT_i/dt = Q_int_i(t) + Q_ext_i(t)             [internal + external loads]
%                    + sum_j G_ij*(T_j - T_i)               [conduction]
%                    + sum_j sigma*eps_ij*A_ij*(T_j^4-T_i^4) [internal radiation]
%                    - sigma*eps_i*A_rad_i*(T_i^4 - T_space^4) [radiation to space]
%
%    Uses MATLAB ode15s (stiff solver) for robustness.
%    Integrates until quasi-steady-state (temperature swing < tol between orbits).
%
%  INPUTS:
%    nodes     — struct array of thermal nodes
%    C_vec     — thermal capacitance vector [J/K]
%    G_matrix  — conductive conductance matrix [W/K]
%    geom      — geometry struct
%    surf      — surface optical properties
%    orbit     — orbit parameters
%    env       — environment (solar flux, Earth IR, albedo)
%    scenario  — hot or cold case parameters
%    components — component definitions
%    power     — power struct (indices, totals)
%
%  OUTPUTS:
%    results   — struct with:
%      .t_s        — time vector [s]
%      .T_K        — temperature matrix [K], N_nodes x N_timesteps
%      .T_C        — temperature matrix [°C]
%      .T_peak_C   — peak temperature per node [°C]
%      .T_min_C    — min temperature per node [°C]
%      .T_mean_C   — mean temperature per node [°C]
%      .hottest_node  — name of hottest node
%      .coldest_node  — name of coldest node
%      .margin_hot    — margin to max op temp [°C]
%      .margin_cold   — margin to min op temp [°C]
%      .Q_ext_hist — external heat flux history
%
%  TODO: Validate solver convergence — check energy balance error < 1%
%  TODO: If stiffness causes issues, try ode23s or reduce max step size
% =========================================================================

sigma = 5.6704e-8;   % W/m^2/K^4
T_space_K = 2.7;     % K — deep space background

N_nodes    = length(nodes);
N_comp     = length(components);
i_shell    = N_nodes;   % shell is last node

%% ---- Radiation-to-space parameters per node ----------------------------
%  Each node radiates to space through the exterior shell or directly
%  For exterior components: use their own area and emissivity
%  For interior components: radiation to space is via shell (no direct path)
%
%  A_rad_i = effective area radiating to space [m^2]
%  eps_rad_i = effective emissivity for that area

A_rad   = zeros(N_nodes, 1);
eps_rad = zeros(N_nodes, 1);

% Shell radiates from its available side area
A_rad(i_shell)   = geom.A_radiator_available;
eps_rad(i_shell) = surf.eps_body;

% Solar panels — back side radiates to space
idx_panels = 2;
A_rad(idx_panels)   = geom.A_panel_one_face * 2;   % both panel backs
eps_rad(idx_panels) = 0.85;   % TODO: confirm panel back emissivity

% Exterior components radiate small amounts directly
% TODO: Add exterior component areas once mechanical layout is finalized
% For now, assume interior components radiate only through shell

%% ---- Integration parameters --------------------------------------------
N_orbits   = 5;            % integrate for N orbits (check SS convergence)
T_end      = N_orbits * orbit.period;
t_span     = [0, T_end];

% Initial conditions: all nodes at 20°C
T0 = ones(N_nodes, 1) * 293.15;   % K

% ODE options
opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-3, ...
              'MaxStep', orbit.period / 200, ...   % ~200 pts per orbit minimum
              'Stats', 'off');

%% ---- ODE function (nested) --------------------------------------------
% This captures all inputs via closure
    function dTdt = thermal_ode(t, T)
        dTdt = zeros(N_nodes, 1);

        %% External heat fluxes at current time
        Q_ext = external_heat_fluxes(geom, surf, env, scenario, t, orbit);

        %% Distribute external fluxes to nodes
        %  Shell node absorbs the bulk of external radiation
        %  TODO: Apportion to exterior components based on their exposed area
        Q_ext_node = zeros(N_nodes, 1);
        Q_ext_node(i_shell) = Q_ext.total_sum_W;   % all to shell for now

        %% Internal heat dissipation at current time
        Q_int = zeros(N_nodes, 1);
        orbit_frac = mod(t, orbit.period) / orbit.period;
        in_eclipse = Q_ext.in_eclipse;

        for ii = 1:N_comp
            % Power dissipation: use nominal during sunlit, allow reduction in eclipse
            % TODO: Add operating mode logic — currently all components run at full power
            %       unless scenario.power_fracs turns them off (cold/safe mode)
            P = components(ii).P_nom_W * scenario.power_fracs(ii);
            Q_int(ii) = P;
        end

        %% Conduction between nodes
        for ii = 1:N_nodes
            for jj = 1:N_nodes
                if G_matrix(ii,jj) > 0
                    dTdt(ii) = dTdt(ii) + G_matrix(ii,jj) * (T(jj) - T(ii));
                end
            end
        end

        %% Radiation to space (from radiating nodes)
        for ii = 1:N_nodes
            if A_rad(ii) > 0
                Q_rad_out = sigma * eps_rad(ii) * A_rad(ii) * ...
                            (T(ii)^4 - T_space_K^4);
                dTdt(ii) = dTdt(ii) - Q_rad_out;
            end
        end

        %% Add internal generation and external loads
        for ii = 1:N_nodes
            dTdt(ii) = dTdt(ii) + Q_int(ii) + Q_ext_node(ii);
        end

        %% Divide by capacitance
        for ii = 1:N_nodes
            if C_vec(ii) > 0
                dTdt(ii) = dTdt(ii) / C_vec(ii);
            else
                dTdt(ii) = 0;
            end
        end
    end

%% ---- Run ODE solver ----------------------------------------------------
fprintf('  Integrating over %.1f orbits (%.1f min each)...\n', ...
    N_orbits, orbit.period/60);

[t_sol, T_sol] = ode15s(@thermal_ode, t_span, T0, opts);

fprintf('  Integration complete: %d time steps\n', length(t_sol));

%% ---- Check quasi-steady-state convergence ------------------------------
%  Compare last two orbits: if max ΔT < 0.5°C, call it converged
T_last_orbit   = T_sol(t_sol > (N_orbits-1)*orbit.period, :);
T_second_last  = T_sol(t_sol > (N_orbits-2)*orbit.period & ...
                        t_sol < (N_orbits-1)*orbit.period, :);

delta_T_max = max(abs(max(T_last_orbit) - max(T_second_last)));
if delta_T_max > 0.5
    warning('THERMAL: Quasi-steady state not reached. ΔT = %.2f°C between last orbits.', ...
        delta_T_max);
    fprintf('  WARNING: Consider increasing N_orbits.\n');
else
    fprintf('  Quasi-steady state reached (ΔT = %.3f°C)\n', delta_T_max);
end

%% ---- Extract last orbit for analysis -----------------------------------
idx_last = t_sol > (N_orbits-1)*orbit.period;
t_last   = t_sol(idx_last) - (N_orbits-1)*orbit.period;
T_last   = T_sol(idx_last, :);   % K

%% ---- Compile results ---------------------------------------------------
results.t_s     = t_last;
results.T_K     = T_last';       % N_nodes x N_timesteps
results.T_C     = T_last' - 273.15;

results.T_peak_C = max(results.T_C, [], 2);
results.T_min_C  = min(results.T_C, [], 2);
results.T_mean_C = mean(results.T_C, 2);

[~, idx_hot]  = max(results.T_peak_C);
[~, idx_cold] = min(results.T_min_C);
results.hottest_node  = nodes(idx_hot).name;
results.coldest_node  = nodes(idx_cold).name;

% Temperature margins
results.margin_hot  = zeros(N_nodes, 1);   % positive = margin to max
results.margin_cold = zeros(N_nodes, 1);   % positive = margin to min

for i = 1:N_nodes
    results.margin_hot(i)  = nodes(i).T_max_op_K - 273.15 - results.T_peak_C(i);
    results.margin_cold(i) = results.T_min_C(i) - (nodes(i).T_min_op_K - 273.15);
end

results.nodes  = nodes;
results.N_nodes = N_nodes;
results.full_t_s = t_sol;
results.full_T_C = T_sol' - 273.15;

end
