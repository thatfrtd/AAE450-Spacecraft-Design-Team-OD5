function plot_results(results_hot, results_cold, radiator, heater, components, orbit, ...
                      results_hot_ctrl, results_cold_ctrl)
% PLOT_RESULTS  Generate all figures for thermal analysis deliverable
%
%  DESCRIPTION:
%    Produces publication-quality figures suitable for a presentation:
%      Fig 1 — Temperature vs time (hot case, all nodes, 1 orbit)
%      Fig 2 — Temperature vs time (cold case, all nodes, 1 orbit)
%      Fig 3 — Hot case peak temperatures with operating limits
%      Fig 4 — Cold case min temperatures with operating limits
%      Fig 5 — Thermal margins bar chart (hot and cold)
%      Fig 6 — Radiator sizing summary
%      Fig 7 — Heater budget summary
%      Fig 8 — Per-node before/after: hot case with radiators applied
%      Fig 9 — Per-node before/after: cold case with heaters applied
%      Fig 10 — Controlled margin comparison (before vs after, both cases)
%
%  INPUTS:
%    results_hot, results_cold        — uncontrolled transient results
%    radiator, heater                 — sizing structs
%    components, orbit                — spacecraft definitions
%    results_hot_ctrl, results_cold_ctrl — controlled transient results
%                                         (from run_transient_controlled)
%
%  TODO: Adjust color scheme to match your team's presentation template.
% =========================================================================

N_nodes = results_hot.N_nodes;
t_min   = results_hot.t_s / 60;   % convert to minutes
t_min_c = results_cold.t_s / 60;

% Color scheme
colors = lines(N_nodes);
node_names = {results_hot.nodes.name};
% Clean up names for labels
node_labels = strrep(node_names, '_', ' ');

%% ---- Figure 1: Hot Case Temperatures -----------------------------------
fig1 = figure('Name','Hot Case Temperature History','Position',[100 100 1200 500]);
hold on; grid on;
for i = 1:N_nodes
    plot(t_min, results_hot.T_C(i,:), 'Color', colors(i,:), ...
         'LineWidth', 1.5, 'DisplayName', node_labels{i});
end
xlabel('Time in Orbit [min]', 'FontSize', 12);
ylabel('Temperature [°C]', 'FontSize', 12);
title('HOT CASE — Node Temperatures (1 Orbit)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 9);
% Add eclipse shading
eclipse_start = orbit.sunlit_dur / 60;
eclipse_end   = orbit.period / 60;
patch([eclipse_start eclipse_end eclipse_end eclipse_start], ...
      [min(ylim) min(ylim) max(ylim) max(ylim)], ...
      [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
      'DisplayName', 'Eclipse');


%% ---- Figure 2: Cold Case Temperatures ----------------------------------
fig2 = figure('Name','Cold Case Temperature History','Position',[100 100 1200 500]);
hold on; grid on;
for i = 1:N_nodes
    plot(t_min_c, results_cold.T_C(i,:), 'Color', colors(i,:), ...
         'LineWidth', 1.5, 'DisplayName', node_labels{i});
end
xlabel('Time in Orbit [min]', 'FontSize', 12);
ylabel('Temperature [°C]', 'FontSize', 12);
title('COLD CASE — Node Temperatures (1 Orbit)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'eastoutside', 'FontSize', 9);
eclipse_start = (1 - 0.40) * orbit.period / 60;
eclipse_end   = orbit.period / 60;
patch([eclipse_start eclipse_end eclipse_end eclipse_start], ...
      [min(ylim) min(ylim) max(ylim) max(ylim)], ...
      [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
      'DisplayName', 'Eclipse');


%% ---- Figure 3: Temperature vs Limits (Hot Case) -------------------------
fig3 = figure('Name','Hot Case — Temps vs Limits','Position',[100 100 1000 600]);
x = 1:N_nodes;
bar(x, results_hot.T_peak_C, 0.5, 'FaceColor', [0.85 0.33 0.10]);
hold on;
T_max_op = arrayfun(@(n) n.T_max_op_K - 273.15, results_hot.nodes);
plot(x, T_max_op, 'r^', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Max Op Limit');
grid on;
xlabel('Node', 'FontSize', 12);
ylabel('Temperature [°C]', 'FontSize', 12);
title('HOT CASE — Peak Temperatures vs Maximum Operating Limits', ...
      'FontSize', 13, 'FontWeight', 'bold');
xticks(x);
xticklabels(node_labels);
xtickangle(30);
legend({'Peak Temp', 'Max Op Limit'}, 'Location', 'northeast');


%% ---- Figure 4: Temperature vs Limits (Cold Case) ------------------------
fig4 = figure('Name','Cold Case — Temps vs Limits','Position',[100 100 1000 600]);
bar(x, results_cold.T_min_C, 0.5, 'FaceColor', [0.13 0.47 0.71]);
hold on;
T_min_op = arrayfun(@(n) n.T_min_op_K - 273.15, results_cold.nodes);
plot(x, T_min_op, 'bv', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Min Op Limit');
T_min_surv = arrayfun(@(n) n.T_min_surv_K - 273.15, results_cold.nodes);
plot(x, T_min_surv, 'k*', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Survival Limit');
grid on;
xlabel('Node', 'FontSize', 12);
ylabel('Temperature [°C]', 'FontSize', 12);
title('COLD CASE — Minimum Temperatures vs Operating & Survival Limits', ...
      'FontSize', 13, 'FontWeight', 'bold');
xticks(x);
xticklabels(node_labels);
xtickangle(30);
legend({'Min Temp', 'Min Op Limit', 'Survival Limit'}, 'Location', 'southeast');


%% ---- Figure 5: Thermal Margin Summary ----------------------------------
fig5 = figure('Name','Thermal Margins','Position',[100 100 1100 500]);

subplot(1,2,1);
margin_hot = results_hot.margin_hot;
bar_hot = bar(x, margin_hot, 0.6);
bar_hot.FaceColor = 'flat';
for ii = 1:length(margin_hot)
    if margin_hot(ii) >= 0
        bar_hot.CData(ii,:) = [0.2 0.7 0.2];
    else
        bar_hot.CData(ii,:) = [0.9 0.1 0.1];
    end
end
yline(0, 'k--', 'LineWidth', 1.5);
grid on; xticks(x); xticklabels(node_labels); xtickangle(35);
ylabel('Margin [°C]'); title('HOT CASE — Margin to Max Op Temp');

subplot(1,2,2);
margin_cold = results_cold.margin_cold;
bar_cold = bar(x, margin_cold, 0.6);
bar_cold.FaceColor = 'flat';
for ii = 1:length(margin_cold)
    if margin_cold(ii) >= 0
        bar_cold.CData(ii,:) = [0.2 0.7 0.2];
    else
        bar_cold.CData(ii,:) = [0.9 0.1 0.1];
    end
end
yline(0, 'k--', 'LineWidth', 1.5);
grid on; xticks(x); xticklabels(node_labels); xtickangle(35);
ylabel('Margin [°C]'); title('COLD CASE — Margin to Min Op Temp');

sgtitle('Thermal Margin Summary (Green = Safe, Red = Violation)', 'FontSize', 13);


%% ---- Figure 6: Radiator Sizing -----------------------------------------
fig6 = figure('Name','Radiator Summary','Position',[100 100 700 500]);

categories  = {'Required', 'Side Body Available', 'Panel Backs Available'};
areas       = [radiator.area_required_m2, radiator.area_available_m2, radiator.A_panel_backs];
bar_colors  = [[0.85 0.33 0.10]; [0.20 0.60 0.20]; [0.12 0.47 0.71]];
b = bar(areas, 0.5);
b.FaceColor = 'flat';
for ii = 1:3; b.CData(ii,:) = bar_colors(ii,:); end
xticks(1:3); xticklabels(categories);
ylabel('Radiator Area [m^2]'); grid on;
title(sprintf('Radiator Area Summary\nRequired: %.2f m^2 | Available (body): %.2f m^2', ...
    radiator.area_required_m2, radiator.area_available_m2), 'FontSize', 12);
if radiator.body_sufficient
    text(1, radiator.area_required_m2*1.05, '\checkmark BODY SUFFICIENT', ...
         'Color', 'green', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
else
    text(1, radiator.area_required_m2*1.05, '✗ NEED MORE AREA', ...
         'Color', 'red', 'FontSize', 13, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
end


%% ---- Figure 7: Heater Budget -------------------------------------------
if heater.n_heaters > 0
    fig7 = figure('Name','Heater Budget','Position',[100 100 800 500]);
    barh(heater.power_per_component_W, 0.5, 'FaceColor', [0.94 0.50 0.18]);
    yticks(1:heater.n_heaters);
    yticklabels(strrep(heater.components_needing_heat, '_', ' '));
    xlabel('Required Heater Power [W]');
    title(sprintf('Heater Budget Summary\nTotal: %.1f W', heater.total_power_W), ...
          'FontSize', 13);
    grid on;
    xline(heater.total_power_W/heater.n_heaters, 'r--', 'Average');
end

%% =========================================================================
%  FIGURES 8–10: CONTROLLED CASES (before vs after thermal controls)
% ==========================================================================

%% ---- Figure 8: Per-node, Hot Case — Uncontrolled vs With Radiators -----
%  One subplot per node. Each subplot shows:
%    - Red dashed line  : uncontrolled temperature (no extra radiator)
%    - Blue solid line  : controlled temperature   (radiator applied)
%    - Horizontal line  : max operating limit
%  Grey shading marks eclipse window.

N_nodes     = results_hot.N_nodes;
node_labels = strrep({results_hot.nodes.name}, '_', ' ');
t_min_h     = results_hot.t_s      / 60;
t_min_hc    = results_hot_ctrl.t_s / 60;
eclipse_start_h = orbit.sunlit_dur / 60;
eclipse_end_h   = orbit.period     / 60;

n_cols = 3;
n_rows = ceil(N_nodes / n_cols);
fig8   = figure('Name','Hot Case — Before vs After Radiators', ...
                'Position', [50 50 1400 250*n_rows]);

for i = 1:N_nodes
    subplot(n_rows, n_cols, i);
    hold on; grid on;

    % Eclipse shading
    yl_lo = min([results_hot.T_C(i,:), results_hot_ctrl.T_C(i,:)]) - 5;
    yl_hi = max([results_hot.T_C(i,:), results_hot_ctrl.T_C(i,:)]) + 5;
    patch([eclipse_start_h eclipse_end_h eclipse_end_h eclipse_start_h], ...
          [yl_lo yl_lo yl_hi yl_hi], [0.85 0.85 0.85], ...
          'FaceAlpha', 0.4, 'EdgeColor', 'none');

    % Temperature curves
    plot(t_min_h,  results_hot.T_C(i,:),      'r--', 'LineWidth', 1.4, ...
         'DisplayName', 'Uncontrolled');
    plot(t_min_hc, results_hot_ctrl.T_C(i,:), 'b-',  'LineWidth', 1.8, ...
         'DisplayName', 'With radiator');

    % Max operating limit
    T_lim = results_hot.nodes(i).T_max_op_K - 273.15;
    yline(T_lim, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Max op limit');

    ylim([yl_lo, max(yl_hi, T_lim + 5)]);
    title(node_labels{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time [min]', 'FontSize', 8);
    ylabel('Temp [°C]',  'FontSize', 8);

    if i == 1
        legend('Location', 'northeast', 'FontSize', 7);
    end
end
sgtitle('HOT CASE — Uncontrolled (red dashed) vs With Radiators (blue solid)', ...
        'FontSize', 13, 'FontWeight', 'bold');

%% ---- Figure 9: Per-node, Cold Case — Uncontrolled vs With Heaters ------
%  Same layout as Fig 8 but for cold case.
%  Shows heater power duty on a secondary y-axis where applicable.

t_min_c    = results_cold.t_s      / 60;
t_min_cc   = results_cold_ctrl.t_s / 60;
eclipse_start_c = (1 - 0.40) * orbit.period / 60;
eclipse_end_c   = orbit.period / 60;

fig9 = figure('Name','Cold Case — Before vs After Heaters', ...
              'Position', [50 50 1400 250*n_rows]);

for i = 1:N_nodes
    subplot(n_rows, n_cols, i);
    hold on; grid on;

    yl_lo = min([results_cold.T_C(i,:), results_cold_ctrl.T_C(i,:)]) - 5;
    yl_hi = max([results_cold.T_C(i,:), results_cold_ctrl.T_C(i,:)]) + 5;

    % Eclipse shading
    patch([eclipse_start_c eclipse_end_c eclipse_end_c eclipse_start_c], ...
          [yl_lo yl_lo yl_hi yl_hi], [0.75 0.85 1.0], ...
          'FaceAlpha', 0.4, 'EdgeColor', 'none');

    % Temperature curves
    plot(t_min_c,  results_cold.T_C(i,:),      'r--', 'LineWidth', 1.4, ...
         'DisplayName', 'Uncontrolled');
    plot(t_min_cc, results_cold_ctrl.T_C(i,:), 'b-',  'LineWidth', 1.8, ...
         'DisplayName', 'With heater');

    % Min operating limit
    T_op_lim   = results_cold.nodes(i).T_min_op_K   - 273.15;
    T_surv_lim = results_cold.nodes(i).T_min_surv_K - 273.15;
    yline(T_op_lim,   'b:', 'LineWidth', 1.2, 'DisplayName', 'Min op limit');
    yline(T_surv_lim, 'k-', 'LineWidth', 1.0, 'DisplayName', 'Survival limit');

    % Heater power on secondary axis (if this node has a heater)
    if results_cold_ctrl.ctrl.has_heater(i)
        yyaxis right
        area(t_min_cc, results_cold_ctrl.heater_Q_hist(i,:)', ...
             'FaceColor', [1.0 0.7 0.2], 'FaceAlpha', 0.35, ...
             'EdgeColor', 'none', 'DisplayName', 'Heater power [W]');
        ylabel('Heater [W]', 'FontSize', 7, 'Color', [0.8 0.5 0.0]);
        yyaxis left
    end

    ylim([min(yl_lo, T_surv_lim - 5), yl_hi]);
    title(node_labels{i}, 'FontSize', 10, 'FontWeight', 'bold');
    xlabel('Time [min]', 'FontSize', 8);
    ylabel('Temp [°C]',  'FontSize', 8);

    if i == 1
        legend('Location', 'southeast', 'FontSize', 7);
    end
end
sgtitle('COLD CASE — Uncontrolled (red dashed) vs With Heaters (blue solid)', ...
        'FontSize', 13, 'FontWeight', 'bold');

%% ---- Figure 10: Controlled Margin Comparison (before vs after) ---------
%  Four grouped bar charts side by side:
%    Left panel  : hot case margins before and after radiators
%    Right panel : cold case margins before and after heaters
%  Dashed line at zero = pass/fail boundary.
%  A bar above zero = compliant; below = violation.

fig10 = figure('Name','Thermal Margins — Before vs After Controls', ...
               'Position', [50 50 1300 550]);

x = 1:N_nodes;

subplot(1, 2, 1);
margins_before_hot = results_hot.margin_hot;
margins_after_hot  = results_hot_ctrl.margin_hot;
b = bar(x, [margins_before_hot, margins_after_hot], 0.75);
b(1).FaceColor = [0.85 0.33 0.10];   % burnt orange — before
b(2).FaceColor = [0.18 0.60 0.18];   % green        — after
yline(0, 'k--', 'LineWidth', 1.5);
grid on;
xticks(x); xticklabels(node_labels); xtickangle(35);
ylabel('Margin to T_{max,op} [°C]');
title('HOT CASE — Margin Before/After Radiators', 'FontWeight', 'bold');
legend({'Before (no radiator)', 'After (with radiator)'}, 'Location', 'northeast');

subplot(1, 2, 2);
margins_before_cold = results_cold.margin_cold;
margins_after_cold  = results_cold_ctrl.margin_cold;
b2 = bar(x, [margins_before_cold, margins_after_cold], 0.75);
b2(1).FaceColor = [0.13 0.47 0.71];   % blue   — before
b2(2).FaceColor = [0.18 0.60 0.18];   % green  — after
yline(0, 'k--', 'LineWidth', 1.5);
grid on;
xticks(x); xticklabels(node_labels); xtickangle(35);
ylabel('Margin to T_{min,op} [°C]');
title('COLD CASE — Margin Before/After Heaters', 'FontWeight', 'bold');
legend({'Before (no heater)', 'After (with heater)'}, 'Location', 'northeast');

sgtitle('Thermal Margin Summary — Effect of Thermal Controls (green = compliant)', ...
        'FontSize', 13, 'FontWeight', 'bold');

%% ---- Print control specs to console ------------------------------------
print_control_specs(results_hot_ctrl, results_cold_ctrl, heater, radiator, components);

end

%% =========================================================================
%  HELPER: Print thermal control specifications
% ==========================================================================
function print_control_specs(results_hot_ctrl, results_cold_ctrl, heater, radiator, components)
% PRINT_CONTROL_SPECS  Print a formatted spec sheet for all thermal controls
%
%  Called automatically by plot_results. Also useful to call standalone.

fprintf('\n');
fprintf('============================================================\n');
fprintf('  THERMAL CONTROL SPECIFICATIONS\n');
fprintf('============================================================\n');

%% ---- Radiator specs ----------------------------------------------------
fprintf('\n--- RADIATORS (Hot Case) ---\n');
fprintf('  Strategy:         Body-mounted side surface radiation\n');
fprintf('  Required area:    %.3f m^2\n',   radiator.area_required_m2);
fprintf('  Available area:   %.3f m^2\n',   radiator.area_available_m2);
fprintf('  Coverage:         %.1f%% of available side area\n', ...
    (radiator.area_required_m2 / radiator.area_available_m2) * 100);
fprintf('  Operating temp:   %.1f degC  (%.1f K)\n', ...
    radiator.T_rad_C, radiator.T_rad_K);
fprintf('  Rejection cap:    %.1f W/m^2\n', radiator.Q_per_m2);
fprintf('  Total rejection:  %.1f W\n',     radiator.Q_reject_W);
fprintf('  Recommended coat: White paint (eps~0.85-0.90, alpha~0.20)\n');
fprintf('                    OR Optical Solar Reflectors (OSR)\n');
fprintf('                    TODO: Confirm with thermal coating vendor\n');
if radiator.body_sufficient
    fprintf('  Decision:         PASSIVE BODY RADIATORS SUFFICIENT\n');
else
    fprintf('  Decision:         WARNING — additional area needed\n');
end

% Per-node radiator area apportionment
fprintf('\n  Radiator area per node (proportional to heat dissipation):\n');
N_nodes   = results_hot_ctrl.N_nodes;
ctrl      = results_hot_ctrl.ctrl;
total_P   = sum([components.P_nom_W]);
for i = 1:length(components)
    frac   = components(i).P_nom_W / max(total_P, 1);
    A_node = radiator.area_required_m2 * frac;
    if A_node > 0.001
        fprintf('    %-25s  %.4f m^2  (%.0f W dissipated)\n', ...
            components(i).name, A_node, components(i).P_nom_W);
    end
end

% Hot case margins after radiators
fprintf('\n  Hot case margins AFTER radiators applied:\n');
fprintf('    %-25s  %8s  %8s  %8s\n', 'Node', 'Peak°C', 'Limit°C', 'Margin°C');
fprintf('    %s\n', repmat('-', 1, 60));
for i = 1:N_nodes
    T_pk  = results_hot_ctrl.T_peak_C(i);
    T_lim = results_hot_ctrl.nodes(i).T_max_op_K - 273.15;
    marg  = T_lim - T_pk;
    flag  = '';
    if marg < 0;  flag = '  *** VIOLATION'; end
    if marg < 10; flag = '  [low margin]';  end
    fprintf('    %-25s  %8.1f  %8.1f  %8.1f%s\n', ...
        results_hot_ctrl.nodes(i).name, T_pk, T_lim, marg, flag);
end

%% ---- Heater specs ------------------------------------------------------
fprintf('\n--- HEATERS (Cold Case) ---\n');
fprintf('  Type:             Kapton film heaters (recommended)\n');
fprintf('  Control:          Thermostat ON/OFF with %.1f degC hysteresis\n', ...
    results_cold_ctrl.ctrl.heater_hyst_K);
fprintf('  Total budget:     %.1f W\n', heater.total_power_W);
fprintf('  Redundancy:       TODO — confirm A/B redundant circuits with EPS\n\n');

fprintf('  %-25s  %8s  %8s  %10s  %12s\n', ...
    'Component', 'Power(W)', 'Setpt(°C)', 'Type', 'Heater area');
fprintf('  %s\n', repmat('-', 1, 78));
for h = 1:heater.n_heaters
    % Estimate heater area: Kapton heaters typically rated at ~0.5–2 W/cm^2
    % Use 1.0 W/cm^2 as nominal surface power density  % TODO: confirm with heater vendor
    W_per_cm2    = 1.0;
    A_heater_cm2 = heater.power_per_component_W(h) / W_per_cm2;
    fprintf('  %-25s  %8.1f  %8.1f  %10s  %8.1f cm^2\n', ...
        heater.components_needing_heat{h}, ...
        heater.power_per_component_W(h), ...
        heater.setpoints_C(h), ...
        heater.type{h}, ...
        A_heater_cm2);
end

% Cold case margins after heaters
fprintf('\n  Cold case margins AFTER heaters applied:\n');
fprintf('    %-25s  %8s  %8s  %8s\n', 'Node', 'Min°C', 'Limit°C', 'Margin°C');
fprintf('    %s\n', repmat('-', 1, 60));
for i = 1:N_nodes
    T_mn  = results_cold_ctrl.T_min_C(i);
    T_lim = results_cold_ctrl.nodes(i).T_min_op_K - 273.15;
    marg  = T_mn - T_lim;
    flag  = '';
    if marg < 0;  flag = '  *** VIOLATION'; end
    if marg < 10; flag = '  [low margin]';  end
    fprintf('    %-25s  %8.1f  %8.1f  %8.1f%s\n', ...
        results_cold_ctrl.nodes(i).name, T_mn, T_lim, marg, flag);
end

fprintf('\n  NOTE: Eclipse duration used for heater sizing: %.0f min\n', ...
    0.40*5765/60);
fprintf('  TODO: Confirm eclipse battery capacity covers %.1f W for %.0f min\n', ...
    heater.total_power_W, 0.40*5765/60);
fprintf('============================================================\n\n');

end
