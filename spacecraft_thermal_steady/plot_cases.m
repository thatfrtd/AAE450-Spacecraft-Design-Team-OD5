function plot_cases(comps, T_hot, T_cold, cold_limit, ...
                    hot_viol, cold_viol, ...
                    A_rad, A_rad_avail, ...
                    Q_heat, A_heat_cm2, T_setpt, ...
                    Q_heat_total, margin_C, ...
                    T_skin_hot_C, T_skin_cold_C, Q_leak_in, Q_leak_out)
% PLOT_CASES  Four-figure thermal analysis summary
%
%  Fig 1 — Hot case component temperatures vs max operating limits
%  Fig 2 — Cold case component temperatures vs applicable cold limits
%  Fig 3 — Required radiator area per violating component
%  Fig 4 — Required heater power per violating component
%
%  Color convention:
%    Green  = within thermal limits (pass)
%    Red    = violation (fail)
%    Orange = radiator requirements
%    Blue   = heater requirements
% ==========================================================================

N = length(comps);
x = 1:N;

% Shorten names for cleaner x-axis labels
short = cell(1, N);
for i = 1:N
    s = comps(i).name;
    s = strrep(s, ' [CONFIRM w/ ADCS]', '');
    s = strrep(s, ' [CONFIRM]', '');
    s = strrep(s, '(NEXT-C)', '');
    s = strrep(s, ' (2x)', ' x2');
    s = strrep(s, ' (16x)', ' x16');
    s = strrep(s, ' (4x)', ' x4');
    short{i} = strtrim(s);
end

%% ---- Figure 1: Hot Case ------------------------------------------------
fig1 = figure('Name','Fig 1 — Hot Case Temperatures','NumberTitle','off', ...
              'Position',[50 550 1200 500]);
hold on; grid on;

for i = 1:N
    if hot_viol(i)
        clr = [0.85 0.20 0.15];
    else
        clr = [0.22 0.63 0.30];
    end
    bar(i, T_hot(i), 0.6, 'FaceColor', clr, 'EdgeColor', 'none');
end

T_max_all = arrayfun(@(c) c.T_op_max_C, comps);
plot(x, T_max_all, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);

% MLI skin temperature reference line
yline(T_skin_hot_C, 'b--', 'LineWidth', 1.4, ...
      'Label', sprintf('MLI skin = %.1f C  (Q_{in} = %.1f W)', T_skin_hot_C, Q_leak_in), ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');

xlabel('Component', 'FontSize', 12);
ylabel('Temperature [C]', 'FontSize', 12);
title(sprintf('HOT CASE — Component Temperatures vs Max Operating Limits  (+%d C margin)', ...
    margin_C), 'FontSize', 12, 'FontWeight', 'bold');
xticks(x); xticklabels(short); xtickangle(32);

h1 = bar(NaN, NaN, 'FaceColor', [0.22 0.63 0.30], 'EdgeColor', 'none');
h2 = bar(NaN, NaN, 'FaceColor', [0.85 0.20 0.15], 'EdgeColor', 'none');
h3 = plot(NaN, NaN, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
legend([h1 h2 h3], {'Pass', 'Fail -> radiator needed', 'Max op. limit'}, ...
       'Location', 'northwest', 'FontSize', 9);

%% ---- Figure 2: Cold Case -----------------------------------------------
fig2 = figure('Name','Fig 2 — Cold Case Temperatures','NumberTitle','off', ...
              'Position',[50 50 1200 500]);
hold on; grid on;

for i = 1:N
    if cold_viol(i)
        clr = [0.85 0.20 0.15];
    else
        clr = [0.22 0.63 0.30];
    end
    bar(i, T_cold(i), 0.6, 'FaceColor', clr, 'EdgeColor', 'none');
end

plot(x, cold_limit, 'b^', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);

% MLI skin temperature reference line
yline(T_skin_cold_C, 'b--', 'LineWidth', 1.4, ...
      'Label', sprintf('MLI skin = %.1f C  (Q_{out} = %.1f W)', T_skin_cold_C, Q_leak_out), ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', ...
      'HandleVisibility', 'off');

xlabel('Component', 'FontSize', 12);
ylabel('Temperature [C]', 'FontSize', 12);
title(sprintf('COLD CASE — Component Temperatures vs Min Limits  (-%d C margin)', ...
    margin_C), 'FontSize', 12, 'FontWeight', 'bold');
xticks(x); xticklabels(short); xtickangle(32);

h4 = bar(NaN, NaN, 'FaceColor', [0.22 0.63 0.30], 'EdgeColor', 'none');
h5 = bar(NaN, NaN, 'FaceColor', [0.85 0.20 0.15], 'EdgeColor', 'none');
h6 = plot(NaN, NaN, 'b^', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
legend([h4 h5 h6], {'Pass', 'Fail -> heater needed', 'Min limit (op or survival)'}, ...
       'Location', 'northeast', 'FontSize', 9);

%% ---- Figure 3: Radiator Requirements -----------------------------------
hot_idx = find(A_rad > 0);
n_rad   = length(hot_idx);

fig3 = figure('Name','Fig 3 — Radiator Requirements','NumberTitle','off', ...
              'Position',[1270 550 700 500]);

if n_rad == 0
    text(0.5, 0.5, 'No radiators required', ...
         'Units','normalized','HorizontalAlignment','center','FontSize',14);
    axis off;
else
    hold on; grid on;
    rad_names = short(hot_idx);
    rad_areas = A_rad(hot_idx);

    bar(1:n_rad, rad_areas, 0.6, 'FaceColor', [0.91 0.49 0.19], 'EdgeColor', 'none');

    yline(A_rad_avail, 'g--', 'LineWidth', 2.0, ...
          'Label', sprintf('Available = %.1f m^2', A_rad_avail), ...
          'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');

    for k = 1:n_rad
        text(k, rad_areas(k)*1.04, sprintf('%.3f m^2', rad_areas(k)), ...
             'HorizontalAlignment','center','FontSize',9);
    end

    xticks(1:n_rad); xticklabels(rad_names); xtickangle(28);
    ylabel('Required Radiator Area [m^2]', 'FontSize', 12);
    title(sprintf('HOT CASE — Radiator Area per Component\nTotal: %.3f m^2  |  Available: %.2f m^2', ...
        sum(rad_areas), A_rad_avail), 'FontSize', 11, 'FontWeight', 'bold');

    if sum(rad_areas) <= A_rad_avail
        str = sprintf('PASS: Body radiators sufficient (%.0f%% used)', ...
            100*sum(rad_areas)/A_rad_avail);
        bg = [0.85 1.0 0.85];
    else
        str = sprintf('WARN: Need +%.2f m^2 (deployable?)', sum(rad_areas)-A_rad_avail);
        bg = [1.0 0.85 0.85];
    end
    annotation('textbox', [0.13 0.80 0.40 0.08], 'String', str, ...
        'BackgroundColor', bg, 'FontSize', 10, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'EdgeColor', [0.4 0.4 0.4]);
end

%% ---- Figure 4: Heater Requirements ------------------------------------
cold_idx = find(Q_heat > 0);
n_heat   = length(cold_idx);

fig4 = figure('Name','Fig 4 — Heater Requirements','NumberTitle','off', ...
              'Position',[1270 50 700 500]);

if n_heat == 0
    text(0.5, 0.5, 'No component heaters required', ...
         'Units','normalized','HorizontalAlignment','center','FontSize',14);
    axis off;
else
    hold on; grid on;
    heat_names = short(cold_idx);
    heat_power = Q_heat(cold_idx);

    bar(1:n_heat, heat_power, 0.6, 'FaceColor', [0.20 0.55 0.85], 'EdgeColor', 'none');

    for k = 1:n_heat
        text(k, heat_power(k)*1.04, sprintf('%.1f W', heat_power(k)), ...
             'HorizontalAlignment','center','FontSize',9);
        % Setpoint label inside bar
        text(k, heat_power(k)*0.45, sprintf('setpt\n%.0f C', T_setpt(cold_idx(k))), ...
             'HorizontalAlignment','center','FontSize',8,'Color','white','FontWeight','bold');
    end

    % Total battery load annotation (includes Q_leak_out)
    yline(Q_heat_total, 'r--', 'LineWidth', 1.8, ...
          'Label', sprintf('Battery load = %.1f W (incl. %.1f W MLI loss)', ...
          Q_heat_total, Q_leak_out), ...
          'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');

    xticks(1:n_heat); xticklabels(heat_names); xtickangle(28);
    ylabel('Required Heater Power [W]', 'FontSize', 12);
    title(sprintf('COLD CASE — Heater Power per Component  (20%% margin incl.)\nTotal battery load in eclipse: %.1f W', ...
        Q_heat_total), 'FontSize', 11, 'FontWeight', 'bold');
end

end
