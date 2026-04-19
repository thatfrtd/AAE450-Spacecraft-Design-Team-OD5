%% Change History- 4/18
%  CHANGE SUMMARY

%  CHANGE 1 — Figure 2: MLI energy balance annotation added
%    BEFORE: Figure 2 showed cold case component temperatures only.
%    AFTER:  A coloured annotation box added at the bottom of Figure 2
%            showing whether internal cold-case dissipation covers Q_leak_out.
%            Green = surplus (no bulk heater needed).
%            Red   = deficit (bulk interior heater required, with magnitude).
%    REASON: When no component heaters are required, Q_leak_out — the most
%            important cold-case number — was not visible in any plot.
%            This annotation makes the cold-case energy balance visible
%            on the temperature figure where the reader is already focused.
%
%  CHANGE 2 — Figure 4: Full heater requirements figure overhauled
%    BEFORE: If no component heaters required, figure showed only the text
%            "No component heaters required" with no further information.
%            If component heaters existed, bar chart showed those heaters
%            with no reference to the interior energy balance.
%    AFTER:  Three distinct display cases:
%
%            CASE 1 — No component heaters, no deficit:
%              Confirms internal dissipation covers MLI loss and displays
%              surplus value. No action required.
%
%            CASE 2 — No component heaters, deficit exists:
%              Clearly flags that a bulk interior heater is required,
%              states the deficit magnitude, and explains that without
%              it T_int drops below 20 C, invalidating all component
%              temperature predictions in the model.
%
%            CASE 3 — Component heaters exist:
%              Shows original bar chart plus a green/red annotation box
%              at the bottom indicating whether component heaters plus
%              internal dissipation together cover the MLI loss, or
%              whether an additional bulk interior heater is still needed.
%
%    REASON: The previous figure was silent on the most important cold-case
%            question — whether the interior can sustain 20 C — in the
%            common case where no individual component heaters are triggered.
%            The overhaul ensures Figure 4 always gives a complete and
%            self-contained answer to the cold-case thermal control question.
%
%  CHANGE 3 — Figure 4: Bug fix — Q_heater_budget_total replaced by
%             Q_heater_total in calling argument (fix in main_thermal.m
%             Section 12, reflected here for traceability)
%    BEFORE: Q_heat_total received Q_heater_budget_total from main_thermal.m,
%            which already included Q_leak_out. This caused Q_leak_out to
%            cancel out of the surplus calculation inside this function,
%            making has_deficit always evaluate to false.
%    AFTER:  Q_heat_total now correctly receives Q_heater_total (component
%            heaters only), allowing the surplus check to correctly compare
%            internal dissipation plus component heaters against Q_leak_out.
%    REASON: Deficit flag in Figure 2 was correctly showing a deficit while
%            Figure 4 was simultaneously and incorrectly reporting no bulk
%            heater needed. This fix makes the two figures consistent.

function plot_cases_4_19_edit_most_recent(comps, T_hot, T_cold, cold_limit, ...
    hot_viol, cold_viol, A_rad, A_rad_avail, Q_heat_W, A_heat, ...
    T_setpt, Q_heat_total, margin, T_skin_hot, T_skin_cold, ...
    Q_leak_in, Q_leak_out, Q_int_cold_plot, Q_surplus_plot);
% function plot_cases_4_18_edit(comps, T_hot, T_cold, cold_limit, ...
%                     hot_viol, cold_viol, ...
%                     A_rad, A_rad_avail, ...
%                     Q_heat, A_heat_cm2, T_setpt, ...
%                     Q_heat_total, margin_C, ...
%                     T_skin_hot_C, T_skin_cold_C, Q_leak_in, Q_leak_out)
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
fig1 = figure('Name','Fig 1 — Hot Case Temperatures','NumberTitle','off');
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
yline(T_skin_hot, 'b--', 'LineWidth', 1.4, ...
      'Label', sprintf('MLI skin = %.1f C  (Q_{in} = %.1f W)', T_skin_hot, Q_leak_in), ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', ...
      'HandleVisibility', 'off');

xlabel('Component', 'FontSize', 12);
ylabel('Temperature [C]', 'FontSize', 12);
title(sprintf('HOT CASE — Component Temperatures vs Max Operating Limits  (+%d C margin)', ...
    margin), 'FontSize', 12, 'FontWeight', 'bold');
xticks(x); xticklabels(short); xtickangle(32);

h1 = bar(NaN, NaN, 'FaceColor', [0.22 0.63 0.30], 'EdgeColor', 'none');
h2 = bar(NaN, NaN, 'FaceColor', [0.85 0.20 0.15], 'EdgeColor', 'none');
h3 = plot(NaN, NaN, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
legend([h1 h2 h3], {'Pass', 'Fail -> radiator needed', 'Max op. limit'}, ...
       'Location', 'northwest', 'FontSize', 9);

%% ---- Figure 2: Cold Case -----------------------------------------------
fig2 = figure('Name','Fig 2 — Cold Case Temperatures','NumberTitle','off');
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
yline(T_skin_cold, 'b--', 'LineWidth', 1.4, ...
      'Label', sprintf('MLI skin = %.1f C  (Q_{out} = %.1f W)', T_skin_cold, Q_leak_out), ...
      'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top', ...
      'HandleVisibility', 'off');

xlabel('Component', 'FontSize', 12);
ylabel('Temperature [C]', 'FontSize', 12);
title(sprintf('COLD CASE — Component Temperatures vs Min Limits  (-%d C margin)', ...
    margin), 'FontSize', 12, 'FontWeight', 'bold');
xticks(x); xticklabels(short); xtickangle(32);

h4 = bar(NaN, NaN, 'FaceColor', [0.22 0.63 0.30], 'EdgeColor', 'none');
h5 = bar(NaN, NaN, 'FaceColor', [0.85 0.20 0.15], 'EdgeColor', 'none');
h6 = plot(NaN, NaN, 'b^', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
legend([h4 h5 h6], {'Pass', 'Fail -> heater needed', 'Min limit (op or survival)'}, ...
       'Location', 'northeast', 'FontSize', 9);
% MLI loss vs internal dissipation annotation
Q_int_cold_plot = sum(arrayfun(@(c) c.Q_cold_W * ~c.decoupled, comps));
surplus = Q_int_cold_plot - Q_leak_out;
if surplus >= 0
    bal_str = sprintf('Internal dissipation (%.0f W) covers MLI loss (%.0f W) — surplus %.0f W', ...
        Q_int_cold_plot, Q_leak_out, surplus);
    bg_col = [0.85 1.0 0.85];
else
    bal_str = sprintf('DEFICIT: Internal dissipation (%.0f W) < MLI loss (%.0f W) by %.0f W', ...
        Q_int_cold_plot, Q_leak_out, abs(surplus));
    bg_col = [1.0 0.85 0.85];
end
annotation('textbox', [0.13 0.01 0.75 0.07], 'String', bal_str, ...
    'BackgroundColor', bg_col, 'FontSize', 9, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'EdgeColor', [0.4 0.4 0.4]);
%% ---- Figure 3: Radiator Requirements -----------------------------------
hot_idx = find(A_rad > 0);
n_rad   = length(hot_idx);

fig3 = figure('Name','Fig 3 — Radiator Requirements','NumberTitle','off');

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
cold_idx = find(Q_heat_W > 0);
n_heat   = length(cold_idx);

% fig4 = figure('Name','Fig 4 — Heater Requirements','NumberTitle','off');
% 
% if n_heat == 0
%     text(0.5, 0.5, 'No component heaters required', ...
%          'Units','normalized','HorizontalAlignment','center','FontSize',14);
%     axis off;
% else

%     hold on; grid on;
%     heat_names = short(cold_idx);
%     heat_power = Q_heat(cold_idx);
% 
%     bar(1:n_heat, heat_power, 0.6, 'FaceColor', [0.20 0.55 0.85], 'EdgeColor', 'none');
% 
%     for k = 1:n_heat
%         text(k, heat_power(k)*1.04, sprintf('%.1f W', heat_power(k)), ...
%              'HorizontalAlignment','center','FontSize',9);
%         % Setpoint label inside bar
%         text(k, heat_power(k)*0.45, sprintf('setpt\n%.0f C', T_setpt(cold_idx(k))), ...
%              'HorizontalAlignment','center','FontSize',8,'Color','white','FontWeight','bold');
%     end
% 
%     % Total battery load annotation (includes Q_leak_out)
%     yline(Q_heat_total, 'r--', 'LineWidth', 1.8, ...
%           'Label', sprintf('Battery load = %.1f W (incl. %.1f W MLI loss)', ...
%           Q_heat_total, Q_leak_out), ...
%           'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');
% 
%     xticks(1:n_heat); xticklabels(heat_names); xtickangle(28);
%     ylabel('Required Heater Power [W]', 'FontSize', 12);
%     title(sprintf('COLD CASE — Heater Power per Component  (20%% margin incl.)\nTotal battery load in eclipse: %.1f W', ...
%         Q_heat_total), 'FontSize', 11, 'FontWeight', 'bold');
% end
% 
% end

fig4 = figure('Name','Fig 4 — Heater Requirements','NumberTitle','off');

% Compute interior energy balance for annotation
Q_int_cold_plot = sum(arrayfun(@(c) c.Q_cold_W * ~c.decoupled, comps));
Q_surplus_plot  = Q_int_cold_plot + Q_heat_total - Q_leak_out;
has_deficit     = Q_surplus_plot < 0;

if n_heat == 0 && ~has_deficit
    % No component heaters and interior is self-sufficient
    text(0.5, 0.60, 'No component heaters required', ...
         'Units','normalized','HorizontalAlignment','center', ...
         'FontSize',14,'FontWeight','bold','Color',[0.22 0.63 0.30]);
    text(0.5, 0.45, sprintf('Interior dissipation (%.1f W) covers MLI loss (%.1f W)', ...
         Q_int_cold_plot, Q_leak_out), ...
         'Units','normalized','HorizontalAlignment','center','FontSize',11);
    text(0.5, 0.35, sprintf('Surplus = %.1f W  —  no bulk heater needed', ...
         Q_surplus_plot), ...
         'Units','normalized','HorizontalAlignment','center','FontSize',11);
    axis off;

elseif n_heat == 0 && has_deficit
    % No component heaters but interior cannot sustain 20 C on its own
    text(0.5, 0.72, 'No individual component heaters required', ...
         'Units','normalized','HorizontalAlignment','center', ...
         'FontSize',13,'FontWeight','bold','Color',[0.22 0.63 0.30]);
    text(0.5, 0.57, 'BUT — bulk interior heater required:', ...
         'Units','normalized','HorizontalAlignment','center', ...
         'FontSize',13,'FontWeight','bold','Color',[0.85 0.20 0.15]);
    text(0.5, 0.44, sprintf('Internal dissipation (%.1f W) < MLI heat loss (%.1f W)', ...
         Q_int_cold_plot, Q_leak_out), ...
         'Units','normalized','HorizontalAlignment','center','FontSize',11);
    text(0.5, 0.33, sprintf('Deficit = %.1f W  —  bulk interior heater must supply this', ...
         abs(Q_surplus_plot)), ...
         'Units','normalized','HorizontalAlignment','center', ...
         'FontSize',12,'FontWeight','bold','Color',[0.85 0.20 0.15]);
    text(0.5, 0.20, sprintf('Without this heater, T_{int} drops below 20 C,'), ...
         'Units','normalized','HorizontalAlignment','center','FontSize',10, ...
         'Color',[0.4 0.4 0.4]);
    text(0.5, 0.12, 'invalidating all component temperature predictions.', ...
         'Units','normalized','HorizontalAlignment','center','FontSize',10, ...
         'Color',[0.4 0.4 0.4]);
    axis off;

else
    % Component heaters exist — show bar chart plus deficit annotation if needed
    hold on; grid on;
    heat_names = short(cold_idx);
    heat_power = Q_heat_W(cold_idx);

    bar(1:n_heat, heat_power, 0.6, 'FaceColor', [0.20 0.55 0.85], 'EdgeColor', 'none');

    for k = 1:n_heat
        text(k, heat_power(k)*1.04, sprintf('%.1f W', heat_power(k)), ...
             'HorizontalAlignment','center','FontSize',9);
        text(k, heat_power(k)*0.45, sprintf('setpt\n%.0f C', T_setpt(cold_idx(k))), ...
             'HorizontalAlignment','center','FontSize',8, ...
             'Color','white','FontWeight','bold');
    end

    yline(Q_heat_total, 'r--', 'LineWidth', 1.8, ...
          'Label', sprintf('Component heaters = %.1f W', Q_heat_total), ...
          'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom');

    xticks(1:n_heat); xticklabels(heat_names); xtickangle(28);
    ylabel('Required Heater Power [W]', 'FontSize', 12);
    title(sprintf('COLD CASE — Heater Power per Component  (20%% margin incl.)\nTotal battery load in eclipse: %.1f W', ...
        Q_heat_total), 'FontSize', 11, 'FontWeight', 'bold');

    % Deficit/surplus annotation at bottom of bar chart
    if has_deficit
        def_str = sprintf('DEFICIT: Internal dissipation (%.0f W) + component heaters (%.0f W) < MLI loss (%.0f W) by %.0f W — bulk interior heater needed', ...
            Q_int_cold_plot, Q_heat_total, Q_leak_out, abs(Q_surplus_plot));
        bg_col = [1.0 0.85 0.85];
    else
        def_str = sprintf('OK: Internal dissipation (%.0f W) + component heaters (%.0f W) covers MLI loss (%.0f W) — surplus %.0f W', ...
            Q_int_cold_plot, Q_heat_total, Q_leak_out, Q_surplus_plot);
        bg_col = [0.85 1.0 0.85];
    end
    annotation('textbox', [0.10 0.01 0.85 0.07], 'String', def_str, ...
        'BackgroundColor', bg_col, 'FontSize', 9, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'EdgeColor', [0.4 0.4 0.4]);
end