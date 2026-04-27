function Final_Plots_Version_4_26( ...
        comps, T_hot, T_cold, cold_limit, hot_viol, cold_viol, nd, T_bus_max_C, ...
        A_rad_total, A_rad_avail, Q_heat_W, A_heat, T_setpt, Q_heat_total, ...
        margin, T_skin_hot, T_skin_cold, Q_leak_in, Q_leak_out, ...
        Q_int_cold, Q_batt_load, Q_rad_cold_loss, Q_panel_drain, Q_rad_hot, ...
        T_int_hot_C, T_int_cold_C, idx_bind)
%  FINAL_4_21_PLOTS  Four thermal analysis figures
%
%  Fig 1  HOT CASE   — component temperatures vs T_op_max limits
%  Fig 2  COLD CASE  — component temperatures vs cold limits
%  Fig 3  RADIATOR   — per-component area attribution summing to total required,
%                      with available area as a horizontal limit line
%  Fig 4  HEATERS    — per-heater-needing-component power + total

% ---- Setup: body components only ----------------------------------------
idx_s = find(nd);  Ns = numel(idx_s);  xs = 1:Ns;

names_s = cell(1,Ns);
for k = 1:Ns
    s = comps(idx_s(k)).name;
    s = strrep(s,' (NEXT-C)',''); s = strrep(s,' [CONFIRM]','');
    s = strrep(s,' (2x)',' x2'); s = strrep(s,' (16x)',' x16');
    s = strrep(s,' (4x)',' x4');
    names_s{k} = strtrim(s);
end

T_hot_s    = T_hot(idx_s);
T_cold_s   = T_cold(idx_s);
cold_lim_s = cold_limit(idx_s);
cold_viol_s= cold_viol(idx_s);
Tmax_s     = [comps(idx_s).T_op_max_C];
Q_heat_s   = Q_heat_W(idx_s);
Q_hot_s    = [comps(idx_s).Q_hot_W];
margin_s   = Tmax_s - T_hot_s;

bind_local = find(idx_s == idx_bind, 1);

clr_ok   = [0.25 0.60 0.28];
clr_warn = [0.93 0.70 0.08];
clr_fail = [0.80 0.18 0.14];

%% FIG 1 — HOT CASE --------------------------------------------------------
figure('Name','Fig 1 — Hot Case','NumberTitle','off', ...
       'Color','w','Position',[60 510 900 460]);
hold on; grid on;

for k = 1:Ns
    % Inside your Figure 1 loop...
    if contains(comps(idx_s(k)).name, 'Ion Thruster Head')
        % Place a text note if the value is off-chart
        text(k, 140, sprintf('Actual: ~%.0f°C', 300 ), ...
        'HorizontalAlignment','center','Color','r');
    end
    if     margin_s(k) < 0;   clr = clr_fail;
    elseif margin_s(k) <= 10; clr = clr_warn;
    else;                      clr = clr_ok;
    end
    bar(k, T_hot_s(k), 0.55, 'FaceColor',clr,'EdgeColor','k','LineWidth',0.5);
end
plot(xs, Tmax_s, 'rv','MarkerSize',8,'MarkerFaceColor','r','LineWidth',1.5);
ylim([0 150]); % Focus the graph on the sensitive electronics
% yline(T_int_hot_C,'b-','LineWidth',1.8);
    % 'Label',sprintf('Thermal Bus = %.1f°C',T_int_hot_C), ...
    % 'LabelHorizontalAlignment','left','FontSize',8);
% if ~isempty(bind_local)
%     text(bind_local, T_hot_s(bind_local)+2,'★','HorizontalAlignment','center', ...
%         'FontSize',12,'FontWeight','bold','Color','k');
% end
legend([bar(NaN,NaN,'FaceColor',clr_ok,'EdgeColor','none'), ...
        bar(NaN,NaN,'FaceColor',clr_warn,'EdgeColor','none'), ...
        plot(NaN,NaN,'rv','MarkerFaceColor','r')], ...
    {'OK (>10°C margin)','Close (0–10°C)','T_{op,max} limit'}, ...
    'Location','best','FontSize',10,'Box','off');
set(gca,'XTick',xs,'XTickLabel',names_s,'XTickLabelRotation',35,'FontSize',8);
ylabel('Temperature [°C]');
title(sprintf('HOT CASE Temperatures'),'FontWeight','bold','FontSize',11);

%% FIG 2 — COLD CASE -------------------------------------------------------
figure('Name','Fig 2 — Cold Case','NumberTitle','off', ...
       'Color','w','Position',[60 30 900 460]);
hold on; grid on;

for k = 1:Ns
    clr = ternary(cold_viol_s(k), clr_fail, [0.20 0.48 0.78]);
    bar(k, T_cold_s(k), 0.55, 'FaceColor',clr,'EdgeColor','k','LineWidth',0.5);
end
plot(xs, cold_lim_s, 'b^','MarkerSize',8,'MarkerFaceColor','b','LineWidth',1.5);
% yline(T_int_cold_C,'b-','LineWidth',1.8);
    % 'Label',sprintf('Bus = %.1f°C',T_int_cold_C), ...
    % 'LabelHorizontalAlignment','left','LabelVerticalAlignment','top','FontSize',8);
legend([bar(NaN,NaN,'FaceColor',[0.20 0.48 0.78],'EdgeColor','none'), ...
        bar(NaN,NaN,'FaceColor',clr_fail,'EdgeColor','none'), ...
        plot(NaN,NaN,'b^','MarkerFaceColor','b')], ...
    {'OK','Needs heater','Cold limit'}, ...
    'Location','best','FontSize',10,'Box','off');
set(gca,'XTick',xs,'XTickLabel',names_s,'XTickLabelRotation',35,'FontSize',8);
ylabel('Temperature [°C]');
n_fail = sum(cold_viol_s);
title(sprintf('COLD CASE Temperatures  (%d of %d Components Need Heaters)', ...
     n_fail, Ns),'FontWeight','bold','FontSize',11);

%% FIG 3 — RADIATOR (per-component attribution) ----------------------------
%
%  Each component's share of total radiator area is proportional to its
%  fraction of Q_rad_hot.  Q_rad_hot = Q_int_hot + Q_panel_boom + Q_MLI_leak.
%  Components contribute Q_int_hot; the remainder (panel + MLI) is lumped
%  into a single "Thermal overhead" bar.
%
%  Available area (A_rad_avail) shown as a red dashed limit line.
%  A_rad_avail = (24 - 4) * (1 - 0.25) = 15 m^2 for a 2x2x2 m cube.
%
figure('Name','Fig 3 — Radiator Area by Component','NumberTitle','off', ...
       'Color','w','Position',[980 510 900 460]);
hold on; grid on;

if isinf(A_rad_total) || A_rad_total == 0
    axis off;
    if isinf(A_rad_total)
        text(0.5,0.55,'RADIATOR INFEASIBLE','Units','normalized', ...
            'HorizontalAlignment','center','FontSize',16,'FontWeight','bold','Color',clr_fail);
        text(0.5,0.38,sprintf('T_{int,hot} ≤ T_{sink} — increase G_{W/K} for: %s', ...
            comps(idx_bind).name),'Units','normalized','HorizontalAlignment','center','FontSize',11);
    else
        text(0.5,0.55,'No radiator required','Units','normalized', ...
            'HorizontalAlignment','center','FontSize',16,'FontWeight','bold','Color',clr_ok);
    end
else
    % Per-component area proportional to Q_hot contribution to Q_rad_hot
    A_comp_s = zeros(1, Ns);
    for k = 1:Ns
        A_comp_s(k) = A_rad_total * (Q_hot_s(k) / Q_rad_hot);
    end
    Q_overhead    = Q_rad_hot - sum(Q_hot_s);   % panel boom + MLI leak-in
    A_overhead    = A_rad_total * (Q_overhead / Q_rad_hot);

    % All bars: components + overhead
    all_labels = [names_s, {sprintf('Overhead(MLI+ solar panel boom)')}];
    all_areas  = [A_comp_s, A_overhead];
    Nb = length(all_areas);
    xb = 1:Nb;

    for k = 1:Nb
        if k == bind_local
            clr3 = clr_fail;      % binding component — red
        elseif k == Nb
            clr3 = [0.60 0.60 0.60];  % overhead — grey
        else
            clr3 = [0.30 0.58 0.78];  % other components — blue
        end
        bar(k, all_areas(k), 0.55, 'FaceColor',clr3,'EdgeColor','k','LineWidth',0.5);
        if all_areas(k) > 0.001
            text(k, all_areas(k) + A_rad_total*0.015, ...
                sprintf('%.3f m²', all_areas(k)), ...
                'HorizontalAlignment','center','FontSize',8);
        end
    end

    % Total required line
    yline(A_rad_total,'r-','LineWidth',2.0, ...
        'Label',sprintf('Total required = %.3f m²', A_rad_total), ...
        'LabelHorizontalAlignment','right','LabelVerticalAlignment','top','FontSize',9);
    % 
    % % Available area limit line
    % yline(A_rad_avail,'k--','LineWidth',1.8, ...
    %     'Label',sprintf('Available = %.2f m²', A_rad_avail), ...
    %     'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',9);

    % Binding marker
    if ~isempty(bind_local)
        text(bind_local, all_areas(bind_local) + A_rad_total*0.06, '★ Most Restrictive', ...
            'HorizontalAlignment','center','FontSize',8,'FontWeight','bold','Color',clr_fail);
    end

    set(gca,'XTick',xb,'XTickLabel',all_labels,'XTickLabelRotation',35,'FontSize',8);
    ylabel('Radiator Area Attribution [m²]','FontSize',11);

    % Status
    if A_rad_total <= A_rad_avail
        % status = sprintf('SUFFICIENT — %.0f%% of available area used', ...
        %     100*A_rad_total/A_rad_avail);
        tcol = clr_ok;
    else
        % status = sprintf('INSUFFICIENT — need +%.3f m² beyond limit', A_rad_total-A_rad_avail);
        tcol = clr_fail;
    end
    title(sprintf('Radiator Area by Component  (total = %.3f m²,  rejects %.0f W at bus = %.1f°C)\n%s', ...
        A_rad_total, Q_rad_hot, T_int_hot_C), ...
        'FontWeight','bold','FontSize',14, 'Color', tcol);
end

%% FIG 4 — HEATERS ---------------------------------------------------------
need_htr = find(Q_heat_s > 0);
n_heat   = numel(need_htr);

figure('Name','Fig 4 — Heater Power','NumberTitle','off', ...
       'Color','w','Position',[980 30 900 460]);

if n_heat == 0
    axis off;
    text(0.5,0.65,'NO HEATERS REQUIRED','Units','normalized', ...
        'HorizontalAlignment','center','FontSize',16,'FontWeight','bold','Color',clr_ok);
    text(0.5,0.45,sprintf('Bus = %.1f°C in eclipse\nAll components above cold limit', ...
        T_int_cold_C),'Units','normalized','HorizontalAlignment','center','FontSize',11);
    text(0.5,0.25,sprintf('Eclipse battery load: %.1f W\n(MLI leak %.1f W  +  panel drain %.1f W)', ...
        Q_batt_load, Q_leak_out, Q_panel_drain), ...
        'Units','normalized','HorizontalAlignment','center','FontSize',10,'Color',[0.4 0.4 0.4]);
else
    hold on; grid on;
    htr_power = Q_heat_s(need_htr);
    xh = 1:n_heat;
    bar(xh, htr_power, 0.55, 'FaceColor',[0.88 0.50 0.10],'EdgeColor','k','LineWidth',0.8);
    for k = 1:n_heat
        text(k, htr_power(k) + max(htr_power)*0.03, sprintf('%.1f W',htr_power(k)), ...
            'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
    end
    yline(Q_heat_total,'r--','LineWidth',2.0, ...
        'Label',sprintf('Total: %.1f W',Q_heat_total), ...
        'LabelHorizontalAlignment','right','FontSize',9);
    set(gca,'XTick',xh,'XTickLabel',names_s(need_htr),'XTickLabelRotation',35,'FontSize',8);
    ylabel('Heater Power [W]  (20% margin incl.)','FontSize',11);
    title(sprintf('Component Heaters  (%.0f W total  |  eclipse battery load: %.0f W)', ...
        Q_heat_total, Q_batt_load),'FontWeight','bold','FontSize',11);
end

end % function

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end

