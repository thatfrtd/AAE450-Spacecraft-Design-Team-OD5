function print_summary_table(results_hot, results_cold, components, heater, radiator)
% PRINT_SUMMARY_TABLE  Print formatted thermal summary to console
%
%  Displays a table of all components with:
%    - Peak temp (hot case) vs max operating limit
%    - Min temp (cold case) vs min operating limit
%    - Hot margin, cold margin
%    - Heater required flag
%    - Overall pass/fail

N = length(components);
node_names = {results_hot.nodes(1:N).name};

fprintf('\n');
fprintf('%-25s | %8s | %8s | %8s | %8s | %7s | %7s | %6s\n', ...
    'Component', 'T_peak°C', 'Tmax_op', 'T_min°C', 'Tmin_op', ...
    'HotMrgn', 'CldMrgn', 'Status');
fprintf('%s\n', repmat('-', 1, 95));

for i = 1:N
    c      = components(i);
    T_pk   = results_hot.T_peak_C(i);
    T_mn   = results_cold.T_min_C(i);
    T_mxop = c.T_op_max_C;
    T_mnop = c.T_op_min_C;
    mh     = T_mxop - T_pk;   % hot margin (positive = OK)
    mc     = T_mn - T_mnop;   % cold margin (positive = OK)

    pass = mh >= 0 && mc >= 0;
    if pass
        status = 'PASS';
    else
        status = 'FAIL ***';
    end

    fprintf('%-25s | %8.1f | %8.0f | %8.1f | %8.0f | %7.1f | %7.1f | %s\n', ...
        strrep(c.name,'_',' '), T_pk, T_mxop, T_mn, T_mnop, mh, mc, status);
end

fprintf('%s\n', repmat('-', 1, 95));
if radiator.body_sufficient
    rad_status = 'BODY SUFFICIENT';
else
    rad_status = 'NEED MORE AREA';
end
fprintf('\nRadiator:  Required %.2f m^2  |  Available %.2f m^2  |  %s\n', ...
    radiator.area_required_m2, radiator.area_available_m2, rad_status);
fprintf('Heaters:   Total budget = %.1f W  (%d components)\n', ...
    heater.total_power_W, heater.n_heaters);
fprintf('\nNOTE: All values from placeholder parameters — replace with actual data.\n');

end
