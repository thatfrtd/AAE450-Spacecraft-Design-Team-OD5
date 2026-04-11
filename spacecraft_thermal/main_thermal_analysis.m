%% =========================================================================
%  SPACECRAFT THERMAL ANALYSIS — MAIN ENTRY POINT
%  =========================================================================
%  Project    : [YOUR SPACECRAFT NAME]               % TODO: fill in
%  Orbit      : LEO, 1000–1200 km altitude
%  Author     : [YOUR NAME]                          % TODO: fill in
%  Date       : 2025
%  Version    : 1.0
%
%  DESCRIPTION:
%    Full transient thermal model for a LEO spacecraft. Computes:
%      1. External heat fluxes (solar, albedo, Earth IR)
%      2. Internal heat dissipation from all components
%      3. Node temperatures via lumped-capacitance ODE integration
%      4. Worst-case hot and cold scenarios
%      5. Radiator sizing
%      6. Heater sizing
%
%  UNITS: SI throughout (W, m, m^2, K, J, kg)
%
%  USAGE:
%    Run this script. Results are saved to 'thermal_results.mat' and
%    figures are exported to /figures/
%
%  DEPENDENCIES:
%    spacecraft_geometry.m     — geometry, areas, view factors
%    orbit_model.m             — orbital mechanics, beta angle, eclipse
%    external_heat_fluxes.m   — solar, albedo, Earth IR
%    component_heat_map.m     — internal dissipation per component
%    thermal_network.m        — conductance matrix, node definitions
%    radiator_sizing.m        — required radiator area calculation
%    heater_sizing.m          — required heater power calculation
%    run_transient.m          — ODE45 integration loop
%    plot_results.m           — all figures for presentation
% =========================================================================

clear; clc; close all;
addpath(genpath('.'));

fprintf('====================================================\n');
fprintf('  SPACECRAFT THERMAL ANALYSIS — STARTING\n');
fprintf('====================================================\n\n');

%% -------------------------------------------------------------------------
%  SECTION 1: MISSION & ORBIT PARAMETERS
% --------------------------------------------------------------------------
orbit = struct();
orbit.altitude_km   = 1100;          % km — mid of 1000–1200 range
orbit.altitude_m    = orbit.altitude_km * 1e3;
orbit.R_earth       = 6371e3;        % m — mean Earth radius
orbit.R_orbit       = orbit.R_earth + orbit.altitude_m;
orbit.mu_earth      = 3.986e14;      % m^3/s^2 — Earth gravitational param
orbit.period        = 2*pi * sqrt(orbit.R_orbit^3 / orbit.mu_earth); % s
orbit.eclipse_frac  = 0.40;          % 40% eclipse fraction (given)
orbit.eclipse_dur   = orbit.eclipse_frac * orbit.period;             % s
orbit.sunlit_dur    = (1 - orbit.eclipse_frac) * orbit.period;       % s

fprintf('Orbital period:   %.1f min\n', orbit.period/60);
fprintf('Eclipse duration: %.1f min\n', orbit.eclipse_dur/60);
fprintf('Sunlit duration:  %.1f min\n\n', orbit.sunlit_dur/60);

%% -------------------------------------------------------------------------
%  SECTION 2: SPACECRAFT GEOMETRY
% --------------------------------------------------------------------------
fprintf('Loading spacecraft geometry...\n');
geom = spacecraft_geometry();

%% -------------------------------------------------------------------------
%  SECTION 3: THERMAL ENVIRONMENT
% --------------------------------------------------------------------------
fprintf('Computing thermal environment...\n');
env = struct();
env.S_solar         = 1361;          % W/m^2 — solar constant (1 AU)
env.albedo_earth    = 0.30;          % Earth albedo (visible)   % TODO: verify
env.q_earth_IR      = 237;           % W/m^2 — Earth outgoing LW radiation % TODO: verify for altitude
env.T_space         = 2.7;           % K — deep space background

%% -------------------------------------------------------------------------
%  SECTION 4: SURFACE OPTICAL PROPERTIES (Aluminum + coatings)
% --------------------------------------------------------------------------
%  NOTE: These are placeholder values for bare Al + typical coatings.
%  TODO: Replace with your actual surface finish / MLI spec from materials team.
surf = struct();
surf.alpha_solar_bare   = 0.15;   % absorptivity — bare anodized Al
surf.eps_IR_bare        = 0.05;   % emissivity — bare polished Al
surf.alpha_solar_paint  = 0.92;   % absorptivity — white paint (cold case)  % TODO: confirm coating
surf.eps_IR_paint       = 0.88;   % emissivity — white paint
surf.alpha_solar_MLI    = 0.05;   % MLI outer layer (aluminized Kapton)
surf.eps_IR_MLI         = 0.05;   % MLI outer emissivity

% Default body surface: assume anodized Al with partial MLI
% TODO: Assign per-face coating once final surface finish is decided
surf.alpha_body  = 0.20;   % effective body absorptivity        % TODO: update
surf.eps_body    = 0.80;   % effective body emissivity (incl. paint patches) % TODO: update
surf.alpha_panel = 0.80;   % solar panel front (GaAs cells)
surf.eps_panel   = 0.85;   % solar panel front emissivity

%% -------------------------------------------------------------------------
%  SECTION 5: COMPONENT DEFINITIONS & POWER BUDGET
% --------------------------------------------------------------------------
fprintf('Building component thermal map...\n');
[components, power] = component_heat_map();

%% -------------------------------------------------------------------------
%  SECTION 6: DEFINE WORST-CASE SCENARIOS
% --------------------------------------------------------------------------
% HOT CASE:  max solar flux, max albedo, max internal power, zero eclipse
%            beta angle ~ 90 deg (continuous sun), summer solstice geometry
% COLD CASE: min solar flux, min albedo, max eclipse, min internal power
%            (e.g. safe-mode — only survival heaters + computer on)

scenarios = struct();

% --- HOT CASE ---
scenarios.hot.beta_deg       = 75;    % deg — near-continuous illumination  % TODO: confirm max beta
scenarios.hot.solar_flux     = 1414;  % W/m^2 — perihelion (Earth closest to Sun, ~Jan)
scenarios.hot.albedo         = 0.35;  % upper bound                         % TODO: verify
scenarios.hot.eclipse_frac   = 0.05;  % near-zero eclipse at high beta
scenarios.hot.power_fracs    = ones(1, length(components)); % all components at full power

% --- COLD CASE ---
scenarios.cold.beta_deg      = 0;     % deg — maximum eclipse geometry
scenarios.cold.solar_flux    = 1322;  % W/m^2 — aphelion                   % TODO: verify
scenarios.cold.albedo        = 0.25;  % lower bound
scenarios.cold.eclipse_frac  = 0.40;  % maximum eclipse
scenarios.cold.power_fracs   = zeros(1, length(components)); % safe mode — all off except...
scenarios.cold.power_fracs(power.idx_flight_computer) = 1.0; % flight computer always on
scenarios.cold.power_fracs(power.idx_imu)             = 1.0; % IMU always on

%% -------------------------------------------------------------------------
%  SECTION 7: THERMAL NETWORK SETUP
% --------------------------------------------------------------------------
fprintf('Building thermal network...\n');
[nodes, C_matrix, G_matrix] = thermal_network(geom, components, surf);

%% -------------------------------------------------------------------------
%  SECTION 8: TRANSIENT SIMULATION — HOT CASE
% --------------------------------------------------------------------------
fprintf('\nRunning HOT CASE transient simulation...\n');
results_hot = run_transient(nodes, C_matrix, G_matrix, ...
                             geom, surf, orbit, env, ...
                             scenarios.hot, components, power);
fprintf('  Peak temperature:  %.1f C (node: %s)\n', ...
    max(results_hot.T_peak_C), results_hot.hottest_node);

%% -------------------------------------------------------------------------
%  SECTION 9: TRANSIENT SIMULATION — COLD CASE
% --------------------------------------------------------------------------
fprintf('Running COLD CASE transient simulation...\n');
results_cold = run_transient(nodes, C_matrix, G_matrix, ...
                              geom, surf, orbit, env, ...
                              scenarios.cold, components, power);
fprintf('  Min temperature:   %.1f C (node: %s)\n', ...
    min(results_cold.T_min_C), results_cold.coldest_node);

%% -------------------------------------------------------------------------
%  SECTION 10: RADIATOR SIZING (Hot Case Driven)
% --------------------------------------------------------------------------
fprintf('\nSizing radiators for hot case...\n');
radiator = radiator_sizing(results_hot, components, surf, env);

fprintf('  Required radiator area: %.3f m^2\n', radiator.area_required_m2);
fprintf('  Available body area:    %.3f m^2\n', radiator.area_available_m2);
if radiator.body_sufficient
    fprintf('  >> Body-mounted radiators SUFFICIENT\n');
else
    fprintf('  >> WARNING: Deployable/dedicated radiators may be required\n');
end

%% -------------------------------------------------------------------------
%  SECTION 11: HEATER SIZING (Cold Case Driven)
% --------------------------------------------------------------------------
fprintf('\nSizing heaters for cold case...\n');
heater = heater_sizing(results_cold, components);

fprintf('  Total heater power required: %.1f W\n', heater.total_power_W);
fprintf('  Components needing heaters:\n');
for i = 1:length(heater.components_needing_heat)
    fprintf('    - %s: %.1f W\n', ...
        heater.components_needing_heat{i}, heater.power_per_component_W(i));
end

%% -------------------------------------------------------------------------
%  SECTION 12: CONTROLLED SIMULATION — HOT CASE WITH RADIATORS
% --------------------------------------------------------------------------
fprintf('\nRunning HOT CASE with radiators applied...\n');
results_hot_ctrl = run_transient_controlled(nodes, C_matrix, G_matrix, ...
                                             geom, surf, orbit, env, ...
                                             scenarios.hot, components, power, ...
                                             heater, radiator, 'hot');
fprintf('  Controlled peak: %.1f C (was %.1f C)\n', ...
    max(results_hot_ctrl.T_peak_C), max(results_hot.T_peak_C));

%% -------------------------------------------------------------------------
%  SECTION 13: CONTROLLED SIMULATION — COLD CASE WITH HEATERS
% --------------------------------------------------------------------------
fprintf('Running COLD CASE with heaters applied...\n');
results_cold_ctrl = run_transient_controlled(nodes, C_matrix, G_matrix, ...
                                              geom, surf, orbit, env, ...
                                              scenarios.cold, components, power, ...
                                              heater, radiator, 'cold');
fprintf('  Controlled min:  %.1f C (was %.1f C)\n', ...
    min(results_cold_ctrl.T_min_C), min(results_cold.T_min_C));

%% -------------------------------------------------------------------------
%  SECTION 14: RESULTS SUMMARY TABLE
% --------------------------------------------------------------------------
fprintf('\n====================================================\n');
fprintf('  THERMAL SUMMARY — ALL COMPONENTS\n');
fprintf('====================================================\n');
print_summary_table(results_hot, results_cold, components, heater, radiator);

%% -------------------------------------------------------------------------
%  SECTION 15: SAVE & PLOT
% --------------------------------------------------------------------------
fprintf('\nSaving results...\n');
save('thermal_results.mat', 'results_hot', 'results_cold', ...
     'results_hot_ctrl', 'results_cold_ctrl', ...
     'radiator', 'heater', 'components', 'orbit', 'geom', 'surf');

fprintf('Generating figures...\n');
plot_results(results_hot, results_cold, radiator, heater, components, orbit, ...
             results_hot_ctrl, results_cold_ctrl);

fprintf('\n====================================================\n');
fprintf('  ANALYSIS COMPLETE\n');
fprintf('====================================================\n');
