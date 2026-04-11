function [components, power] = component_heat_map()
% COMPONENT_HEAT_MAP  Define all spacecraft components, power dissipation,
%                     temperature limits, and thermal mass.
%
%  DESCRIPTION:
%    Populates a struct array of components with:
%      - Name and location on spacecraft
%      - Power dissipation (nominal and peak)
%      - Operating temperature limits [T_min, T_max] in °C
%      - Survival temperature limits  [T_surv_min, T_surv_max] in °C
%      - Thermal mass (Cp * m) in J/K
%      - Node index in thermal network
%
%  POWER BUDGET (from provided data):
%    Ion thruster  (NEXT-C) : 6900 W
%    Solar panels  (2x)     : generates 3600 W (not dissipated internally)
%    Battery                :   75 W dissipated
%    Camera                 :   20 W
%    Star tracker           : 7.5 W (midpoint of 5.6–8.9)
%    QH Antennas (2x)       :   10 W total (assumed per pair)
%    CMG (4x)               : 27.5 W total (midpoint of 20–35)
%    Flight computer        :    5 W
%    IMU                    :   12 W
%    RCS thrusters (16x)    :  ~50 W total (estimate)   % TODO: get from propulsion
%
%  NOTE ON EFFICIENCY:
%    All input electrical power is assumed to become waste heat EXCEPT
%    for the solar panels (output goes to bus) and the ion thruster
%    (majority is kinetic — see efficiency note below).
%
%  TODO: Replace all mass estimates with actual mass budget values.
%  TODO: Confirm survival temp limits with component vendors.
%  TODO: Confirm RCS thruster dissipation with propulsion team.
%  TODO: Confirm exact duty cycles for each operating mode.
% =========================================================================

sigma = 5.6704e-8;   % Stefan-Boltzmann (not used here, just for context)

%% ---- Helper: Build a single component struct ----------------------------
%  make_comp(name, location, P_nom_W, P_peak_W, T_op_min, T_op_max,
%             T_surv_min, T_surv_max, mass_kg, Cp_J_kgK, node_idx)

function c = make_comp(name, location, P_nom_W, P_peak_W, ...
                        T_op_min, T_op_max, T_surv_min, T_surv_max, ...
                        mass_kg, Cp_J_kgK, node_idx)
    c.name          = name;
    c.location      = location;
    c.P_nom_W       = P_nom_W;     % W — nominal dissipation
    c.P_peak_W      = P_peak_W;    % W — peak dissipation
    c.T_op_min_C    = T_op_min;    % °C — min operating temp
    c.T_op_max_C    = T_op_max;    % °C — max operating temp
    c.T_surv_min_C  = T_surv_min;  % °C — min survival (heater setpoint reference)
    c.T_surv_max_C  = T_surv_max;  % °C — max survival
    c.mass_kg       = mass_kg;     % kg
    c.Cp_J_kgK      = Cp_J_kgK;   % J/(kg·K)
    c.C_thermal_J_K = mass_kg * Cp_J_kgK;  % J/K — thermal capacitance
    c.node_idx      = node_idx;
end

%% ==========================================================================
%  COMPONENT DEFINITIONS
%  =========================================================================

i = 0;  % node counter

%% 1. NEXT-C Ion Thruster --------------------------------------------------
%  Power:  6900 W input to thruster
%  Thruster efficiency (NEXT-C): ~70% goes to thrust, ~30% becomes heat
%  Thermal dissipation estimate: 0.30 * 6900 = ~2070 W    % TODO: verify with Aerojet Rocketdyne datasheet
%  Location: bottom of spacecraft
%  Temp limits: operating -20 to 300°C (discharge chamber can be very hot)
%  TODO: Get exact thermal map from NEXT-C ICD — thruster has multiple zones
i = i+1;
thruster_efficiency = 0.70;  % ~70% to thrust  % TODO: confirm from NEXT-C spec
P_thruster_heat = 6900 * (1 - thruster_efficiency);

comp(i) = make_comp('NEXT_C_Ion_Thruster', 'bottom_exterior', ...
    P_thruster_heat, P_thruster_heat * 1.1, ...
    -20, 300, -30, 350, ...   % TODO: confirm temp limits from datasheet
    12.0, 500, i);            % TODO: mass from NEXT-C ICD (estimate 12 kg)
                              % TODO: Cp — thruster is complex; use Fe/Al avg

%% 2. Solar Panels (2x, as a single thermal node) -------------------------
%  Input solar: up to ~1414 W/m^2 * 2 * 3.2*9.0 * efficiency
%  Electrical output: 3600 W (given)
%  Heat dissipated = absorbed solar - electrical output
%  Cell efficiency (GaAs triple junction): ~30%              % TODO: confirm
%  Total absorbed (assuming alpha=0.80): 1361 * 2*28.8 * 0.80 = ~62,800 W absorbed?
%  But only ~3600 W extracted as electricity; the rest is re-radiated + conducted
%  Thermal dissipation to structure: estimate 200 W per panel  % TODO: from panel vendor
i = i+1;
comp(i) = make_comp('Solar_Panels', 'wing_exterior', ...
    20, 40, ...   % W — heat into structure; panels radiate most directly  % TODO: update
    28, 140, -40, 150, ...  % op: peak power near 28°C; max 140°C
    162.0, 900, i);   % TODO: mass from panel vendor (estimate 162 kg total both panels)

%% 3. Battery ---------------------------------------------------------------
%  Mass estimate for Li-ion battery in a 0.15m cube:
%  Volume = 0.15^3 = 0.003375 m^3; Li-ion density ~2500 J/kg/K, rho ~2400 kg/m^3
%  Estimated mass: ~12 kg  % TODO: get from EPS team
i = i+1;
comp(i) = make_comp('Battery', 'interior_lower', ...
    75, 100, ...     % W nominal / peak
    10, 25, 0, 45, ...  % tight thermal limits — Li-ion sensitive!
    12.0, 960, i);  % TODO: confirm mass and Cp from EPS/battery datasheet

%% 4. Camera ---------------------------------------------------------------
i = i+1;
comp(i) = make_comp('Camera', 'top_exterior', ...
    20, 20, ...
    0, 40, -10, 50, ...   % TODO: confirm survival from payload team
    1.5, 900, i);          % TODO: mass estimate — depends on camera model

%% 5. Star Tracker ---------------------------------------------------------
i = i+1;
comp(i) = make_comp('Star_Tracker', 'exterior_side', ...
    7.5, 8.9, ...   % midpoint of 5.6–8.9 W
    -30, 60, -40, 70, ...  % TODO: confirm from vendor datasheet
    0.8, 900, i);   % TODO: mass estimate

%% 6. Quadrifilar Helix Antennas (2x, as one node) -------------------------
%  Total power: 10 W (both antennas combined)
i = i+1;
comp(i) = make_comp('QH_Antennas', 'exterior_top', ...
    10, 10, ...
    -40, 80, -50, 90, ...   % TODO: confirm survival from antenna vendor
    0.5, 900, i);

%% 7. CMGs (4x, modeled as one lumped node) --------------------------------
%  Total power: midpoint = 27.5 W
%  CMGs are typically mounted to a rigid structure — model as single node
i = i+1;
comp(i) = make_comp('CMGs', 'interior_central', ...
    27.5, 35, ...
    -20, 70, -30, 85, ...   % TODO: confirm from CMG vendor
    6.0, 500, i);   % TODO: mass (4x CMG units, estimate 1.5 kg each = 6 kg)

%% 8. Flight Computer ------------------------------------------------------
i = i+1;
comp(i) = make_comp('Flight_Computer', 'interior_central', ...
    5, 8, ...
    -20, 70, -30, 85, ...   % TODO: confirm from avionics team
    0.5, 900, i);

%% 9. IMU ------------------------------------------------------------------
i = i+1;
comp(i) = make_comp('IMU', 'interior_central', ...
    12, 12, ...
    -62, 85, -62, 95, ...   % TODO: confirm survival temps
    0.3, 900, i);

%% 10. RCS Thrusters (16x, as one node) -----------------------------------
%  TODO: Get actual power dissipation from propulsion team
%  Estimate: monoprop or biprop RCS; each thruster fires intermittently
%  Rough estimate: 50 W total average dissipation during operation
%  Duty cycle assumed low (attitude corrections)  % TODO: confirm with propulsion
i = i+1;
comp(i) = make_comp('RCS_Thrusters', 'exterior_side', ...
    50, 200, ...    % TODO: wildly uncertain — get from propulsion team
    -40, 120, -50, 150, ...  % TODO: confirm temp limits
    3.0, 500, i);   % TODO: mass estimate for 16 thrusters + feed lines

%% 11. Spacecraft Structure (shell node — lumped) --------------------------
%  Shell node represents conduction to/from all exterior surfaces
i = i+1;
comp(i) = make_comp('Structure_Shell', 'exterior', ...
    0, 0, ...       % no internal dissipation
    -60, 120, -80, 150, ...  % structural limits  % TODO: confirm from structures
    50.0, 896, i);  % TODO: mass from mass budget; Cp for Al 6061

%% ---- Assemble output ---------------------------------------------------
components = comp;

% Index shortcuts for use in main script
power.idx_flight_computer = 8;
power.idx_imu             = 9;
power.idx_battery         = 3;
power.idx_ion_thruster    = 1;

% Total dissipation (nominal)
power.total_nom_W = sum([comp.P_nom_W]);
power.total_peak_W = sum([comp.P_peak_W]);

fprintf('  [Power] Total nominal internal dissipation: %.1f W\n', power.total_nom_W);
fprintf('  [Power] Total peak   internal dissipation: %.1f W\n', power.total_peak_W);

end
