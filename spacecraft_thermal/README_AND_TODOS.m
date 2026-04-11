% =========================================================================
%  SPACECRAFT THERMAL ANALYSIS — README & TODO TRACKER
% =========================================================================
%
%  HOW TO RUN:
%    >> main_thermal_analysis
%
%  FILE STRUCTURE:
%    main_thermal_analysis.m  — entry point, top-level script
%    spacecraft_geometry.m    — geometry, areas, view factors
%    component_heat_map.m     — power budget, temp limits, thermal mass
%    external_heat_fluxes.m  — solar / albedo / Earth IR computation
%    thermal_network.m        — node definitions, conductance matrix
%    run_transient.m          — ODE integration (ode15s)
%    radiator_sizing.m        — radiator area calculation
%    heater_sizing.m          — heater power calculation
%    plot_results.m           — all presentation figures
%    print_summary_table.m    — console summary table
%
% =========================================================================
%  OUTSTANDING TODO LIST — Items to Replace with Real Values
% =========================================================================
%
%  FROM STRUCTURES TEAM:
%    [ ] Actual spacecraft mass budget (kg per component)
%    [ ] Wall thickness and material spec (currently Al 6061, 3mm)
%    [ ] Mounting interface contact areas (m^2) per component
%    [ ] Confirm arm geometry on top face (blocks top radiator?)
%    [ ] Confirm RCS thruster blockage fraction on side (currently 30%)
%    [ ] Surface finish / coating spec per face (currently Al + paint)
%
%  FROM PROPULSION TEAM:
%    [ ] NEXT-C ion thruster thermal dissipation (% of 6.9kW becomes heat)
%        Currently assumed: 30% = 2070 W. Get from Aerojet Rocketdyne ICD.
%    [ ] RCS thruster average power dissipation (currently 50W, very rough)
%    [ ] RCS duty cycle in nominal and safe mode
%    [ ] Confirm ion thruster temperature limits (operating and survival)
%
%  FROM POWER / EPS TEAM:
%    [ ] Battery exact mass and thermal capacitance
%    [ ] Battery exact mounting interface area
%    [ ] Available battery energy in eclipse (Wh) — must cover heaters
%    [ ] Solar panel back-surface coating spec (for back-side radiation)
%    [ ] Solar panel thermal conductance to bus (hinge/boom spec)
%
%  FROM ADCS TEAM:
%    [ ] Confirm max/typical beta angle range
%    [ ] Pointing modes (nadir, Sun-pointing, safe mode) and durations
%    [ ] CMG exact mass (currently 4 x 1.5 kg = 6 kg estimated)
%
%  FROM PAYLOAD TEAM:
%    [ ] Camera exact model, mass, operating duty cycle
%    [ ] Camera survival temperature limits
%    [ ] Star tracker exact model and mounting location
%
%  FROM THERMAL VENDOR / HERITAGE DATA:
%    [ ] Thermal interface material (TIM) conductance — currently 2000 W/m^2/K
%    [ ] MLI specification (layers, seams, attachment)
%    [ ] Exact surface coating absorptivity/emissivity (alpha, eps) per face
%        Currently: alpha_body=0.20, eps_body=0.80 (rough anodized+paint)
%
%  MODEL IMPROVEMENTS (future work):
%    [ ] Replace simplified view factors with ESARAD/Thermal Desktop model
%    [ ] Add heater ON/OFF logic inside ODE (currently sized by steady-state)
%    [ ] Add attitude-dependent solar/albedo view factors
%    [ ] Add ion plume heating effect on bottom surface
%    [ ] Add thermal strap/heat pipe conductances once routing is finalized
%    [ ] Validate with ESATAN-TMS or Thermal Desktop for CDR
%    [ ] Add Monte Carlo uncertainty analysis on key parameters
%
% =========================================================================
%  THEORY SUMMARY (for reference)
% =========================================================================
%
%  GOVERNING EQUATION (each node i):
%
%  C_i * dT_i/dt = Q_solar_i + Q_albedo_i + Q_earthIR_i    [environment]
%                + Q_int_i                                   [internal dissipation]
%                + sum_j G_ij*(T_j - T_i)                  [conduction]
%                - sigma * eps_i * A_rad_i * (T_i^4 - T_space^4)  [radiation loss]
%                + Q_heater_i                               [heater, if active]
%
%  RADIATOR SIZING (steady-state hot case):
%    A_rad = Q_reject / [eps * sigma * (T_rad^4 - T_sink^4) - q_absorbed]
%
%  HEATER SIZING (eclipse cold case):
%    P_heater = C_thermal * delta_T / t_eclipse
%
%  KEY CONSTANTS:
%    sigma = 5.6704e-8 W/m^2/K^4
%    S_solar = 1361 W/m^2 (1 AU mean)
%    T_space = 2.7 K
%    R_Earth = 6371 km
%    mu_Earth = 3.986e14 m^3/s^2
%
% =========================================================================
