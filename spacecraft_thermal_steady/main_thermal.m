%% =========================================================================
%  SPACECRAFT STEADY-STATE THERMAL ANALYSIS
%  AAE450 — Combined External Heat Flux & Component Thermal Analysis
%  =========================================================================
%
%  METHOD OVERVIEW
%  ---------------
%  STEP 1  (Sections 2-4): EXTERIOR SKIN TEMPERATURE  (teammate model)
%    Solves for exterior MLI skin temperature and net heat leak
%    through the MLI+honeycomb wall, using fzero.
%
%    At steady state on the skin:
%      Q_ext  =  ε·σ·A·(Ts^4 - Tsp^4)  +  A·(Ts - T_int)/R_spec
%
%    HOT:  Q_leak_in  = heat leaking IN  (adds to interior thermal load)
%    COLD: Q_leak_out = heat leaking OUT (must be replaced by heaters)
%
%  STEP 2  (Sections 5-8): COMPONENT TEMPERATURES
%    Interior held at T_int = 20°C (bus requirement).
%    Each component runs at:
%      T_comp = T_int + Q_gen / G_contact  (+/- design margin)
%
%  STEP 3  (Sections 9-11): THERMAL CONTROL SIZING
%    Radiators sized for components exceeding T_op_max.
%    Heaters sized for components below their cold limit,
%    plus budget to replace Q_leak_out through MLI.
%
%  SOURCES
%  -------
%    [S1] SMAD (Space Mission Engineering: The New SMAD, 2026 Ref)
%    [S2] NASA-TP-2002-209914 (NASA Thermal Control Handbook)
%    [S3] NASA GSFC-STD-7000
%    [S4] ECSS-E-ST-31C (Space Engineering: Thermal Control)
%
%  UNITS: SI throughout — W, m, m^2, K, C
% =========================================================================

clear; clc; close all;

sigma = 5.6704e-8;   % [W/m^2/K^4]  Stefan-Boltzmann constant

%% =========================================================================
%  SECTION 1 — SMAD ENVIRONMENTAL PARAMETERS
%  Values provided by teammate from SMAD Tables 11-45A and 11-48A.
%  TODO: Re-verify all values if orbit altitude changes significantly.
% ==========================================================================

% ---- WORST-CASE HOT ----
% Full sunlight, beta=90 deg, max projected area (Y-face, 6 m^2), eclipse=0
hot.S              = 1420;   % [W/m^2]  solar flux                 [S1, Table 11-48A]
hot.q_IR           = 244;    % [W/m^2]  Earth IR flux              [S1, Table 11-45A]
hot.albedo_skin    = 0.30;   % [-]      Earth albedo (teammate used average in skin calc)
hot.F_earth        = 0.15;   % [-]      view factor s/c -> Earth at ~700 km LEO [S2]

% ---- WORST-CASE COLD ----
% Fully eclipsed, beta=0 deg, only Earth IR, entire surface radiates
cold.q_IR          = 218;    % [W/m^2]  Earth IR flux              [S1, Table 11-45A]
cold.F_earth       = 0.15;   % [-]      view factor (same orbit)   [S2]

% ---- MLI OPTICAL PROPERTIES (teammate EOL degraded values) ----
% Aluminized Kapton outer cover, ~4.2 year mission, AO + UV degraded
alpha_MLI = 0.50;   % [-]  absorptivity (degraded from 0.40)  [S2]
eps_MLI   = 0.80;   % [-]  emissivity                         [S2]

%% =========================================================================
%  SECTION 2 — SPACECRAFT GEOMETRY  (cuboid, updated from teammate)
%  The design has changed from cylinder to box.
%  TODO: Confirm final dimensions with structures team.
% ==========================================================================

Lx = 3.0;   % [m]
Ly = 2.0;   % [m]
Lz = 2.0;   % [m]

A_Xface = Ly * Lz;              % = 4.0 m^2
A_Yface = Lx * Lz;              % = 6.0 m^2  <- max sun-facing face (hot case)
A_Zface = Lx * Ly;              % = 6.0 m^2
A_total = 2*(A_Xface + A_Yface + A_Zface);   % = 40 m^2 total exterior surface

% Projected areas for flux calculations
A_hot_face  = A_Yface;   % 6 m^2  max face broadside to sun (beta=90 deg)
A_cold_proj = A_Xface;   % 4 m^2  min projected area toward Earth (beta=0 deg, end-on)
%   NOTE: Cold-case uses A_total for full-surface radiation loss (all faces radiate).

fprintf('=================================================================\n');
fprintf('  GEOMETRY  (Box: %.1f x %.1f x %.1f m)\n', Lx, Ly, Lz);
fprintf('=================================================================\n');
fprintf('  Hot  case sun-facing area   : %.1f m^2\n', A_hot_face);
fprintf('  Total exterior surface      : %.1f m^2\n', A_total);

%% =========================================================================
%  SECTION 3 — MLI + HONEYCOMB WALL RESISTANCE
%  Reproduces teammate model exactly.
%  C_eff = effective thermal conductance of MLI blanket [W/m^2 K]
%  TODO: Verify C_eff values with structures/thermal team once MLI spec finalised.
%  TODO: Confirm honeycomb thickness and effective conductivity.
% ==========================================================================

C_eff_hot  = 0.050;   % [W/m^2 K]  MLI conductance, hot case   [S2]  % TODO: verify
C_eff_cold = 0.030;   % [W/m^2 K]  MLI conductance, cold case  [S2]  % TODO: verify
t_Al       = 0.040;   % [m]        honeycomb panel thickness    % TODO: confirm with structures
k_Al       = 2.5;     % [W/m K]    honeycomb effective k        % TODO: confirm

R_spec_hot  = (1/C_eff_hot)  + (t_Al/k_Al);   % [m^2 K/W]
R_spec_cold = (1/C_eff_cold) + (t_Al/k_Al);   % [m^2 K/W]

T_int   = 293.15;   % [K]  interior bus target = 20 C  % TODO: confirm with EPS/avionics
T_space = 3.0;      % [K]  deep space background

fprintf('  MLI wall resistance (hot)   : %.3f m^2 K/W\n', R_spec_hot);
fprintf('  MLI wall resistance (cold)  : %.3f m^2 K/W\n', R_spec_cold);

%% =========================================================================
%  SECTION 4 — EXTERIOR SKIN TEMPERATURE & HEAT LEAK  (teammate model)
%
%  Energy balance on exterior skin surface:
%    Q_ext = eps*sigma*A*(Ts^4 - Tspace^4) + A*(Ts - T_int)/R_spec
%
%  Solved with fzero. Results validated against teammate output:
%    Hot:  T_skin = 360.70 K (87.55 C),  Q_leak_in   = 20.25 W
%    Cold: T_skin = 160.55 K (-112.60 C), Q_leak_out = 127.24 W
% ==========================================================================

% ---- HOT CASE SKIN ----
Q_ext_hot = (alpha_MLI * hot.S               * A_hot_face) ...
          + (alpha_MLI * hot.S * hot.albedo_skin * hot.F_earth * A_hot_face) ...
          + (eps_MLI   * hot.q_IR * hot.F_earth  * A_hot_face);

func_hot     = @(Ts) Q_ext_hot ...
    - eps_MLI * sigma * A_hot_face * (Ts^4 - T_space^4) ...
    - (Ts - T_int) / (R_spec_hot / A_hot_face);
T_skin_hot_K = fzero(func_hot, 350);
T_skin_hot_C = T_skin_hot_K - 273.15;
Q_leak_in    = (T_skin_hot_K - T_int) / (R_spec_hot / A_hot_face);

% ---- COLD CASE SKIN ----
Q_ext_cold = eps_MLI * cold.q_IR * cold.F_earth * A_total;

func_cold     = @(Ts) Q_ext_cold ...
    - eps_MLI * sigma * A_total * (Ts^4 - T_space^4) ...
    - (Ts - T_int) / (R_spec_cold / A_total);
T_skin_cold_K = fzero(func_cold, 150);
T_skin_cold_C = T_skin_cold_K - 273.15;
Q_leak_out    = (T_int - T_skin_cold_K) / (R_spec_cold / A_total);

fprintf('\n=================================================================\n');
fprintf('  MLI SKIN TEMPERATURES  (validated against teammate output)\n');
fprintf('=================================================================\n');
fprintf('  HOT  skin temperature       : %7.2f K  (%6.2f C)\n', T_skin_hot_K, T_skin_hot_C);
fprintf('  HOT  heat leak INTO s/c     : %7.2f W\n', Q_leak_in);
fprintf('  COLD skin temperature       : %7.2f K  (%6.2f C)\n', T_skin_cold_K, T_skin_cold_C);
fprintf('  COLD heat leak OUT of s/c   : %7.2f W  (heaters must replace this)\n', Q_leak_out);
fprintf('  Reference (teammate output) : 360.70 K / 20.25 W  |  160.55 K / 127.24 W\n');

%% =========================================================================
%  SECTION 5 — SOLAR PANEL THERMAL LOADS
%  TODO: Confirm 3624 W is waste heat to structure (not total waste heat).
%  TODO: Confirm panel dimensions and back-surface emissivity with vendor.
%  TODO: Confirm G_boom_W_K with mechanical team.
% ==========================================================================

panel.width_m    = 2.904;   % [m]   % TODO: confirm EPS
panel.length_m   = 9.607;   % [m]   % TODO: confirm EPS
panel.A_total    = 2 * panel.width_m * panel.length_m;   % 55.8 m^2, both panels
panel.eps_back   = 0.85;    % [-]   back-surface emissivity  % TODO: vendor
panel.G_boom_W_K = 5.0;     % [W/K] deployment boom conductance  % TODO: mech team
panel.VF_earth   = 0.12;    % [-]   back-face view factor to Earth  % TODO: verify
panel.Q_hot_W    = 3624;    % [W]   waste heat to structure in full sun  (team value)

%% =========================================================================
%  SECTION 6 — COMPONENT DEFINITIONS
%  All specs in components.m — edit that file when team data arrives.
%  To add a component: copy one block in components.m, increment i, fill fields.
% ==========================================================================

comps = components();
N = length(comps);

Q_int_hot  = sum(arrayfun(@(c) c.Q_hot_W,  comps));
Q_int_cold = sum(arrayfun(@(c) c.Q_cold_W, comps));

fprintf('\n=================================================================\n');
fprintf('  INTERNAL HEAT LOADS\n');
fprintf('=================================================================\n');
fprintf('  Hot  case (all systems on)  : %7.1f W\n', Q_int_hot);
fprintf('  Cold case (safe mode)       : %7.1f W\n', Q_int_cold);
fprintf('  Cold case MLI heat leak out : %7.1f W\n', Q_leak_out);
fprintf('  Total cold-case heat need   : %7.1f W\n', Q_int_cold + Q_leak_out);

%% =========================================================================
%  SECTION 7 — COMPONENT TEMPERATURES
%
%  T_comp = T_int + Q_gen / G_contact  (+margin for hot, -margin for cold)
%
%  T_int held at 20 C (bus requirement).
%  G_contact [W/K] is thermal conductance component-to-structure.
%  Design margin = 5 C per ECSS-E-ST-31C preliminary design.
%  TODO: Tighten to 3 C at CDR after test validation.
% ==========================================================================

margin_C = 5;   % [C]  % TODO: reduce to 3 C at CDR

T_int_C = T_int - 273.15;   % 20 C
T_comp_hot_C  = zeros(N, 1);
T_comp_cold_C = zeros(N, 1);

for i = 1:N
    if i == 2
        % Solar panels computed from their own radiation balance
        Q_panel_net = alpha_MLI * hot.S * panel.A_total - panel.Q_hot_W;
        T_panel_hot = (Q_panel_net / (panel.eps_back * sigma * panel.A_total))^0.25 - 273.15;
        T_comp_hot_C(i) = T_panel_hot + margin_C;

        f_pan = @(Tp) panel.eps_back * sigma * 2*panel.A_total * Tp^4 ...
                    + panel.G_boom_W_K * (Tp - T_int) ...
                    - panel.eps_back * cold.q_IR * panel.A_total * panel.VF_earth;
        T_panel_cold_K = fzero(f_pan, [10, 400]);
        T_comp_cold_C(i) = T_panel_cold_K - 273.15 - margin_C;
    else
        T_comp_hot_C(i)  = T_int_C + comps(i).Q_hot_W  / comps(i).G_W_K + margin_C;
        T_comp_cold_C(i) = T_int_C + comps(i).Q_cold_W / comps(i).G_W_K - margin_C;
    end
end

%% =========================================================================
%  SECTION 8 — VIOLATION CHECK
% ==========================================================================

hot_viol   = false(N, 1);
cold_viol  = false(N, 1);
cold_limit = zeros(N, 1);

for i = 1:N
    if T_comp_hot_C(i) > comps(i).T_op_max_C
        hot_viol(i) = true;
    end
    if comps(i).always_on
        cold_limit(i) = comps(i).T_op_min_C;
    else
        cold_limit(i) = comps(i).T_surv_min_C;
    end
    if T_comp_cold_C(i) < cold_limit(i)
        cold_viol(i) = true;
    end
end

%% =========================================================================
%  SECTION 9 — RADIATOR SIZING
%  A_rad = Q_gen / [ eps_rad * sigma * (T_max^4 - T_sink^4) ]
%  Radiator: white paint or OSR, eps=0.85, faces away from sun.
%  TODO: Confirm coating with thermal coatings vendor.
%  TODO: Confirm blockage fraction with propulsion/mech team.
% ==========================================================================

eps_rad      = 0.85;    % radiator emissivity  % TODO: confirm
VF_earth_rad = 0.15;    % view factor radiator->Earth
T_sink_K     = (eps_rad * hot.q_IR * VF_earth_rad / (eps_rad * sigma))^0.25;

A_rad_m2 = zeros(N, 1);
for i = 1:N
    if hot_viol(i)
        T_max_K      = comps(i).T_op_max_C + 273.15;
        q_net_W_m2   = eps_rad * sigma * (T_max_K^4 - T_sink_K^4);
        A_rad_m2(i)  = comps(i).Q_hot_W / q_net_W_m2;
    end
end
A_rad_total = sum(A_rad_m2);

frac_blocked    = 0.25;   % fraction of non-sun faces blocked  % TODO: confirm
A_rad_available = (A_total - A_hot_face) * (1 - frac_blocked);

%% =========================================================================
%  SECTION 10 — HEATER SIZING
%  Q_heater = G * (T_limit - T_int) - Q_gen_cold,  plus 20% margin.
%  Total battery load in eclipse = component heaters + Q_leak_out.
%  TODO: Confirm battery can supply total load for full eclipse duration.
% ==========================================================================

heater_density = 1.0;   % [W/cm^2]  % TODO: confirm with Kapton heater vendor
Q_heater_W  = zeros(N, 1);
T_setpt_C   = zeros(N, 1);
A_heat_cm2  = zeros(N, 1);

for i = 1:N
    if cold_viol(i)
        T_target = cold_limit(i) + 273.15;
        Q_h = comps(i).G_W_K * (T_target - T_int) - comps(i).Q_cold_W;
        Q_h = max(Q_h, 0) * 1.20;
        Q_heater_W(i) = Q_h;
        T_setpt_C(i)  = cold_limit(i) + 2;
        A_heat_cm2(i) = Q_h / heater_density;
    end
end

Q_heater_total        = sum(Q_heater_W);
Q_heater_budget_total = Q_heater_total + Q_leak_out;

%% =========================================================================
%  SECTION 11 — CONSOLE REPORT
% ==========================================================================

fprintf('\n=================================================================\n');
fprintf('  HOT CASE ENERGY BUDGET\n');
fprintf('=================================================================\n');
fprintf('  Solar direct (6 m^2 face)   : %8.1f W\n', alpha_MLI*hot.S*A_hot_face);
fprintf('  Albedo (F=%.2f)             : %8.1f W\n', hot.F_earth, alpha_MLI*hot.S*hot.albedo_skin*hot.F_earth*A_hot_face);
fprintf('  Earth IR (F=%.2f)           : %8.1f W\n', hot.F_earth, eps_MLI*hot.q_IR*hot.F_earth*A_hot_face);
fprintf('  Total external on skin      : %8.1f W\n', Q_ext_hot);
fprintf('  Skin temperature            : %8.2f K  (%5.2f C)\n', T_skin_hot_K, T_skin_hot_C);
fprintf('  Heat leak into interior     : %8.2f W\n', Q_leak_in);
fprintf('  Internal dissipation        : %8.1f W\n', Q_int_hot);

fprintf('\n=================================================================\n');
fprintf('  COLD CASE ENERGY BUDGET\n');
fprintf('=================================================================\n');
fprintf('  Earth IR on full surface    : %8.1f W\n', Q_ext_cold);
fprintf('  Skin temperature            : %8.2f K  (%6.2f C)\n', T_skin_cold_K, T_skin_cold_C);
fprintf('  Heat leak out of interior   : %8.2f W\n', Q_leak_out);
fprintf('  Internal dissipation        : %8.1f W\n', Q_int_cold);
fprintf('  Total battery load needed   : %8.1f W\n', Q_heater_budget_total);

fprintf('\n=================================================================\n');
fprintf('  COMPONENT TEMPERATURES  (+-%.0f C design margin applied)\n', margin_C);
fprintf('=================================================================\n');
fprintf('  %-28s %8s %8s %6s   %8s %8s %6s\n', ...
    'Component','T_hot C','Tmax C','HOT?','T_cold C','Tlim C','COLD?');
fprintf('  %s\n', repmat('-', 1, 82));
for i = 1:N
    hs = '  OK ';  cs = '  OK ';
    if hot_viol(i);  hs = '**HOT'; end
    if cold_viol(i); cs = '**CLD'; end
    fprintf('  %-28s %8.1f %8.0f %5s   %8.1f %8.0f %5s\n', ...
        comps(i).name, T_comp_hot_C(i), comps(i).T_op_max_C, hs, ...
        T_comp_cold_C(i), cold_limit(i), cs);
end

fprintf('\n=================================================================\n');
fprintf('  RADIATOR SPECIFICATIONS\n');
fprintf('=================================================================\n');
fprintf('  Coating / emissivity     : white paint or OSR  (eps=%.2f)\n', eps_rad);
fprintf('  Sink temperature         : %.1f C\n', T_sink_K - 273.15);
fprintf('  Available area           : %.2f m^2  (%.0f%% blockage applied)\n', ...
    A_rad_available, frac_blocked*100);
fprintf('\n');
fprintf('  %-28s %10s %10s %10s\n', 'Component','Q_hot (W)','T_max (C)','A_rad (m^2)');
fprintf('  %s\n', repmat('-', 1, 63));
for i = 1:N
    if hot_viol(i)
        fprintf('  %-28s %10.1f %10.0f %10.4f\n', ...
            comps(i).name, comps(i).Q_hot_W, comps(i).T_op_max_C, A_rad_m2(i));
    end
end
fprintf('  %s\n', repmat('-', 1, 63));
fprintf('  %-28s %10s %10s %10.4f\n', 'TOTAL', '', '', A_rad_total);
if A_rad_total == 0
    fprintf('  -> No radiators required — all components within hot limits\n');
elseif A_rad_total <= A_rad_available
    fprintf('  -> Body-surface radiators SUFFICIENT (%.1f%% of available area)\n', ...
        100*A_rad_total/A_rad_available);
else
    fprintf('  -> WARNING: Need +%.2f m^2 additional (deployable panels?)\n', ...
        A_rad_total - A_rad_available);
end

fprintf('\n=================================================================\n');
fprintf('  HEATER SPECIFICATIONS\n');
fprintf('=================================================================\n');
fprintf('  Type            : Kapton film heater\n');
fprintf('  Power density   : %.1f W/cm^2  (TODO: confirm with vendor)\n', heater_density);
fprintf('  Margin          : 20%% per ECSS-E-ST-31C\n\n');
fprintf('  %-28s %10s %10s %8s %10s %10s\n', ...
    'Component','Type','Limit (C)','P (W)','Setpt (C)','Area (cm^2)');
fprintf('  %s\n', repmat('-', 1, 80));
for i = 1:N
    if cold_viol(i)
        htype = 'survival';
        if comps(i).always_on; htype = 'operat.'; end
        fprintf('  %-28s %10s %10.0f %8.1f %10.1f %10.1f\n', ...
            comps(i).name, htype, cold_limit(i), Q_heater_W(i), T_setpt_C(i), A_heat_cm2(i));
    end
end
fprintf('  %s\n', repmat('-', 1, 80));
fprintf('  %-28s %10s %10s %8.1f\n', 'Component heaters', '', '', Q_heater_total);
fprintf('  %-28s %10s %10s %8.1f  (replace MLI loss in eclipse)\n', 'MLI loss replacement', '', '', Q_leak_out);
fprintf('  %-28s %10s %10s %8.1f  (total battery draw in eclipse)\n', 'TOTAL BATTERY LOAD', '', '', Q_heater_budget_total);
fprintf('  TODO: confirm battery supports %.1f W for %.0f min eclipse\n', ...
    Q_heater_budget_total, 0.38*5765/60);
fprintf('=================================================================\n\n');

%% =========================================================================
%  SECTION 12 — PLOTS
% ==========================================================================
plot_cases(comps, T_comp_hot_C, T_comp_cold_C, cold_limit, ...
           hot_viol, cold_viol, ...
           A_rad_m2, A_rad_available, ...
           Q_heater_W, A_heat_cm2, T_setpt_C, ...
           Q_heater_budget_total, margin_C, ...
           T_skin_hot_C, T_skin_cold_C, Q_leak_in, Q_leak_out);
