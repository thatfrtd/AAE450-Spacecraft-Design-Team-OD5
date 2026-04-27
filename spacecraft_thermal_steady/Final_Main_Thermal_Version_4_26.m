%% =========================================================================
%  SPACECRAFT STEADY-STATE THERMAL ANALYSIS
%  AAE450 — Worst-Case Hot & Cold Sizing  (2 x 2 x 2 m)
%  =========================================================================
%
%  HOT CASE — Sunlit orbit, broadside to sun.
%    Solar + Earth IR + albedo on two faces; all components operating.
%    Goal: size the radiator to keep every component ≤ T_op_max.
%
%  COLD CASE — Eclipse, safe mode.
%    Earth IR on one small face only; most components off.
%    Goal: size component heaters to keep every component ≥ T_limit.
%
%  ARCHITECTURE
%    Single shared thermal bus (structural cold plate).
%    All body components connect to bus via their G_W_K conductance.
%    Hot: T_int_hot is set so the tightest component just meets its limit.
%         A single radiator holds the bus at T_int_hot.
%    Cold: T_int_cold is solved from the global energy balance.
%          Component heaters (Kapton film, sized here) maintain each
%          component above its cold limit. No bulk interior heater.
%    Radiator thermal switch isolates the radiator from the bus in eclipse.
%    Panel thermal switch isolates the boom path from the bus in eclipse
%    (prevents cold panels from draining hundreds of watts from the bus).
%
%  SOURCES
%    [S1] Wertz et al., Space Mission Engineering: The New SMAD (2011)
%    [S2] NASA-TP-2002-209914, Spacecraft Thermal Control Handbook
%    [S3] ECSS-E-ST-31C, Space Engineering: Thermal Control
%
%  UNITS: SI throughout (W, m, m^2, K, C)
% =========================================================================

clear; clc; close all;
sigma = 5.6704e-8;   % [W/m^2/K^4]

%% =========================================================================
%  SECTION 1 — ENVIRONMENT
% =========================================================================

% Worst-case hot (perihelion, sub-solar point, max Earth IR)
hot.S       = 1420;  % [W/m^2]  solar flux (perihelion ~1414, rounded up)   [S1]
hot.q_IR    = 258;   % [W/m^2]  max Earth IR (equatorial, sub-solar)        [S1 Tbl 11-48]
hot.albedo  = 0.40;  % [-]      max Earth albedo (clouds/polar ice)         [S1 Tbl 11-48]
hot.F_earth = 0.15;  % [-]      spacecraft-to-Earth view factor             [S2]
                     %          NOTE: F = (R_e/(R_e+h))^2 ≈ 0.15 @ h~2500km
                     %          Adjust if orbit altitude changes.

% Worst-case cold (eclipse, min Earth IR)
cold.q_IR   = 218;   % [W/m^2]  min Earth IR (polar night)                  [S1 Tbl 11-48]
cold.F_earth = 0.15; % [-]      same orbit geometry

% MLI outer surface optical properties (end-of-life, UV-degraded Kapton)
alpha_MLI = 0.50;    % solar absorptivity — EOL degraded   [S2; fresh ~0.10]
eps_MLI   = 0.80;    % IR emissivity — EOL degraded        [S2; fresh ~0.05 for aluminized]
                     % NOTE: These are for the outer MLI/paint surface, not bulk MLI.
                     % A high-emissivity outer coat (e.g. black Kapton) is assumed.
                     % TODO: confirm outer surface finish with thermal coatings team.
T_space = 3.0;       % [K] deep-space background temperature

%% =========================================================================
%  SECTION 2 — GEOMETRY  (2 x 2 x 2 m cube, 6 faces each 4 m^2)
%  TODO: Update if structure dimensions change.
% =========================================================================

A_face  = 4.0;    % [m^2] each face area (cube: L*L = 2*2 = 4)
A_total = 24.0;   % [m^2] total exterior surface (6 faces * 4 m^2)

% Face assignments — fixed for both cases
A_sun_face    = A_face;                              % +Y: broadside to sun (hot)
A_earth_hot   = A_face;                              % +Z: nadir face, IR+albedo (hot)
A_shadow_hot  = A_total - A_sun_face - A_earth_hot;  % remaining 4 faces = 16 m^2 (hot)
A_earth_cold  = A_face;                              % +X: smallest face, IR only (cold)
A_shadow_cold = A_total - A_earth_cold;              % remaining 5 faces = 20 m^2 (cold)

fprintf('GEOMETRY: %.0fx%.0fx%.0f m cube  |  A_total = %.0f m^2\n', ...
    sqrt(A_face), sqrt(A_face), sqrt(A_face), A_total);

%% =========================================================================
%  SECTION 3 — MLI WALL RESISTANCE
%
%  Skin-to-interior resistance = MLI resistance + honeycomb wall resistance.
%  With R_hot ≈ 20 m^2 K/W, the MLI dominates (honeycomb ~0.01 m^2 K/W).
%  Skin temperatures change <1 K for ±30 C shifts in T_int, so using a
%  fixed T_int_ref = 20 C as the fzero BC is valid (<5% error in Q_leak).
%  TODO: confirm C_eff and panel thickness with thermal/structures team.
% =========================================================================

C_eff_hot  = 0.050;  % [W/m^2 K]  MLI effective conductance, hot  [S2; ~0.05 for 20-layer EOL]
C_eff_cold = 0.030;  % [W/m^2 K]  MLI effective conductance, cold [S2; tighter layers, colder]
t_Al = 0.025;        % [m]  honeycomb panel thickness (2.5 cm — typical small sat)
k_Al = 2.5;          % [W/m K]  effective conductivity of Al honeycomb panel

R_hot  = 1/C_eff_hot  + t_Al/k_Al;   % [m^2 K/W]
R_cold = 1/C_eff_cold + t_Al/k_Al;   % [m^2 K/W]

T_int_ref = 293.15;  % [K] = 20 C — used ONLY as fzero BC for skin temps (approximation)

fprintf('MLI resistance: hot = %.3f m^2K/W  |  cold = %.3f m^2K/W\n', R_hot, R_cold);

%% =========================================================================
%  SECTION 4 — PER-FACE SKIN TEMPERATURES
%
%  Each face group has its own radiation balance solved with fzero:
%    Q_absorbed_external - Q_radiated_to_space - Q_conducted_to_interior = 0
%  The area-weighted average skin temp is more accurate than using only
%  the hot face (shadow faces are ~150 K colder and drain heat from interior).
%  Skin temperatures are used to build Q_leak as a function of T_int,
%  which the cold-case fzero needs.
% =========================================================================

% HOT: sun face
f_sun = @(T) alpha_MLI*hot.S*A_sun_face ...
    - eps_MLI*sigma*A_sun_face*(T^4 - T_space^4) ...
    - A_sun_face*(T - T_int_ref)/R_hot;
T_sun_K = fzero(f_sun, 350);

% HOT: Earth face (IR + albedo)
Q_earth_hot_in = (eps_MLI*hot.q_IR + alpha_MLI*hot.S*hot.albedo) * hot.F_earth * A_earth_hot;
f_earth_hot = @(T) Q_earth_hot_in ...
    - eps_MLI*sigma*A_earth_hot*(T^4 - T_space^4) ...
    - A_earth_hot*(T - T_int_ref)/R_hot;
T_earth_hot_K = fzero(f_earth_hot, 250);

% HOT: shadow faces (no external input)
f_shadow_hot = @(T) -eps_MLI*sigma*A_shadow_hot*(T^4 - T_space^4) ...
                    -A_shadow_hot*(T - T_int_ref)/R_hot;
T_shadow_hot_K = fzero(f_shadow_hot, 200);

% Q_leak_in as a function of bus temperature Ti [K]
Q_leak_in_fn = @(Ti) A_sun_face*(T_sun_K - Ti)/R_hot ...
                   + A_earth_hot*(T_earth_hot_K - Ti)/R_hot ...
                   + A_shadow_hot*(T_shadow_hot_K - Ti)/R_hot;
T_skin_hot_C = (A_sun_face*T_sun_K + A_earth_hot*T_earth_hot_K + ...
                A_shadow_hot*T_shadow_hot_K)/A_total - 273.15;

% COLD: Earth face (IR only)
Q_earth_cold_in = eps_MLI*cold.q_IR*cold.F_earth*A_earth_cold;
f_earth_cold = @(T) Q_earth_cold_in ...
    - eps_MLI*sigma*A_earth_cold*(T^4 - T_space^4) ...
    - A_earth_cold*(T - T_int_ref)/R_cold;
T_earth_cold_K = fzero(f_earth_cold, 200);

% COLD: all other faces
f_shadow_cold = @(T) -eps_MLI*sigma*A_shadow_cold*(T^4 - T_space^4) ...
                     -A_shadow_cold*(T - T_int_ref)/R_cold;
T_shadow_cold_K = fzero(f_shadow_cold, 150);

% Q_leak_out as a function of bus temperature Ti [K]
Q_leak_out_fn = @(Ti) A_earth_cold*(Ti - T_earth_cold_K)/R_cold ...
                    + A_shadow_cold*(Ti - T_shadow_cold_K)/R_cold;
T_skin_cold_C = (A_earth_cold*T_earth_cold_K + ...
                 A_shadow_cold*T_shadow_cold_K)/A_total - 273.15;

fprintf('HOT  skin: sun=%.0fC  earth=%.0fC  shadow=%.0fC  avg=%.1fC\n', ...
    T_sun_K-273.15, T_earth_hot_K-273.15, T_shadow_hot_K-273.15, T_skin_hot_C);
fprintf('COLD skin: earth=%.0fC  shadow=%.0fC  avg=%.1fC\n', ...
    T_earth_cold_K-273.15, T_shadow_cold_K-273.15, T_skin_cold_C);

%% =========================================================================
%  SECTION 5 — SOLAR PANEL TEMPERATURE  (decoupled from body MLI)
%
%  Panels solve their own radiation balance. Not inside the MLI envelope.
%  HOT: alpha*S*A = eps*sigma*A*Tp^4 + G_boom*(Tp-T_bus) + Q_elec
%  COLD: eps*q_IR*A*VF + G_boom*(T_bus-Tp) = eps*sigma*A*Tp^4
%
%  G_boom = 5.0 W/K: standard deployment hinge conductance.
%    This causes panels (~110°C hot) to push ~400W into the bus in the
%    hot case (real, must be rejected by radiator).
%    It also drains hundreds of watts in eclipse if left connected.
%    use_panel_thermal_switch = true cuts the cold-case boom path,
%    representing a boom thermal switch or panel rotation out of eclipse.
%
%  TODO: confirm panel area, G_boom with mech team; Q_elec with EPS team.
%  TODO: confirm whether 3624 W is electrical output or waste heat.
% =========================================================================

panel.A      = 2 * 2.904 * 11.657;  % [m²] both panels — updated to team PDF (OD5, April 2026)
                                      % PDF states each panel: 2.904 m × 11.657 m = 33.84 m² each
                                      % Total = 67.7 m²; mass 201.42 kg matches ~3 kg/m² GaAs arrays
                                      % NOTE: code previously used 9.607 m — verify with EPS team
                                      % if 3.624 kW is peak sunlit output vs time-averaged output.
panel.alpha  = 0.80;   % GaAs cell solar absorptivity
panel.eps    = 0.85;   % back-surface IR emissivity (painted)
panel.G_boom = 5.0;    % [W/K] boom conductance — standard hinge; TODO: confirm with mech
panel.VF_IR  = 0.12;   % back-face view factor to Earth
panel.Q_elec = 3624;   % [W] electrical power extracted (team value)

% Panel thermal switch: if true, boom is isolated in eclipse (cold case)
use_panel_thermal_switch = true;   % TODO: confirm switch or equivalent isolation in design

f_pan_hot = @(Tp) panel.alpha*hot.S*panel.A ...
    - panel.eps*sigma*panel.A*Tp^4 ...
    - panel.G_boom*(Tp - T_int_ref) - panel.Q_elec;
T_panel_hot_K = fzero(f_pan_hot, [200 800]);
T_panel_hot_C = T_panel_hot_K - 273.15;

f_pan_cold = @(Tp) panel.eps*cold.q_IR*panel.A*panel.VF_IR ...
    + panel.G_boom*(T_int_ref - Tp) ...
    - panel.eps*sigma*panel.A*Tp^4;
T_panel_cold_K = fzero(f_pan_cold, [10 300]);
T_panel_cold_C = T_panel_cold_K - 273.15;

% Panel boom drain as function of T_int [K] — zero if switch is ON
if use_panel_thermal_switch
    Q_panel_drain_fn = @(Ti) 0;
else
    Q_panel_drain_fn = @(Ti) panel.G_boom * max(0, Ti - T_panel_cold_K);
end

fprintf('Panels: HOT=%.1fC (limit:140C)  |  COLD=%.1fC (survival:-150C)\n', ...
    T_panel_hot_C, T_panel_cold_C);

%% =========================================================================
%  SECTION 6 — COMPONENT DEFINITIONS
% =========================================================================

comps = Final_Components_Version_4_26();
N  = length(comps);
nd = ~[comps.decoupled];   % true = body component, false = decoupled (panels)

Q_int_hot  = sum([comps(nd).Q_hot_W ]);   % body components only
Q_int_cold = sum([comps(nd).Q_cold_W]);

fprintf('Body heat: HOT=%.1fW  |  SAFE-MODE=%.1fW\n', Q_int_hot, Q_int_cold);

%% =========================================================================
%  SECTION 7 — INTERIOR BUS TEMPERATURE (SOLVED)
%
%  Design margin: 5 C per ECSS-E-ST-31C (preliminary phase). [S3]
%  TODO: tighten to 3 C after CDR test data.
%
%  7a HOT: each component constrains the bus:
%     T_bus ≤ T_op_max - Q_hot/G - margin
%     T_int_hot = min over all body components (tightest wins).
%     Radiator is sized to hold bus at exactly this temperature.
%
%  7b COLD: fzero on the energy balance:
%     Q_int_cold + ΣQ_heaters(Ti) = Q_leak_out(Ti) + Q_panel_drain(Ti)
%     Heater power for component i:
%       max(0, G_i*(T_limit_i - Ti) - Q_cold_i)   [raw, before margin]
%     This is self-consistent: fzero includes heaters, so T_int_cold
%     already "knows" what heaters will be needed.
% =========================================================================

margin_C = 5;   % [C] design margin — applied in both directions

% Cold limit per component (operating if always_on, survival if not)
cold_limit_C = zeros(1, N);
cold_limit_K = zeros(1, N);
for i = 1:N
    if comps(i).always_on
        cold_limit_C(i) = comps(i).T_op_min_C;
    else
        cold_limit_C(i) = comps(i).T_surv_min_C;
    end
    cold_limit_K(i) = cold_limit_C(i) + 273.15;
end

% 7a: HOT — derive T_int_hot
T_bus_max_C = zeros(1, N);
for i = 1:N
    if comps(i).decoupled
        T_bus_max_C(i) = Inf;
    else
        T_bus_max_C(i) = comps(i).T_op_max_C - comps(i).Q_hot_W/comps(i).G_W_K - margin_C;
    end
end
[T_int_hot_C, idx_bind] = min(T_bus_max_C);
T_int_hot_K = T_int_hot_C + 273.15;

if T_int_hot_K <= 0
    error('T_int_hot = %.1f K — below absolute zero. Check T_op_max and G_W_K values.', T_int_hot_K);
end

% Hot-case component temperatures
T_comp_hot_C = zeros(1, N);
for i = 1:N
    if comps(i).decoupled
        T_comp_hot_C(i) = T_panel_hot_C + margin_C;
    else
        T_comp_hot_C(i) = T_int_hot_C + comps(i).Q_hot_W/comps(i).G_W_K + margin_C;
    end
end

fprintf('T_int_hot = %.2f C  [binding: %s, T_op_max=%.0fC, Q=%.0fW, G=%.1fW/K]\n', ...
    T_int_hot_C, comps(idx_bind).name, comps(idx_bind).T_op_max_C, ...
    comps(idx_bind).Q_hot_W, comps(idx_bind).G_W_K);

% Panel boom contribution at T_int_hot (re-evaluated at solved bus temp)
Q_panel_to_sc = panel.G_boom * (T_panel_hot_K - T_int_hot_K);

% 7b: COLD — solve T_int_cold
Q_heater_raw = @(i,Ti) nd(i) * max(0, comps(i).G_W_K*(cold_limit_K(i)-Ti) - comps(i).Q_cold_W);
Q_heaters_fn = @(Ti) sum(arrayfun(@(i) Q_heater_raw(i,Ti), 1:N));

f_cold = @(Ti) Q_int_cold + Q_heaters_fn(Ti) - Q_leak_out_fn(Ti) - Q_panel_drain_fn(Ti);
T_int_cold_K = fzero(f_cold, [150 450]);
T_int_cold_C = T_int_cold_K - 273.15;

fprintf('T_int_cold = %.2f C  [solved from eclipse energy balance]\n', T_int_cold_C);

%% =========================================================================
%  SECTION 8 — VIOLATION CHECK
% =========================================================================

hot_viol      = false(1, N);
cold_viol     = false(1, N);
T_comp_cold_C = zeros(1, N);

for i = 1:N
    if T_comp_hot_C(i) > comps(i).T_op_max_C + 0.01
        hot_viol(i) = true;   % should not occur by construction
    end
    if comps(i).decoupled
        T_comp_cold_C(i) = T_panel_cold_C - margin_C;
    else
        T_comp_cold_C(i) = T_int_cold_C + comps(i).Q_cold_W/comps(i).G_W_K - margin_C;
    end
    if T_comp_cold_C(i) < cold_limit_C(i)
        cold_viol(i) = true;
    end
end

%% ========================================================================
%  SECTION 9 — RADIATOR SIZING  (single shared body-mounted radiator)
%
%  The radiator rejects all spacecraft heat in the hot case:
%    Q_rad = Q_int_hot + Q_panel_to_sc + Q_leak_in(T_int_hot)
%
%  Radiator equation: Q_rad = eps_rad * sigma * A_rad * (T_bus^4 - T_sink^4)
%
%  T_sink: effective radiative sink the radiator sees (Earth IR + reflected
%  solar + deep space). Formula derived from eps_rad*sigma*T_sink^4 =
%  eps_rad*q_IR*VF + alpha_rad*S*albedo*VF + eps_rad*sigma*T_space^4.
%
%  Available area = non-sun-facing faces minus hardware blockage.
%  TODO: confirm eps_rad, alpha_rad with coatings vendor.
%  TODO: confirm frac_blocked with mech/propulsion team.
% =========================================================================

eps_rad      = 0.85;   % radiator emissivity (white paint or OSR)
alpha_rad    = 0.12;   % radiator solar absorptivity (white paint ~0.15, OSR ~0.07)
VF_earth_rad = 0.15;   % radiator view factor to Earth (same orbit geometry)
frac_blocked = 0.25;   % fraction of non-sun area blocked by hardware

T_sink_K = (hot.q_IR*VF_earth_rad/sigma ...
          + (alpha_rad/eps_rad)*hot.S*hot.albedo*VF_earth_rad/sigma ...
          + T_space^4)^0.25;

Q_leak_in_at_hot = Q_leak_in_fn(T_int_hot_K);
Q_rad_hot = Q_int_hot + Q_panel_to_sc + Q_leak_in_at_hot;

if Q_rad_hot <= 0
    A_rad_total = 0;
elseif T_int_hot_K <= T_sink_K
    fprintf('!! CRITICAL: T_int_hot (%.1fC) ≤ T_sink (%.1fC) — radiator infeasible.\n', ...
        T_int_hot_C, T_sink_K-273.15);
    fprintf('   Increase G_W_K for: %s\n', comps(idx_bind).name);
    A_rad_total = Inf;
else
    A_rad_total = Q_rad_hot / (eps_rad*sigma*(T_int_hot_K^4 - T_sink_K^4));
end

A_rad_avail = (A_total - A_sun_face) * (1 - frac_blocked);   % = 15 m^2

% Cold-case radiator thermal switch (isolates radiator from bus in eclipse)
use_rad_thermal_switch = true;   % strongly recommended — cuts eclipse radiation loss
Q_rad_cold_loss = 0;             % = 0 when switch ON; set below if switch OFF
if ~use_rad_thermal_switch && ~isinf(A_rad_total) && A_rad_total > 0
    Q_rad_cold_loss = eps_rad * sigma * A_rad_total * (T_int_cold_K^4 - T_space^4);
    f_cold2 = @(Ti) Q_int_cold + Q_heaters_fn(Ti) ...
        - Q_leak_out_fn(Ti) - Q_panel_drain_fn(Ti) - Q_rad_cold_loss;
    T_int_cold_K = fzero(f_cold2, [150 450]);
    T_int_cold_C = T_int_cold_K - 273.15;
    for i = 1:N
        if ~comps(i).decoupled
            T_comp_cold_C(i) = T_int_cold_C + comps(i).Q_cold_W/comps(i).G_W_K - margin_C;
            cold_viol(i)     = T_comp_cold_C(i) < cold_limit_C(i);
        end
    end
end

%% =========================================================================
%  SECTION 10 — COMPONENT HEATER SIZING
%
%  Each cold-violating component gets a Kapton film heater sized to close
%  the gap between what the cold bus provides and what the component needs:
%    Q_heater_i = max(0, G_i*(T_limit_i - T_int_cold) - Q_cold_i) * 1.20
%  The 20% margin is per ECSS-E-ST-31C. [S3]
%  Thermostat setpoint is 2 C above the cold limit.
%  Heater area = power / density (1 W/cm^2 typical Kapton — TODO: confirm).
% =========================================================================

heater_density = 1.0;   % [W/cm^2]  TODO: confirm with Kapton heater vendor
Q_heater_W = zeros(1, N);
T_setpt_C  = zeros(1, N);
A_heat_cm2 = zeros(1, N);

for i = 1:N
    if comps(i).decoupled; continue; end
    if cold_viol(i)
        Q_heater_W(i) = Q_heater_raw(i, T_int_cold_K) * 1.20;
        T_setpt_C(i)  = cold_limit_C(i) + 2;
        A_heat_cm2(i) = Q_heater_W(i) / heater_density;
    end
end

Q_heater_total  = sum(Q_heater_W);
Q_leak_out_cold = Q_leak_out_fn(T_int_cold_K);
Q_panel_drain_c = Q_panel_drain_fn(T_int_cold_K);
Q_batt_load     = Q_heater_total + Q_leak_out_cold + Q_panel_drain_c + Q_rad_cold_loss;

%% =========================================================================
%  SECTION 11 — CONSOLE REPORT
% =========================================================================

fprintf('\n=== HOT CASE ENERGY BUDGET ==========================================\n');
fprintf('  Solar on sun face             : %8.1f W\n', alpha_MLI*hot.S*A_sun_face);
fprintf('  Earth IR + albedo (earth face): %8.1f W\n', Q_earth_hot_in);
fprintf('  Net MLI leak in (at T_int_hot): %8.1f W\n', Q_leak_in_at_hot);
fprintf('  Internal component heat       : %8.1f W\n', Q_int_hot);
fprintf('  Panel boom to bus             : %8.1f W\n', Q_panel_to_sc);
fprintf('  ----------------------------------\n');
fprintf('  Total rejected by radiator    : %8.1f W\n', Q_rad_hot);
fprintf('  T_int_hot (bus target)        : %8.2f C\n', T_int_hot_C);
fprintf('  T_sink (radiator environment) : %8.2f C\n', T_sink_K-273.15);
fprintf('  Required radiator area        : %8.4f m^2\n', A_rad_total);
fprintf('  Available radiator area       : %8.2f m^2  (%.0f%% blockage)\n', ...
    A_rad_avail, frac_blocked*100);
if isinf(A_rad_total)
    fprintf('  STATUS: ** INFEASIBLE — increase G_W_K for %s\n', comps(idx_bind).name);
elseif A_rad_total == 0
    fprintf('  STATUS: No radiator required\n');
elseif A_rad_total <= A_rad_avail
    fprintf('  STATUS: SUFFICIENT (%.1f%% of available area)\n', 100*A_rad_total/A_rad_avail);
else
    fprintf('  STATUS: ** INSUFFICIENT — need +%.3f m^2\n', A_rad_total - A_rad_avail);
end

fprintf('\n=== COLD CASE ENERGY BUDGET =========================================\n');
fprintf('  Safe-mode internal heat       : %8.1f W\n', Q_int_cold);
fprintf('  Component heater power        : %8.1f W\n', Q_heater_total);
fprintf('  MLI heat leak out             : %8.1f W\n', Q_leak_out_cold);
fprintf('  Panel boom drain              : %8.1f W  (%s)\n', Q_panel_drain_c, ...
    ternary_str(use_panel_thermal_switch, 'switch ON', 'switch OFF'));
fprintf('  Radiator cold loss            : %8.1f W  (%s)\n', Q_rad_cold_loss, ...
    ternary_str(use_rad_thermal_switch,  'switch ON', 'switch OFF'));
fprintf('  ----------------------------------\n');
fprintf('  Total eclipse battery load    : %8.1f W\n', Q_batt_load);
fprintf('  T_int_cold (solved)           : %8.2f C\n', T_int_cold_C);

fprintf('\n=== COMPONENT TEMPERATURES  (+/-%dC margin) ==========================\n', margin_C);
fprintf('  %-26s %7s %6s %5s   %7s %6s %5s\n', ...
    'Component','T_hot','Tmax','HOT?','T_cold','Tlim','CLD?');
fprintf('  %s\n', repmat('-',1,75));
for i = 1:N
    hs='  ok'; cs='  ok';
    if hot_viol(i);  hs='**HOT'; end
    if cold_viol(i); cs='**CLD'; end
    fprintf('  %-26s %7.1f %6.0f %5s   %7.1f %6.0f %5s\n', ...
        comps(i).name, T_comp_hot_C(i), comps(i).T_op_max_C, hs, ...
        T_comp_cold_C(i), cold_limit_C(i), cs);
end

fprintf('\n=== HEATERS ==========================================================\n');
fprintf('  %-26s %8s %8s %10s\n', 'Component','P [W]','Setpt[C]','Area[cm^2]');
fprintf('  %s\n', repmat('-',1,60));
n_htr = 0;
for i = 1:N
    if cold_viol(i) && ~comps(i).decoupled
        n_htr = n_htr + 1;
        fprintf('  %-26s %8.1f %8.1f %10.1f\n', ...
            comps(i).name, Q_heater_W(i), T_setpt_C(i), A_heat_cm2(i));
    end
end
if n_htr == 0
    fprintf('  (None required — T_int_cold warm enough)\n');
end
fprintf('  %s\n', repmat('-',1,60));
fprintf('  TOTAL HEATER POWER         : %.1f W\n', Q_heater_total);
fprintf('  TOTAL ECLIPSE BATTERY LOAD : %.1f W\n', Q_batt_load);
fprintf('=====================================================================\n\n');

%% =========================================================================
%  SECTION 12 — PLOTS
% =========================================================================
Final_Plots_Version_4_26( ...
    comps, T_comp_hot_C, T_comp_cold_C, cold_limit_C, ...
    hot_viol, cold_viol, nd, T_bus_max_C, ...
    A_rad_total, A_rad_avail, ...
    Q_heater_W, A_heat_cm2, T_setpt_C, Q_heater_total, ...
    margin_C, T_skin_hot_C, T_skin_cold_C, ...
    Q_leak_in_fn(T_int_hot_K), Q_leak_out_cold, ...
    Q_int_cold, Q_batt_load, Q_rad_cold_loss, ...
    Q_panel_drain_c, Q_rad_hot, T_int_hot_C, T_int_cold_C, idx_bind);

%% =========================================================================
%  SECTION 13 — THERMAL CONTROL SYSTEM (TCS) MASS ESTIMATE
%
%  Parametric mass budget for ALL thermal hardware on the spacecraft.
%  Unit masses sourced from:
%    [M1] Gilmore, "Spacecraft Thermal Control Handbook Vol. I," 2nd ed. (2002)
%         — Tables 5-2 (MLI), 6-1 (radiators), 7-3 (cold plates), 8-1 (heat pipes)
%    [M2] Wertz et al., Space Mission Engineering: The New SMAD (2011) Ch. 11
%    [M3] Bale et al., "Spacecraft Thermal Control Design," ESA SP-1196 (1996)
%    [M4] Technology Applications Inc. (TAI) — thermal strap product data (2023)
%    [M5] Minco Products Inc. — Kapton heater engineering data (2023)
%    [M6] RUAG Space — paraffin thermal switch product sheet (2020)
%    [M7] NASA GSFC Thermal Engineering Branch design guidelines (2019)
% =========================================================================

fprintf('\n=== SECTION 13: THERMAL CONTROL SYSTEM MASS ESTIMATE ================\n');

% ---- 1. MLI Blankets -------------------------------------------------------
%  20-layer aluminized Kapton/Mylar MLI, including fabric, spacers, attachments.
%  Unit mass: 0.65 kg/m² for assembled flight MLI blanket [M1, Table 5-2].
%  Applied to ~85% of body surface (remaining 15% = access panels, vents, ports).
%  The radiator area is NOT covered with MLI; subtract it from body area.
%  Add 2 m² for boom MLI wraps and harness blankets [M7].
% ---------------------------------------------------------------------------
m_MLI_per_m2 = 0.65;    % [kg/m²]  assembled 20-layer MLI [M1 Table 5-2]
A_MLI_rad    = isinf(A_rad_total) * 0 + ~isinf(A_rad_total) * min(max(A_rad_total,0), A_total);
A_MLI_body   = (A_total - A_MLI_rad) * 0.85;
A_MLI_extra  = 2.0;      % [m²]  boom wraps, plumbing blankets, thruster MLI
A_MLI_total  = A_MLI_body + A_MLI_extra;
m_MLI        = A_MLI_total * m_MLI_per_m2;

% ---- 2. Radiator -----------------------------------------------------------
%  Body-mounted passive radiator: Al-honeycomb panel + white paint (AZ-93) or OSR tiles.
%  Unit mass: 4.0 kg/m² for a dedicated flat radiator panel [M1, Table 6-1].
%  Note: if painted onto existing structure face, use 1.5 kg/m² instead.
%  Conservative 4.0 kg/m² used here (assumes dedicated panel).
% ---------------------------------------------------------------------------
m_rad_per_m2 = 4.0;     % [kg/m²]  [M1 Table 6-1]
if isinf(A_rad_total) || A_rad_total == 0
    m_radiator = 0;
else
    m_radiator = min(A_rad_total, A_rad_avail) * m_rad_per_m2;
end

% ---- 3. Kapton Film Heaters ------------------------------------------------
%  One heater circuit per cold-violating component: polyimide (Kapton) film
%  element + thermostat + local wiring loom.
%  Unit mass: 0.25 kg/circuit [M5; Minco HK-series spacecraft heater data].
%  Add a 20% system harness overhead for sensor leads and controller wiring.
% ---------------------------------------------------------------------------
m_per_heater   = 0.25;   % [kg/circuit]  [M5 Minco Kapton heater data]
n_heater_circs = sum(cold_viol & ~[comps.decoupled]);
m_heaters      = n_heater_circs * m_per_heater * 1.20;   % 20% wiring overhead

% ---- 4. Thermal Switches ---------------------------------------------------
%  Paraffin-actuated bi-metallic switches (RUAG Space or equivalent).
%  Radiator thermal switch + solar panel boom switch = 2 total.
%  Unit mass: 0.45 kg/switch [M6; RUAG Space thermal switch datasheet].
% ---------------------------------------------------------------------------
m_per_switch = 0.45;    % [kg/switch]  [M6 RUAG Space datasheet]
n_switches   = double(use_rad_thermal_switch) + double(use_panel_thermal_switch);
m_switches   = n_switches * m_per_switch;

% ---- 5. Cold Plates (PPU + Battery) ----------------------------------------
%  Al-6061 machined cold plate with embedded copper heat-spreader and
%  Bergquist GP3000 TIM (3 W/m·K) at each contact interface.
%  Unit mass: 4.5 kg/m² including TIM and mounting hardware [M1, Table 7-3].
%  PPU footprint: 410×510 mm (team PDF).  Battery footprint: 434×400 mm × 2 (team PDF).
% ---------------------------------------------------------------------------
m_coldplate_per_m2 = 4.5;                    % [kg/m²]  [M1 Table 7-3]
A_PPU_cp           = 0.410 * 0.510;          % [m²]  from team PDF
A_batt_cp          = 2 * (0.434 * 0.400);   % [m²]  2× battery (team PDF)
m_cold_plates      = (A_PPU_cp + A_batt_cp) * m_coldplate_per_m2;

% ---- 6. Thermal Straps -----------------------------------------------------
%  Braided copper heat straps (TAI CHTS series or Thermacore graphite straps)
%  for high-conductance flexible joints (PPU→cold plate, battery→cold plate,
%  star tracker mount, reaction wheels).
%  Unit mass: 0.20 kg/strap, average 30 cm length [M4; TAI product data].
% ---------------------------------------------------------------------------
m_per_strap = 0.20;   % [kg/strap]  [M4 TAI CHTS data]
n_straps    = 5;      % PPU(2) + battery(2) + star tracker(1) — estimate
m_straps    = n_straps * m_per_strap;

% ---- 7. Heat Pipes ---------------------------------------------------------
%  Axial groove or sintered-wick aluminum/ammonia heat pipes for
%  thermal bus equalization (typical: 2 pipes, ~1.5 m each).
%  Unit mass: 0.18 kg/m for 12.7 mm Al/NH3 heat pipe [M1, Table 8-1].
% ---------------------------------------------------------------------------
m_per_m_hp  = 0.18;   % [kg/m]  12.7 mm Al/NH3 heat pipe [M1 Table 8-1]
L_hp_total  = 3.0;    % [m]  estimate: 2 pipes × 1.5 m each along cold-plate bus
m_heat_pipes = L_hp_total * m_per_m_hp;

% ---- 8. Temperature Sensors and Harness ------------------------------------
%  ~2 thermistors/RTDs per component for redundancy; 5 g each [M7].
%  Thermal sensor harness (point-to-point loom to OBC): 2.0 kg estimate [M2].
% ---------------------------------------------------------------------------
n_sensors  = N * 2;           % 2 per component (primary + redundant)
m_sensors  = n_sensors * 0.005;  % [kg]  5 g per thermistor [M7]
m_wiring   = 2.0;             % [kg]  thermal harness to OBC [M2]

% ---- 9. Thermal Interface Material (TIM) -----------------------------------
%  Bergquist GP3000 (3 W/m·K) pads at all bolted joints not covered above.
%  Negligible per joint but add a bulk allowance.
% ---------------------------------------------------------------------------
m_TIM = 0.20;   % [kg]  all remaining TIM pads [M1]

% ---- TOTAL -----------------------------------------------------------------
m_TCS_components = [m_MLI, m_radiator, m_heaters, m_switches, m_cold_plates, ...
                    m_straps, m_heat_pipes, m_sensors, m_wiring, m_TIM];
m_TCS = sum(m_TCS_components);

fprintf('  %-40s : %6.2f kg  (%.1f m², %.2f kg/m²)\n', ...
    'MLI blankets', m_MLI, A_MLI_total, m_MLI_per_m2);
if m_radiator > 0
    fprintf('  %-40s : %6.2f kg  (%.2f m², %.1f kg/m²)\n', ...
        'Radiator panel', m_radiator, min(A_rad_total,A_rad_avail), m_rad_per_m2);
else
    fprintf('  %-40s : %6.2f kg  (not required)\n', 'Radiator panel', m_radiator);
end
fprintf('  %-40s : %6.2f kg  (%d circuits)\n', 'Kapton heaters + wiring', m_heaters, n_heater_circs);
fprintf('  %-40s : %6.2f kg  (%d switches)\n', 'Thermal switches (rad + boom)', m_switches, n_switches);
fprintf('  %-40s : %6.2f kg  (PPU + 2× battery)\n', 'Cold plates + TIM', m_cold_plates);
fprintf('  %-40s : %6.2f kg  (%d straps)\n', 'Thermal straps (copper braid)', m_straps, n_straps);
fprintf('  %-40s : %6.2f kg  (%.1f m of Al/NH3 pipe)\n', 'Heat pipes (bus equalizer)', m_heat_pipes, L_hp_total);
fprintf('  %-40s : %6.2f kg  (%d sensors)\n', 'Temperature sensors (thermistors)', m_sensors, n_sensors);
fprintf('  %-40s : %6.2f kg\n', 'Thermal sensor harness', m_wiring);
fprintf('  %-40s : %6.2f kg\n', 'TIM pads (remaining joints)', m_TIM);
fprintf('  %s\n', repmat('-',1,55));
fprintf('  %-40s : %6.2f kg\n', 'TOTAL TCS MASS (estimate)', m_TCS);
fprintf('\n');
fprintf('  Uncertainty: ±30%% at this design phase (preliminary).\n');
fprintf('  Lower bound (±30%%): %.2f kg | Upper bound: %.2f kg\n', m_TCS*0.70, m_TCS*1.30);
fprintf('\n');
fprintf('  Sources:\n');
fprintf('    MLI / radiator / cold plate unit masses: Gilmore (2002) Tables 5-2, 6-1, 7-3\n');
fprintf('    Heat pipe: Gilmore (2002) Table 8-1  (12.7 mm Al/NH3, 0.18 kg/m)\n');
fprintf('    Heaters: Minco Products HK-series spacecraft heater data (2023)\n');
fprintf('    Thermal switches: RUAG Space paraffin switch datasheet (2020)\n');
fprintf('    Thermal straps: Technology Applications Inc. CHTS data (2023)\n');
fprintf('    Sensors/harness: SMAD (2011) Ch.11; NASA GSFC guidelines (2019)\n');
fprintf('======================================================================\n\n');

%% =========================================================================
%  LOCAL HELPER
% =========================================================================
function s = ternary_str(cond, a, b)
    if cond; s = a; else; s = b; end
end
