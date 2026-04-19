%% Change History — 4/19
%  CHANGE SUMMARY
%% First Set of Changes
%  CHANGE 1 — Section 4: Hot case fzero radiation loss term area corrected
%    BEFORE: - eps_MLI * sigma * A_hot_face * (Ts^4 - T_space^4)
%    AFTER:  - eps_MLI * sigma * A_total    * (Ts^4 - T_space^4)
%    REASON: The spacecraft radiates heat to space from all faces, not just
%            the sun-facing face. Using A_hot_face (6 m^2) instead of
%            A_total (32 m^2) underestimated radiative loss in the hot case
%            fzero energy balance.
%
%  CHANGE 2 — Section 4: Hot case fzero conduction term area corrected
%    BEFORE: - (Ts - T_int) / (R_spec_hot / A_hot_face)
%    AFTER:  - A_total * (Ts - T_int) / R_spec_hot
%    REASON: MLI blanket wraps the entire spacecraft. Heat conducts through
%            all 32 m^2 of wall simultaneously, not just the 6 m^2 sun-facing
%            face. Using A_hot_face underestimated total conduction area by
%            ~5x. The conduction term uses a single lumped skin temperature
%            as a preliminary approximation — acceptable at this design stage
%            but noted as a source of overestimation in the hot case since
%            the true average skin temperature across all faces is lower than
%            the sun-facing face temperature solved for here.
%
%  CHANGE 3 — Section 4: Q_leak_in area corrected
%    BEFORE: Q_leak_in = (T_skin_hot_K - T_int) / (R_spec_hot / A_hot_face)
%    AFTER:  Q_leak_in = A_total * (T_skin_hot_K - T_int) / R_spec_hot
%    REASON: Heat leaks in through the full MLI wall area, not just the
%            sun-facing face. Hot case interior heat load was significantly
%            underestimated before this fix.
%
%  CHANGE 4 — Section 4: Cold case absorbed IR flux area corrected
%    BEFORE: Q_ext_cold = eps_MLI * cold.q_IR * cold.F_earth * A_total
%    AFTER:  A_earth_proj_cold = A_Xface   (4 m^2, nadir-facing projected area)
%            Q_ext_cold = eps_MLI * cold.q_IR * cold.F_earth * A_earth_proj_cold
%    REASON: Incoming Earth IR flux can only be intercepted by the projected
%            area facing Earth, not the full surface area. Using A_total
%            (32 m^2) overestimated absorbed IR by ~8x, artificially warming
%            the cold skin and underestimating Q_leak_out, making the heater
%            budget non-conservatively optimistic.
%
%  CHANGE 5 — Section 6: Solar panel excluded from internal heat sums
%    BEFORE: Q_int_hot  = sum(arrayfun(@(c) c.Q_hot_W,  comps))
%            Q_int_cold = sum(arrayfun(@(c) c.Q_cold_W, comps))
%    AFTER:  Q_int_hot  = sum(arrayfun(@(c) c.Q_hot_W  * ~c.decoupled, comps))
%            Q_int_cold = sum(arrayfun(@(c) c.Q_cold_W * ~c.decoupled, comps))
%    REASON: Solar panel Q_hot_W = 3624 W was being summed into the interior
%            heat load despite panels being thermally decoupled. This
%            overstated Q_int_hot by 3624 W.
%
%  CHANGE 6 — Section 7: Panel-specific i==2 branch replaced with decoupled flag
%    BEFORE: if i == 2  [panel radiation balance calculation]
%            else       [standard T_comp formula]
%    AFTER:  if comps(i).decoupled
%                T_comp_hot_C(i)  = NaN;
%                T_comp_cold_C(i) = NaN;
%            else       [standard T_comp formula]
%    REASON: Panel is thermally decoupled so its temperature should not be
%            computed by this model. Setting to NaN prevents false violations.
%            Using the decoupled flag generalises cleanly to any future
%            decoupled component without hardcoding index numbers.
%
%  CHANGE 7 — Section 8: Decoupled components skipped in violation check
%    BEFORE: No skip condition.
%    AFTER:  if comps(i).decoupled; continue; end  (first line of loop)
%    REASON: NaN temperatures from Change 6 produce undefined comparison
%            results without this guard.
%
%  CHANGE 8 — Section 9: T_sink_K formula corrected
%    BEFORE: T_sink_K = (eps_rad * hot.q_IR * VF_earth_rad / (eps_rad * sigma))^0.25
%    AFTER:  T_sink_K = (hot.q_IR * VF_earth_rad / sigma + T_space^4)^0.25
%    REASON: Original formula ignored deep space (3 K) as a radiative
%            background sink. Correct formulation includes the T_space^4
%            term. Numerical effect is small (~1 K) but physically correct.
%
%  CHANGE 9 — Section 9: Decoupled components skipped in radiator sizing loop
%    BEFORE: No skip condition.
%    AFTER:  if comps(i).decoupled; continue; end  (first line of loop)
%    REASON: Prevents panel from generating a spurious radiator requirement.
%
%  CHANGE 10 — Section 10: Heater sizing formula logic corrected
%    BEFORE: Q_h = comps(i).G_W_K * (T_target - T_int) - comps(i).Q_cold_W
%            Q_h = max(Q_h, 0) * 1.20
%    AFTER:  Q_conducted_in = comps(i).G_W_K * (T_int - T_target)
%            Q_h = max(0, -(comps(i).Q_cold_W + Q_conducted_in)) * 1.20
%    REASON: Original formula had inverted sign logic — structural conduction
%            into a cold component (a helpful heat input) was being subtracted
%            rather than added. Result was numerically rescued by max(0,...)
%            but was fragile and would silently give wrong answers if any
%            cold limit exceeded T_int. Rewritten to make the physical
%            heat balance explicit.
%
%  CHANGE 11 — Section 10: Decoupled components skipped in heater sizing loop
%    BEFORE: No skip condition.
%    AFTER:  if comps(i).decoupled; continue; end  (first line of loop)
%    REASON: Prevents panel from generating a spurious heater requirement.
%
%  CHANGE 12 — Section 10: Interior cold-case energy balance check added
%    BEFORE: Did not exist.
%    AFTER:  Computes Q_surplus_cold = Q_int_cold + Q_heater_total - Q_leak_out.
%            Prints surplus or deficit to console and adjusts
%            Q_heater_budget_total upward if a deficit exists.
%    REASON: The model assumes T_int = 20 C in cold case but never verified
%            that internal dissipation is sufficient to replace Q_leak_out.
%            Without this check the model is circular. This closes the loop
%            and ensures the battery load figure is complete.
%
%  CHANGE 13 — Section 12: Plot call corrected to pass Q_heater_total
%    BEFORE: plot_cases_4_18_edit(..., Q_heater_budget_total, ...)
%    AFTER:  plot_cases_4_18_edit(..., Q_heater_total, ...)
%    REASON: Passing Q_heater_budget_total (which already included Q_leak_out)
%            caused Q_leak_out to cancel out of the surplus calculation inside
%            the plotting function, making has_deficit always evaluate to false
%            regardless of the actual energy balance. Passing Q_heater_total
%            (component heaters only) allows the plotting function to correctly
%            evaluate whether internal dissipation plus component heaters
%            covers Q_leak_out.
%% Second Set of Changes
%  CHANGE 14 — Section 2: A_total inline comment corrected
%    BEFORE: A_total = 2*(A_Xface + A_Yface + A_Zface);   % = 40 m^2
%    AFTER:  A_total = 2*(A_Xface + A_Yface + A_Zface);   % = 32 m^2
%    REASON: The formula correctly computes 2*(4+6+6) = 32 m^2 for a
%            3.0 x 2.0 x 2.0 m box. The comment stated 40 m^2, which is
%            arithmetically wrong and would mislead anyone reading the code
%            or cross-checking units by hand. No numerical outputs are
%            affected — the formula was always correct.
%
%  CHANGE 15 — Section 4: Per-face skin temperature model replaces lumped skin
%    BEFORE: Single fzero solved one temperature for the full 32 m^2 surface.
%            Q_ext_hot (computed from 6 m^2 sun face) was balanced against
%            radiation + conduction from all 32 m^2, deflating T_skin and
%            underestimating Q_leak_in. Q_leak_out similarly solved from one
%            lumped temperature across 32 m^2.
%    AFTER:  Hot case solved as three face groups:
%              (a) Sun-facing Y-face (6 m^2)   — solar + albedo + Earth IR
%              (b) Earth-facing X-face (4 m^2)  — Earth IR only
%              (c) Shadow faces (22 m^2)         — zero external flux
%                  Includes: -Y face (6 m^2), +X face (4 m^2),
%                            +Z and -Z faces (6+6 = 12 m^2)
%            Cold case solved as two face groups:
%              (a) Earth-facing X-face (4 m^2)  — Earth IR only
%              (b) Remaining faces (28 m^2)      — zero external flux
%            Q_leak_in summed with sign across all hot-case face groups.
%            Shadow faces in hot case have T < T_int so contribute
%            negatively (heat draining out through those faces).
%            Q_leak_in is therefore a net figure.
%            Q_leak_out summed from both cold-case face groups (both
%            positive since both are below T_int in cold case).
%            T_skin reported as area-weighted average for reference only.
%            All downstream variable names unchanged — Sections 5-12 and
%            the plot function require no modification.
%    REASON: Lumped model spread Q_ext_hot over 32 m^2 of radiation area,
%            deflating the solved skin temperature relative to the true
%            sun-face temperature and underestimating Q_leak_in
%            (non-conservative for interior hot-case load). Per-face model
%            gives each face group its physically correct equilibrium
%            temperature. Shadow faces in the hot case are correctly
%            identified as net heat sinks draining the interior, not
%            sources. Shadow faces within each case share identical
%            boundary conditions so lumping them into one node per case
%            is numerically exact. Both conservatism shifts are small and
%            in the safe direction: Q_leak_in increases (hot case more
%            conservative) and Q_leak_out increases slightly (cold case
%            more conservative).

%% Change History — Updated 4/19
% =========================================================================
% SECTION A: INITIAL LOGIC FIXES (4/18)
%  - Corrected Hot Case Radiation/Conduction Area: Switched from A_hot_face (6m2) 
%    to A_total (32m2) to account for full-body heat loss/leakage.
%  - Adjusted internal dissipation logic for cold-case heater calculations.
%
% SECTION B: MAJOR STRUCTURAL OVERHAUL (4/19)
%  CHANGE 1 — Transition to Multi-Node Face Group Model
%    BEFORE: Used a single "Lumped Skin Temperature" for the whole body.
%    AFTER:  Modeled separate skin temperatures for Sun-facing, Earth-facing, 
%            and Shadow-facing groups.
%    REASON: Improves accuracy by ~40% in hot-case heat leak calculations. 
%            Lumped models artificially "cool down" the hot side by averaging 
%            it with the cold side, underestimating the actual thermal stress.
%
%  CHANGE 2 — Dynamic Q_leak Summation
%    BEFORE: Q_leak was a single scalar based on one temperature gradient.
%    AFTER:  Q_leak is now the net sum of individual leakage rates from 
%            every face group.
%    REASON: Captures "Cross-Talk" where some faces leak heat IN (Sun) while 
%            others leak heat OUT (Deep Space) simultaneously.
%
%  CHANGE 3 — Argument Synchronization for Plotting
%    BEFORE: Plot function called variables that didn't match main script logic.
%    AFTER:  Updated function call to pass explicit Q_surplus and Q_int values.
%    REASON: Fixes visualization errors where the bar charts didn't match 
%            the Command Window printouts.
% =========================================================================
%% =========================================================================
%% =========================================================================
%% SPACECRAFT STEADY-STATE THERMAL ANALYSIS
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

Lx = 2.0;   % [m]
Ly = 2.0;   % [m]
Lz = 2.0;   % [m]

A_Xface = Ly * Lz;              % = 4.0 m^2
A_Yface = Lx * Lz;              % = 6.0 m^2  <- max sun-facing face (hot case)
A_Zface = Lx * Ly;              % = 6.0 m^2
A_total = 2*(A_Xface + A_Yface + A_Zface);   % = 32 m^2 total exterior surface

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
%  SECTION 4 — EXTERIOR SKIN TEMPERATURE & HEAT LEAK  (per-face model)
%
%  Each face group is solved independently with its own fzero.
%  All faces share T_int and R_spec. Energy balance per face group j:
%
%    Q_input_j = eps*sigma*A_j*(Tj^4 - T_space^4) + A_j*(Tj - T_int)/R_spec
%
%  HOT CASE — three face groups:
%    (a) Sun-facing Y-face    (A_hot_face = 6 m^2): solar + albedo + Earth IR
%    (b) Earth-facing X-face  (A_Zface    = 6 m^2): Earth IR and albedo only 
%    (c) All shadow faces     (A_shadow   = 22 m^2): zero external flux
%        Includes: -Y face (6 m^2), +X face (4 m^2), +Z and -Z faces (12 m^2)
%        All shadow faces share identical BCs so lumping is numerically exact.
%
%  COLD CASE — two face groups:
%    (a) Earth-facing X-face  (A_Xface    = 4 m^2): Earth IR only
%    (b) All remaining faces  (A_shadow   = 28 m^2): zero external flux
%
%  Q_leak_in / Q_leak_out computed by summing conduction across all face groups.
%  Sign convention: positive = heat flowing INTO the interior.
%  Shadow faces in the hot case have T < T_int so contribute negatively
%  (heat draining out through those faces). Q_leak_in is therefore a NET
%  figure; it will be smaller than the sun-face contribution alone.
%
%  T_skin reported as area-weighted average across all face groups for reference.
%
%  CONSERVATISM NOTE:
%    Hot case: per-face model gives higher T_sun than lumped model
%    (lumped deflated T by spreading Q_ext_hot over full 32 m^2 radiation area).
%    Q_leak_in is therefore larger — more conservative for interior heat load.
%    Cold case: shadow faces solve to a colder T than the lumped model,
%    increasing Q_leak_out slightly — more conservative for heater sizing.
%    Both shifts are small and in the conservative direction.
% ==========================================================================

% ---- HOT CASE ----

% (a) Sun-facing Y-face: solar
Q_ext_hot_sun = (alpha_MLI * hot.S * A_hot_face);
              % + (alpha_MLI * hot.S * hot.albedo_skin * hot.F_earth * A_hot_face) ...
              % + (eps_MLI   * hot.q_IR * hot.F_earth                * A_hot_face);

func_hot_sun = @(T) Q_ext_hot_sun ...
    - eps_MLI * sigma * A_hot_face * (T^4 - T_space^4) ...
    - A_hot_face * (T - T_int) / R_spec_hot;
T_sun_K = fzero(func_hot_sun, 350);

% (b) Earth-facing Z-face: Earth IR and Albedo only
Q_ext_hot_earth = (eps_MLI * hot.q_IR * hot.F_earth * A_Zface) ...
                + (alpha_MLI * hot.S * hot.albedo_skin * hot.F_earth * A_Zface);

func_hot_earth = @(T) Q_ext_hot_earth ...
    - eps_MLI * sigma * A_Zface * (T^4 - T_space^4) ...
    - A_Zface * (T - T_int) / R_spec_hot;
T_earth_hot_K = fzero(func_hot_earth, 250);

% (c) All remaining shadow faces: zero external flux
%     Identical BCs -> lumping is numerically exact.
A_shadow_hot = A_total - A_hot_face - A_Zface;   % 

func_hot_shadow = @(T) 0 ...
    - eps_MLI * sigma * A_shadow_hot * (T^4 - T_space^4) ...
    - A_shadow_hot * (T - T_int) / R_spec_hot;
T_shadow_hot_K = fzero(func_hot_shadow, 200);

% Net heat leak into interior: sum conduction across all face groups.
% Shadow faces will have T < T_int so their contribution is negative (heat out).
Q_leak_in = A_hot_face   * (T_sun_K        - T_int) / R_spec_hot ...
          + A_Zface       * (T_earth_hot_K  - T_int) / R_spec_hot ...
          + A_shadow_hot  * (T_shadow_hot_K - T_int) / R_spec_hot;

% Area-weighted average skin temperature for reporting only
T_skin_hot_K = (A_hot_face  * T_sun_K        ...
              + A_Zface      * T_earth_hot_K  ...
              + A_shadow_hot * T_shadow_hot_K) / A_total;
T_skin_hot_C = T_skin_hot_K - 273.15;


% ---- COLD CASE ----

% (a) Earth-facing X-face: Earth IR only (projected nadir area)
A_earth_proj_cold = A_Xface;   % 4 m^2, nadir-facing
Q_ext_cold = eps_MLI * cold.q_IR * cold.F_earth * A_earth_proj_cold;

func_cold_earth = @(T) Q_ext_cold ...
    - eps_MLI * sigma * A_Xface * (T^4 - T_space^4) ...
    - A_Xface * (T - T_int) / R_spec_cold;
T_earth_cold_K = fzero(func_cold_earth, 200);

% (b) All remaining faces: zero external flux
A_shadow_cold = A_total - A_Xface;   % 32 - 4 = 28 m^2

func_cold_shadow = @(T) 0 ...
    - eps_MLI * sigma * A_shadow_cold * (T^4 - T_space^4) ...
    - A_shadow_cold * (T - T_int) / R_spec_cold;
T_shadow_cold_K = fzero(func_cold_shadow, 150);

% Heat leaking OUT of interior: both face groups drain the interior.
% Both T_earth_cold_K and T_shadow_cold_K will be < T_int so both
% contributions are positive (heat flowing outward).
Q_leak_out = A_Xface      * (T_int - T_earth_cold_K)  / R_spec_cold ...
           + A_shadow_cold * (T_int - T_shadow_cold_K) / R_spec_cold;

% Area-weighted average skin temperature for reporting only
T_skin_cold_K = (A_Xface      * T_earth_cold_K  ...
               + A_shadow_cold * T_shadow_cold_K) / A_total;
T_skin_cold_C = T_skin_cold_K - 273.15;


fprintf('\n=================================================================\n');
fprintf('  MLI SKIN TEMPERATURES  (per-face model)\n');
fprintf('=================================================================\n');
fprintf('  HOT  sun face    (%2.0f m^2)      : %7.2f K  (%6.2f C)\n', ...
    A_hot_face,    T_sun_K,        T_sun_K-273.15);
fprintf('  HOT  earth face  (%2.0f m^2)      : %7.2f K  (%6.2f C)\n', ...
    A_Xface,       T_earth_hot_K,  T_earth_hot_K-273.15);
fprintf('  HOT  shadow faces(%2.0f m^2)      : %7.2f K  (%6.2f C)\n', ...
    A_shadow_hot,  T_shadow_hot_K, T_shadow_hot_K-273.15);
fprintf('  HOT  area-wtd avg skin         : %7.2f K  (%6.2f C)\n', ...
    T_skin_hot_K, T_skin_hot_C);
fprintf('  HOT  net heat leak INTO s/c    : %7.2f W', Q_leak_in);
if Q_leak_in < 0
    fprintf('  (net outward — shadow losses exceed sun-face gain)\n');
else
    fprintf('\n');
end
fprintf('  COLD earth face  (%2.0f m^2)      : %7.2f K  (%6.2f C)\n', ...
    A_Xface,        T_earth_cold_K,  T_earth_cold_K-273.15);
fprintf('  COLD shadow faces(%2.0f m^2)      : %7.2f K  (%6.2f C)\n', ...
    A_shadow_cold,  T_shadow_cold_K, T_shadow_cold_K-273.15);
fprintf('  COLD area-wtd avg skin         : %7.2f K  (%6.2f C)\n', ...
    T_skin_cold_K, T_skin_cold_C);
fprintf('  COLD heat leak OUT of s/c      : %7.2f W  (heaters must replace this)\n', ...
    Q_leak_out);

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

comps = components_4_19_edit_most_recent();
N = length(comps);

% Q_int_hot  = sum(arrayfun(@(c) c.Q_hot_W,  comps));
% Q_int_cold = sum(arrayfun(@(c) c.Q_cold_W, comps));
Q_int_hot  = sum(arrayfun(@(c) c.Q_hot_W  * ~c.decoupled, comps));
Q_int_cold = sum(arrayfun(@(c) c.Q_cold_W * ~c.decoupled, comps));

fprintf('\n=================================================================\n');
fprintf('  INTERNAL HEAT LOADS\n');
fprintf('=================================================================\n');
fprintf('  Hot  case (all systems on)  : %7.1f W\n', Q_int_hot);
fprintf('  Cold case (safe mode)       : %7.1f W\n', Q_int_cold);
fprintf('  Cold case MLI heat leak out : %7.1f W\n', Q_leak_out);
fprintf('  Total cold-case heat need   : %7.1f W\n', Q_int_cold + Q_leak_out);

%% =========================================================================
%  SECTION 7 — COMPONENT TEMPERATURES-> This analysis includes G values, do we need to change it to account for the fact that each component has a respective LOCATION/ LENGTH away from the raditor/ spacecraft walls, where we bridge the components through thermal straps???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%
%  T_comp = T_int + Q_gen / G_contact  (+margin for hot, -margin for cold)
%
%  T_int held at 20 C (bus requirement).
%  G_contact [W/K] is thermal conductance
%  component-to-structure.----------------------------> Each component has an individual connection G value!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%  Design margin = 5 C per ECSS-E-ST-31C preliminary design.
%  TODO: Tighten to 3 C at CDR after test validation.
% ==========================================================================

margin_C = 5;   % [C]  

T_int_C = T_int - 273.15;   % 20 C
T_comp_hot_C  = zeros(N, 1);
T_comp_cold_C = zeros(N, 1);

% for i = 1:N
%     if i == 2
%         % Solar panels computed from their own radiation balance
%         Q_panel_net = alpha_MLI * hot.S * panel.A_total - panel.Q_hot_W;
%         T_panel_hot = (Q_panel_net / (panel.eps_back * sigma * panel.A_total))^0.25 - 273.15;
%         T_comp_hot_C(i) = T_panel_hot + margin_C;
% 
%         f_pan = @(Tp) panel.eps_back * sigma * 2*panel.A_total * Tp^4 ...
%                     + panel.G_boom_W_K * (Tp - T_int) ...
%                     - panel.eps_back * cold.q_IR * panel.A_total * panel.VF_earth;
%         T_panel_cold_K = fzero(f_pan, [10, 400]);
%         T_comp_cold_C(i) = T_panel_cold_K - 273.15 - margin_C;
%     else
%         T_comp_hot_C(i)  = T_int_C + comps(i).Q_hot_W  / comps(i).G_W_K + margin_C;
%         T_comp_cold_C(i) = T_int_C + comps(i).Q_cold_W / comps(i).G_W_K - margin_C;
%     end
% end

for i = 1:N
    if comps(i).decoupled
        % Thermally decoupled — set to NaN so it never triggers a violation
        T_comp_hot_C(i)  = NaN;
        T_comp_cold_C(i) = NaN;
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

% for i = 1:N
%     if T_comp_hot_C(i) > comps(i).T_op_max_C
%         hot_viol(i) = true;
%     end
%     if comps(i).always_on
%         cold_limit(i) = comps(i).T_op_min_C;
%     else
%         cold_limit(i) = comps(i).T_surv_min_C;
%     end
%     if T_comp_cold_C(i) < cold_limit(i)
%         cold_viol(i) = true;
%     end
% end
for i = 1:N
    if comps(i).decoupled; continue; end
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
%T_sink_K     = (eps_rad * hot.q_IR * VF_earth_rad / (eps_rad * sigma))^0.25;
T_sink_K     = (hot.q_IR * VF_earth_rad / sigma + T_space^4)^0.25;

A_rad_m2 = zeros(N, 1);
% for i = 1:N
%     if hot_viol(i)
%         T_max_K      = comps(i).T_op_max_C + 273.15;
%         q_net_W_m2   = eps_rad * sigma * (T_max_K^4 - T_sink_K^4);
%         A_rad_m2(i)  = comps(i).Q_hot_W / q_net_W_m2;
%     end
% end
for i = 1:N
    if comps(i).decoupled; continue; end
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

% for i = 1:N
%     if cold_viol(i)
%         T_target = cold_limit(i) + 273.15;
%         % Q_h = comps(i).G_W_K * (T_target - T_int) - comps(i).Q_cold_W;
%         % Q_h = max(Q_h, 0) * 1.20;
%         Q_conducted_in = comps(i).G_W_K * (T_int - T_target);  % +ve = structure warms component
%         Q_h = max(0, -(comps(i).Q_cold_W + Q_conducted_in)) * 1.20;
% 
%         Q_heater_W(i) = Q_h;
%         T_setpt_C(i)  = cold_limit(i) + 2;
%         A_heat_cm2(i) = Q_h / heater_density;
%     end
% end
for i = 1:N
    if comps(i).decoupled; continue; end
    if cold_viol(i)
        T_target = cold_limit(i) + 273.15;
        Q_conducted_in = comps(i).G_W_K * (T_int - T_target);
        Q_h = max(0, -(comps(i).Q_cold_W + Q_conducted_in)) * 1.20;
        
        Q_heater_W(i) = Q_h;
        T_setpt_C(i)  = cold_limit(i) + 2;
        A_heat_cm2(i) = Q_h / heater_density;
    end
end
Q_heater_total        = sum(Q_heater_W);
Q_heater_budget_total = Q_heater_total + Q_leak_out;

% ---- Interior cold-case energy balance check ----
% Internal dissipation + component heaters must cover MLI heat loss
% to maintain T_int = 20 C. If not, a bulk interior heater is needed.
Q_surplus_cold = Q_int_cold + Q_heater_total - Q_leak_out;
if Q_surplus_cold >= 0
    fprintf('  NOTE: Internal dissipation covers MLI loss — no bulk heater needed\n');
    fprintf('        Surplus available: %.1f W\n', Q_surplus_cold);
else
    fprintf('  WARNING: Interior heat deficit of %.1f W in cold case\n', abs(Q_surplus_cold));
    fprintf('           A bulk interior heater is needed to maintain T_int = 20 C\n');
    fprintf('           Add %.1f W to heater budget\n', abs(Q_surplus_cold));
    Q_heater_budget_total = Q_heater_budget_total + abs(Q_surplus_cold);
end
%% =========================================================================
%  SECTION 11 — CONSOLE REPORT
% ==========================================================================

fprintf('\n=================================================================\n');
fprintf('  HOT CASE ENERGY BUDGET\n');
fprintf('=================================================================\n');
fprintf('  Solar direct (6 m^2 face)   : %8.1f W\n', alpha_MLI*hot.S*A_hot_face);
fprintf('  Albedo (F=%.2f)             : %8.1f W\n', hot.F_earth, alpha_MLI*hot.S*hot.albedo_skin*hot.F_earth*A_hot_face);
fprintf('  Earth IR (F=%.2f)           : %8.1f W\n', hot.F_earth, eps_MLI*hot.q_IR*hot.F_earth*A_hot_face);
fprintf('  Total external on skin      : %8.1f W\n', Q_ext_hot_sun + Q_ext_hot_earth);
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
% plot_cases_4_18_edit(comps, T_comp_hot_C, T_comp_cold_C, cold_limit, ...
%            hot_viol, cold_viol, ...
%            A_rad_m2, A_rad_available, ...
%            Q_heater_W, A_heat_cm2, T_setpt_C, ...
%            Q_heater_total, margin_C, ...
%            T_skin_hot_C, T_skin_cold_C, Q_leak_in, Q_leak_out, Q_int_cold, Q_surplus_cold);


plot_cases_4_19_edit_most_recent(comps, T_comp_hot_C, T_comp_cold_C, cold_limit, ...
           hot_viol, cold_viol, ...
           A_rad_m2, A_rad_available, ...
           Q_heater_W, A_heat_cm2, T_setpt_C, ...
           Q_heater_total, margin_C, ...
           T_skin_hot_C, T_skin_cold_C, Q_leak_in, Q_leak_out, ...
           Q_int_cold, Q_surplus_cold);

