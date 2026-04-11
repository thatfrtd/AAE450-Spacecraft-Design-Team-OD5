function radiator = radiator_sizing(results_hot, components, surf, env)
% RADIATOR_SIZING  Size radiators to reject waste heat in worst-case hot scenario
%
%  DESCRIPTION:
%    Computes required radiator area to maintain all components below
%    their maximum operating temperatures in the hot case.
%
%    Radiator sizing equation (steady-state energy balance):
%      Q_reject = eps_rad * sigma * A_rad * (T_rad^4 - T_sink^4)
%
%    where:
%      Q_reject = total heat to be rejected [W]
%      eps_rad  = radiator emissivity (goal: > 0.85 for good rejection)
%      A_rad    = required radiator area [m^2]
%      T_rad    = radiator operating temperature [K]
%      T_sink   = effective sink temperature [K]
%
%    T_sink (effective) accounts for:
%      - Deep space (2.7 K, nearly zero)
%      - Earth IR absorbed by radiator
%      - Albedo absorbed by radiator
%
%    Radiator T is set to the most demanding (hottest allowable) component.
%
%  STRATEGY:
%    1. Body-mounted side radiators (primary — simplest, no deployment)
%    2. Panel back-surface radiation (secondary)
%    3. If insufficient → flag need for deployable radiators
%
%  OUTPUTS:
%    radiator.area_required_m2  — minimum radiator area [m^2]
%    radiator.area_available_m2 — available body area [m^2]
%    radiator.body_sufficient   — boolean
%    radiator.T_rad_K           — assumed radiator operating temperature [K]
%    radiator.Q_reject_W        — heat load to be rejected [W]
%    radiator.Q_per_m2          — W/m^2 specific rejection capacity
%
%  NOTE: This is a steady-state sizing. The transient model verifies that
%        peak temperatures remain bounded with this area.
%
%  TODO: Once attitude profile is fixed, compute actual solar and albedo
%        load on radiator (currently uses conservative estimate).
%  TODO: If deployable radiators are needed, size boom length and panel dims.
% =========================================================================

sigma    = 5.6704e-8;   % W/m^2/K^4
T_space  = 2.7;         % K

%% ---- Determine heat rejection requirement ------------------------------
%  Total internal power dissipation in hot case
%  (all components running at full power)
Q_internal = 0;
for i = 1:length(components)
    Q_internal = Q_internal + components(i).P_nom_W;
end

%  In the hot case, we also absorb external fluxes (already handled by
%  the transient — here we size for the internal dissipation plus a margin)
%  Design margin: 20% on power dissipation                % TODO: confirm margin with team
design_margin = 1.20;
Q_reject = Q_internal * design_margin;

radiator.Q_reject_W = Q_reject;
fprintf('  [Radiator] Required heat rejection: %.1f W (%.1f W + %.0f%% margin)\n', ...
    Q_reject, Q_internal, (design_margin-1)*100);

%% ---- Radiator operating temperature ------------------------------------
%  Set radiator temp to the most restrictive component's max operating temp
%  minus a ΔT gradient for conduction to radiator.
%  ΔT_conduction: assume ~10°C drop from component to radiator surface % TODO: model with strap/pipe
delta_T_conduction = 10;   % °C  % TODO: recalculate with strap sizing

T_max_allowed_C = min(arrayfun(@(c) c.T_op_max_C, components));
T_rad_C  = T_max_allowed_C - delta_T_conduction;
T_rad_K  = T_rad_C + 273.15;

radiator.T_rad_K  = T_rad_K;
radiator.T_rad_C  = T_rad_C;
fprintf('  [Radiator] Radiator operating temp: %.1f°C (%.1f K)\n', T_rad_C, T_rad_K);

%% ---- Effective sink temperature ----------------------------------------
%  Radiators face space, but also absorb some Earth IR and albedo.
%  Net effective sink temperature for the side radiators:
%
%  q_absorbed = alpha_rad * (albedo + solar_reflected) + eps_rad * q_EarthIR * VF_earth
%  Q_net = eps_rad * sigma * A * T_rad^4 - q_absorbed * A - eps_rad * sigma * T_space^4 * A
%
%  Simplified: lump Earth IR absorption into effective T_sink
%  TODO: Compute properly with actual radiator orientation and scenario

eps_rad   = 0.85;   % radiator emissivity (white paint or OSR)          % TODO: confirm coating
alpha_rad = 0.20;   % radiator absorptivity (low for white paint / OSR)  % TODO: confirm

q_EarthIR_on_radiator = env.q_earth_IR * 0.3;  % ~30% of Earth IR flux hits side radiators % TODO
q_solar_leakage       = 0;                      % side radiators don't face Sun in nadir mode

% Effective power absorbed per unit area of radiator
q_absorbed_per_m2 = alpha_rad * env.albedo_earth * 1361 * 0.1 + ...  % albedo leakage (10%)
                    eps_rad * q_EarthIR_on_radiator;

% Net radiative rejection per m^2:
Q_per_m2 = eps_rad * sigma * (T_rad_K^4 - T_space^4) - q_absorbed_per_m2;

radiator.Q_per_m2 = Q_per_m2;
fprintf('  [Radiator] Net rejection capacity: %.1f W/m^2\n', Q_per_m2);

%% ---- Required radiator area --------------------------------------------
A_required = Q_reject / Q_per_m2;
radiator.area_required_m2 = A_required;
fprintf('  [Radiator] Required area: %.3f m^2\n', A_required);

%% ---- Available body area (from geometry) --------------------------------
% Side area available (after RCS blockage)
A_side_avail  = pi * 3.0 * 3.2 * 0.70;   % ~21.1 m^2
% Panel backs available
A_panel_backs = 2 * 3.2 * 9.0;           % 57.6 m^2 (but panels face space)

A_available = A_side_avail;   % conservative: only use side body for now
radiator.area_available_m2 = A_available;
radiator.A_side_avail      = A_side_avail;
radiator.A_panel_backs     = A_panel_backs;

%% ---- Decision ----------------------------------------------------------
if A_required <= A_side_avail
    radiator.body_sufficient = true;
    radiator.coverage_pct    = (A_required / A_side_avail) * 100;
    fprintf('  [Radiator] Body-mounted SUFFICIENT (%.1f%% of side area used)\n', ...
        radiator.coverage_pct);
elseif A_required <= A_side_avail + A_panel_backs * 0.5
    radiator.body_sufficient = false;
    radiator.need_panel_backs = true;
    fprintf('  [Radiator] Need panel back surfaces as additional radiators.\n');
else
    radiator.body_sufficient  = false;
    radiator.need_deployable  = true;
    fprintf('  [Radiator] WARNING: Deployable radiators likely required!\n');
    fprintf('             Shortfall: %.2f m^2\n', A_required - A_side_avail);
end

%% ---- Summary by component ----------------------------------------------
fprintf('\n  [Radiator] Per-component heat rejection requirements:\n');
for i = 1:length(components)
    c = components(i);
    fprintf('    %-25s  P=%.0f W  T_max=%.0f°C\n', c.name, c.P_nom_W, c.T_op_max_C);
end

end
