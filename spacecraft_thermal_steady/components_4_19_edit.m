%% Change History- 4/19
%  CHANGE SUMMARY
%
%  CHANGE 1 — Ion Thruster Q_hot_W corrected
%    BEFORE: comps(i).Q_hot_W = 0.30 * 6900   (2070 W — total thruster waste heat)
%    AFTER:  comps(i).Q_hot_W = 0.05 * 0.30 * 6900  (~103 W — heat into structure only)
%    REASON: The thruster radiates ~95% of its waste heat directly to space
%            from the discharge chamber and bell. Only the fraction conducted
%            through the mounting flange into the spacecraft structure belongs
%            in this model. Using 2070 W produced a physically unrealistic
%            radiator requirement. The 5% figure is a placeholder — get exact
%            value from Aerojet NEXT-C ICD.
%
%  CHANGE 2 — Added .decoupled field to all components
%    BEFORE: Field did not exist.
%    AFTER:  .decoupled = false for all components except Solar Panels.
%            .decoupled = true  for Solar Panels (component 2).
%    REASON: Solar panels are thermally decoupled from the spacecraft body.
%            This flag allows main_thermal.m to exclude the panel from all
%            heat sums, violation checks, radiator sizing, and heater sizing
%            loops using a single consistent flag rather than a hardcoded
%            component index.

function comps = components_4_19_edit()
% COMPONENTS  All spacecraft thermal component definitions
% ==========================================================================
%
%  HOW TO USE THIS FILE
%  --------------------
%  This is the primary file to update as your team provides real data.
%  Each block below defines one component. To add a new component:
%    1. Increment the counter:  i = i + 1;
%    2. Fill in all fields (see field descriptions at the top of each block)
%    3. Run main_thermal.m
%
%  FIELD DESCRIPTIONS
%  ------------------
%  .name          String label (used in plots and reports)
%  .Q_hot_W       Heat dissipated [W] in worst-case hot (all systems on)
%  .Q_cold_W      Heat dissipated [W] in worst-case cold / safe mode
%  .T_op_min_C    Min operating temperature [°C]
%  .T_op_max_C    Max operating temperature [°C]
%  .T_surv_min_C  Minimum SURVIVAL temperature [°C] (checked in cold case)
%  .T_surv_max_C  Maximum survival temperature [°C]
%  .G_W_K         Thermal conductance to spacecraft structure [W/K]
%                 (see guidelines below — replace all values with TODO comment)
%  .always_on     true  = component operates in cold/safe mode (check T_op_min)
%                 false = component is off in cold/safe mode (check T_surv_min)
%
%  CONDUCTANCE GUIDELINES  (G, thermal contact to structure)
%  ----------------------------------------------------------
%  These are estimates. Replace with calculated or measured values.
%    Large bolted flange with TIM (thermal interface material): 3–10 W/K
%    Standard bolted bracket, Al-to-Al:                        1–3  W/K
%    Small exterior bracket, limited contact:                  0.3–1 W/K
%    Hinged deployment boom:                                   3–8  W/K
%    Press-fit or potted in place:                             5–15 W/K
%
%  COLD-CASE POWER NOTE
%  --------------------
%  Q_cold_W = dissipation in eclipse safe mode (not operating power).
%  For components that are OFF in safe mode, Q_cold_W = 0.
%  For components that are ON, use their safe-mode (standby/minimal) power.
%
% ==========================================================================

sigma = 5.6704e-8;   % included so caller can use if needed
i = 0;

%% -------------------------------------------------------------------------
%  1. NEXT-C Ion Thruster
%     Rear of spacecraft. Fires during normal ops; OFF in safe mode.
%     Only ~30% of input power becomes waste heat (rest = thrust + ionisation).
%  TODO: Confirm thermal dissipation fraction from Aerojet Rocketdyne ICD.
%  TODO: Confirm operating and survival temperature limits from NEXT-C ICD.
%  TODO: Confirm whether thruster fires continuously or duty-cycled.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Ion Thruster (NEXT-C)';
comps(i).Q_hot_W      = 0.05* 0.30 * 6900;   % 30% of 6.9 kW input  % TODO: verify efficiency % 4_18 edit--> estimates only a small percentage(~5%) of this heat affects the rest of the structure, the rest of the heat is radiated out to space 
comps(i).Q_cold_W     = 0;             % OFF in safe mode
comps(i).T_op_min_C   = -20;           % TODO: from NEXT-C ICD
comps(i).T_op_max_C   = 300;           % TODO: from NEXT-C ICD (interface temp, not discharge)
comps(i).T_surv_min_C = -30;           % TODO: from NEXT-C ICD
comps(i).T_surv_max_C = 350;           % TODO: from NEXT-C ICD
comps(i).G_W_K        = 4.0;           % large mounting flange  % TODO: calculate from interface area
comps(i).always_on    = false;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  2. Solar Panels  (both panels as one lumped node)
%     Wing-mounted, sun-tracking. Generates 3.624 kW waste heat when illuminated.
%     In eclipse: OFF electrically, acts as cold radiator attached to body.
%  NOTE: Panel temperature is calculated SEPARATELY in main_thermal.m (not via G).
%        This entry is for the panel thermal loads on the spacecraft body only.
%  TODO: Confirm 3.624 kW heat figure is waste heat (not electrical output).
%  TODO: Confirm panel survival temperature with vendor (critical for cold case).
%  TODO: Confirm back-surface emissivity and coating from panel vendor.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Solar Panels (2x)';
comps(i).Q_hot_W      = 3624;     % [W] waste heat into structure  % team value
comps(i).Q_cold_W     = 0;        % eclipse: no solar generation
comps(i).T_op_min_C   = 18;       % 28°C optimal - 10°C margin  % TODO: from vendor
comps(i).T_op_max_C   = 140;      % team confirmed
comps(i).T_surv_min_C = -150;     % GaAs cells cold-tolerant  % TODO: CRITICAL — confirm from vendor
comps(i).T_surv_max_C = 160;      % TODO: confirm from vendor
comps(i).G_W_K        = 5.0;      % hinged deployment boom  % TODO: from mech team
comps(i).always_on    = false;    % off in cold case (eclipse)
comps(i).decoupled    = true;
%% -------------------------------------------------------------------------
%  3. Battery
%     Interior, preferably near solar panels for thermal coupling.
%     130 W heat generation (from charging/discharging losses, team value).
%     Li-ion is temperature-sensitive — this component is a key thermal driver.
%  TODO: Confirm cold-case dissipation with EPS team (battery discharging in eclipse).
%  TODO: Confirm thermal interface material and mounting plate conductance.
%  TODO: Confirm operating range is 10–40°C nominal (some Li-ion = 0–45°C).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Battery';
comps(i).Q_hot_W      = 130;      % team value
comps(i).Q_cold_W     = 40;       % discharging to power safe-mode loads  % TODO: EPS team
comps(i).T_op_min_C   = 10;       % team value (narrow range — key design driver)
comps(i).T_op_max_C   = 40;       % team value
comps(i).T_surv_min_C = 0;        % Li-ion survival  % TODO: confirm chemistry/vendor
comps(i).T_surv_max_C = 60;       % TODO: confirm
comps(i).G_W_K        = 3.0;      % internal shelf mount with TIM  % TODO: mech team
comps(i).always_on    = true;     % must operate in cold case (powers safe-mode loads)
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  4. Camera
%     Front (nadir/forward facing). Off in safe mode.
%  TODO: Confirm survival temperature range with payload team.
%  TODO: Confirm whether 20 W is peak or average (use peak = worst case).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Camera';
comps(i).Q_hot_W      = 20;       % peak, team value
comps(i).Q_cold_W     = 0;        % OFF in safe mode
comps(i).T_op_min_C   = 0;        % team value
comps(i).T_op_max_C   = 40;       % team value
comps(i).T_surv_min_C = -20;      % TODO: confirm from payload team
comps(i).T_surv_max_C = 60;       % TODO: confirm
comps(i).G_W_K        = 1.0;      % exterior bracket  % TODO: mech team
comps(i).always_on    = false;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  5. Star Tracker
%     Exterior, away from thrusters, opposite side from antennas.
%     Always on — attitude determination needed in all modes.
%     Use max power (8.9 W) as conservative estimate.
%  TODO: Confirm continuous vs intermittent operation with ADCS team.
%  TODO: Confirm survival temp from vendor datasheet.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Star Tracker';
comps(i).Q_hot_W      = 8.9;      % max of 5.6–8.9 W range
comps(i).Q_cold_W     = 8.9;      % always on
comps(i).T_op_min_C   = -30;      % team value
comps(i).T_op_max_C   = 60;       % team value
comps(i).T_surv_min_C = -40;      % TODO: confirm from vendor
comps(i).T_surv_max_C = 70;       % TODO: confirm
comps(i).G_W_K        = 1.0;      % exterior bracket  % TODO
comps(i).always_on    = true;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  6. Quadrifilar Helix Antennas  (both antennas as one node)
%     Opposite ends of spacecraft. Always on (TT&C needed in all modes).
%     10 W total power draw; only 1.5 W becomes heat (rest transmitted as RF).
%  TODO: Confirm heat-to-RF split with TT&C team.
%  TODO: Confirm survival temperature from vendor.
%  NOTE: mass = 2 × 6 lb = 2 × 2.72 kg = 5.44 kg total
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'QH Antennas (2x)';
comps(i).Q_hot_W      = 1.5;      % heat generated only (10 W - RF output), team value
comps(i).Q_cold_W     = 1.5;      % always on
comps(i).T_op_min_C   = -40;      % team value
comps(i).T_op_max_C   = 80;       % team value
comps(i).T_surv_min_C = -50;      % TODO: confirm from vendor
comps(i).T_surv_max_C = 90;       % TODO: confirm
comps(i).G_W_K        = 0.5;      % exterior mount  % TODO
comps(i).always_on    = true;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  7. GPS Antennas  (both antennas as one node)
%     Opposite ends of spacecraft. Negligible heat generation per team.
%     Use 0.5 W as a conservative non-zero estimate.
%  TODO: Confirm heat generation (team said "negligible" — use 0 if confirmed).
%  TODO: Confirm continuous operation mode with TT&C team.
%  NOTE: mass = 2 × 80 g = 0.16 kg total
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'GPS Antennas (2x)';
comps(i).Q_hot_W      = 0.5;      % negligible per team — conservative 0.5 W placeholder
comps(i).Q_cold_W     = 0.5;      % always on
comps(i).T_op_min_C   = -55;      % team value
comps(i).T_op_max_C   = 85;       % team value
comps(i).T_surv_min_C = -65;      % TODO: confirm from vendor
comps(i).T_surv_max_C = 95;       % TODO: confirm
comps(i).G_W_K        = 0.5;      % small exterior patch  % TODO
comps(i).always_on    = true;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  8. RCS Thrusters  (all 16 modelled as one node)
%     Mounted on sides near ends. Fire intermittently for attitude control.
%     Average power dissipation while firing is uncertain — use conservative estimate.
%  TODO: Get average dissipation per thruster from propulsion team.
%  TODO: Get temperature limits from propulsion team.
%  NOTE: 16 × 0.64 kg = 10.24 kg total
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'RCS Thrusters (16x)';
comps(i).Q_hot_W      = 50;       % TODO: get from propulsion team (very uncertain placeholder)
comps(i).Q_cold_W     = 0;        % OFF in safe mode
comps(i).T_op_min_C   = -20;      % TODO: from propulsion team
comps(i).T_op_max_C   = 120;      % TODO: from propulsion team
comps(i).T_surv_min_C = -40;      % TODO: from propulsion team
comps(i).T_surv_max_C = 150;      % TODO: from propulsion team
comps(i).G_W_K        = 1.0;      % side brackets  % TODO: mech team
comps(i).always_on    = false;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  9. Arm
%     Front of spacecraft. All specifications are placeholders.
%  TODO: Get mass, power, temperature limits from structures/mechanisms team.
%  TODO: Confirm operating mode (does arm operate in safe mode?).
%  TODO: Confirm arm material (affects thermal mass and conductance).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Arm';
comps(i).Q_hot_W      = 10;       % TODO: actuator/motor power (placeholder)
comps(i).Q_cold_W     = 0;        % TODO: assume stowed/OFF in safe mode
comps(i).T_op_min_C   = -40;      % TODO: from structures/mech team
comps(i).T_op_max_C   = 80;       % TODO
comps(i).T_surv_min_C = -55;      % TODO
comps(i).T_surv_max_C = 100;      % TODO
comps(i).G_W_K        = 2.0;      % main body attachment  % TODO
comps(i).always_on    = false;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  10. Plasma Contactor
%      Rear (near thruster). All specifications are placeholders.
%  TODO: Get all specifications — power, temperature limits, mass.
%  TODO: Confirm whether it operates in safe mode.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Plasma Contactor';
comps(i).Q_hot_W      = 10;       % TODO: placeholder — get from subsystem team
comps(i).Q_cold_W     = 0;        % TODO: assume OFF in safe mode
comps(i).T_op_min_C   = -40;      % TODO
comps(i).T_op_max_C   = 80;       % TODO
comps(i).T_surv_min_C = -50;      % TODO
comps(i).T_surv_max_C = 100;      % TODO
comps(i).G_W_K        = 1.0;      % TODO
comps(i).always_on    = false;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  11. CMGs  (if still in design — confirm with ADCS team)
%      NOTE: CMGs were not listed in the latest team component table.
%            Keeping here as a placeholder — DELETE this block if CMGs are
%            removed from the design.
%  TODO: Confirm CMG inclusion with ADCS team.
%  TODO: Get mass, power, temperature limits from CMG vendor.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'CMGs (4x) [CONFIRM w/ ADCS]';
comps(i).Q_hot_W      = 27.5;     % midpoint of 20–35 W estimate  % TODO
comps(i).Q_cold_W     = 0;        % assume standby/OFF in safe mode  % TODO
comps(i).T_op_min_C   = -20;      % TODO: from vendor
comps(i).T_op_max_C   = 70;       % TODO: from vendor
comps(i).T_surv_min_C = -30;      % TODO
comps(i).T_surv_max_C = 85;       % TODO
comps(i).G_W_K        = 3.0;      % rigid interior mount  % TODO
comps(i).always_on    = false;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  12. Flight Computer
%      NOTE: Not listed in latest team table — confirm with avionics team.
%            Always on (processes all safe-mode commands).
%  TODO: Confirm model and get exact power specs from avionics team.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Flight Computer [CONFIRM]';
comps(i).Q_hot_W      = 5;        % TODO: placeholder
comps(i).Q_cold_W     = 5;        % always on
comps(i).T_op_min_C   = -20;      % TODO: from vendor
comps(i).T_op_max_C   = 70;       % TODO: from vendor
comps(i).T_surv_min_C = -30;      % TODO
comps(i).T_surv_max_C = 85;       % TODO
comps(i).G_W_K        = 2.0;      % interior PCB mount  % TODO
comps(i).always_on    = true;
comps(i).decoupled    = false;
%% -------------------------------------------------------------------------
%  13. IMU
%      NOTE: Not listed in latest team table — confirm with ADCS team.
%            Always on (attitude knowledge needed in safe mode).
%  TODO: Confirm model and get exact specs from ADCS team.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'IMU [CONFIRM]';
comps(i).Q_hot_W      = 12;       % TODO: placeholder from earlier estimate
comps(i).Q_cold_W     = 12;       % always on
comps(i).T_op_min_C   = -40;      % TODO: from vendor (IMUs often rated to -40°C)
comps(i).T_op_max_C   = 85;       % TODO: from vendor
comps(i).T_surv_min_C = -55;      % TODO
comps(i).T_surv_max_C = 95;       % TODO
comps(i).G_W_K        = 2.0;      % interior rigid mount  % TODO
comps(i).always_on    = true;
comps(i).decoupled    = false;
end
