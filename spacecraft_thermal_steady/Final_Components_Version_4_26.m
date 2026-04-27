function comps = Final_Components_Version_4_26()
% COMPONENTS  All spacecraft thermal component definitions
% ==========================================================================
%
%  HOW TO ADD A NEW COMPONENT
%  --------------------------
%  1. Add  i = i + 1;
%  2. Copy any existing block and fill in values.
%  3. Run AAE450_Thermal_Main.m.
%
%  FIELD DESCRIPTIONS
%  ------------------
%  .name          String label (used in plots and console reports)
%  .Q_hot_W       Heat dissipated [W] in worst-case hot (all systems on)
%  .Q_cold_W      Heat dissipated [W] in worst-case cold / safe mode
%  .T_op_min_C    Min operating temperature [C]
%  .T_op_max_C    Max operating temperature [C]
%  .T_surv_min_C  Min SURVIVAL temperature [C]  (checked in cold case if always_on=false)
%  .T_surv_max_C  Max survival temperature [C]
%  .G_W_K         Thermal conductance to spacecraft thermal bus [W/K]
%  .always_on     true  = on in eclipse (T_op_min is the cold-case limit)
%                 false = off in eclipse (T_surv_min is the cold-case limit)
%  .decoupled     true  = solar panels only — temperature solved from own
%                         radiation balance in main, NOT from T_int +/- Q/G
%                 false = standard body component
%
%  THERMAL BUS CONCEPT
%  --------------------
%  G_W_K is the conductance from the component to the common THERMAL BUS
%  (the structural cold plate connected to the radiator).
%  Component temperature: T_comp = T_bus + Q_hot / G_W_K  (hot case)
%                         T_comp = T_bus + Q_cold / G_W_K  (cold case)
%  High G is GOOD for hot-case cooling (keeps T_comp near T_bus).
%  High G is BAD for cold-case heating (drains heat faster).
%  In the cold case the radiator is thermally switched off, so G connects
%  the component to the interior structure temperature only.
%
%  G_W_K DESIGN GUIDELINES
%  -------------------------
%  Cu cold plate with TIM, large contact area  :  10 – 30 W/K
%  Al cold plate with TIM, moderate area       :   5 – 15 W/K
%  Large bolted flange with TIM                :   3 – 10 W/K
%  Standard bolted bracket, Al-to-Al           :   1 –  3 W/K
%  Small exterior bracket, limited contact     : 0.3 –  1 W/K
%  Hinged deployment boom                      :   3 –  8 W/K
%  Titanium strut mount (thermally isolating)  : 0.5 –  2 W/K
%
%  REVISION LOG
%  ------------
%  Rev 1  Original model — single thruster node, bulk heater, fixed T_int.
%  Rev 2  Thruster split into Head + PPU (PPU was entirely missing).
%         T_int solved; component-level heaters; single shared radiator.
%  Rev 3  Geometry corrected to 2x2x2 m.
%  Rev 4  (THIS FILE) — Final consolidated version:
%           - Xenon Tank and Chemical Tanks added from earlier version.
%           - Star Tracker corrected to always_on = true.
%           - Battery G increased to 15 W/K (requires dedicated cold plate).
%           - CMG, Flight Computer, IMU G values tightened.
%           - PPU Q_hot confirmed at 621 W; T_op_max = 70 C.
%           - RCS Q_hot kept conservative at 50 W pending propulsion data.
%           - strap_length / calc_G approach removed (adds false precision
%             without interface resistance data; use G_W_K directly).
%           - All TODO flags updated.
%
% ==========================================================================

i = 0;

%% =========================================================================
%  PROPULSION SYSTEM
%  Ion thruster system has two SEPARATE thermal nodes:
%    (1) Thruster Head  — exterior, structural heat soak only (~103 W)
%    (2) PPU            — interior electronics, full waste heat (~621 W)
%  Do NOT combine these. The PPU is the dominant hot-case heat source.
% =========================================================================

%% -------------------------------------------------------------------------
%  1. NEXT-C Ion Thruster HEAD
%
%  Represents heat that travels from the discharge bell through the gimbal
%  flange into the spacecraft structure. The discharge chamber itself
%  radiates ~70% of its waste heat directly to space.
%
%  Basis:  efficiency ~70%  =>  30% of 6900 W = 2070 W waste at head.
%          Structural soakback fraction ~5% (gimbal/strut isolation).
%          Q_struct = 0.05 * 2070 = 103.5 W.
%
%  TODO: Verify soakback fraction (5%) from Aerojet Rocketdyne NEXT-C ICD.
%  TODO: Confirm T_op limits from NEXT-C ICD (interface, not discharge temp).
%  TODO: Recalculate G_W_K from actual gimbal flange area, TIM, and ICD.
%        Note: low G is intentional to isolate the hot thruster from structure.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Ion Thruster Head';
comps(i).Q_hot_W      = 0.05 * 0.30 * 6900;   % 103.5 W — structural interface only
comps(i).Q_cold_W     = 0;
comps(i).T_op_min_C   = -24;    % Team component layout PDF (Op5 OD5, April 2026)
comps(i).T_op_max_C   = 150;    % 300, Placeholder value
% Actually ~300 C, artificially set to 150 C to appear on graph                                 
                                 % NOTE: at G=4 W/K, ΔT=103.5/4=25.9°C → T_bus_max=19.1°C
                                 % This may be the binding hot-case constraint; see console.
comps(i).T_surv_min_C = -40;    % Team component layout PDF
comps(i).T_surv_max_C = 350;     % Team component layout PDF
comps(i).G_W_K        = 4.0;    % Large mounting flange — low G intentional for isolation
                                 % TODO: calculate from gimbal geometry and TIM spec
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  2. NEXT-C Ion Thruster PPU  (Power Processing Unit)
%
%  CRITICAL: This is the largest single interior heat source (621 W hot).
%  The PPU is a separate interior chassis converting bus power to the
%  discharge/neutralizer voltages. It was missing from earlier model versions.
%
%  Basis:  ~91% PPU efficiency  =>  Q_PPU = (1 - 0.91) * 6900 = 621 W.
%          All 621 W dissipated INSIDE the spacecraft (unlike thruster head).
%
%  With G = 25 W/K:  ΔT = 621/25 = 24.8 C  — borderline, requires cold plate.
%  With G = 15 W/K:  ΔT = 621/15 = 41.4 C  — bus must be < 28.6 C (tight).
%  25 W/K requires a direct-bonded Cu cold plate or embedded heat pipe.
%
%  TODO: Confirm efficiency and Q from Aerojet Rocketdyne / PPU vendor.
%  TODO: Confirm T_op_max from PPU datasheet (70 C assumed typical electronics).
%  TODO: Design cold plate and verify G_W_K achievability.
%  TODO: Confirm PPU is interior-mounted (vs external).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Ion Thruster PPU';
comps(i).Q_hot_W      = (1 - 0.92) * 6900;   % 552 W — 92% efficiency (Brophy et al. 2012,
                                               % "NEXT Ion Propulsion System," AIAA-2012-3833:
                                               % NEXT-class PPU measured at ~92-93% efficiency)
comps(i).Q_cold_W     = 0;                    % thruster off in safe mode
comps(i).T_op_min_C   = -40;                  % Checked
comps(i).T_op_max_C   = 50;                   % Checked
comps(i).T_surv_min_C = -24;                  % Checked
comps(i).T_surv_max_C = 71;                   % Checked
comps(i).G_W_K        = 25.0;   % DESIGN TARGET — requires Cu cold plate or direct strap
                                 % TODO: calculate from cold plate geometry and TIM spec
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  3. Xenon Propellant Tank
%
%  Stores high-pressure xenon for the ion thruster. No self-heating.
%  Thermal limits driven by liner material and pressure vessel design.
%  Upper limit (50 C) is a common conservative bound for composite-lined
%  COPV tanks; confirmed with propulsion team.
%  Mounted on isolated titanium struts to limit structural heat transfer.
%
%  TODO: Confirm T_op limits from tank vendor / propulsion team.
%  TODO: Confirm G_W_K from actual strut geometry.
%  TODO: Confirm whether tank needs a heater to prevent propellant freezing
%        (Xe condenses at -108 C; well below survival limit, so likely OK).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Xenon Tank';
comps(i).Q_hot_W      = 0;
comps(i).Q_cold_W     = 0;
comps(i).T_op_min_C   = 16.1;   % Team component layout PDF (OD5, April 2026)
                                  % WARNING: high T_op_min means tank WILL need heaters in eclipse.
                                  % At G=2 W/K, heater ≈ 2*(16.1 - T_int_cold) W before margin.
comps(i).T_op_max_C   = 50;     % Team component layout PDF (pressure vessel limit)
comps(i).T_surv_min_C = 16.1;     % Margin below op min; above any regulator/valve cold limit
                                  % Xe condenses at -108°C @1atm — not the constraint here
comps(i).T_surv_max_C = 50;     % Team component layout PDF
comps(i).G_W_K        = 2.0;    % Isolated titanium strut mounts
                                 % TODO: calculate from strut geometry
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  4. Chemical Propellant Tanks  (RCS — MMH Fuel + MON-3 Oxidizer)
%
%  Two separate spherical tanks confirmed by team PDF (OD5, April 2026):
%    MMH  (fuel):      D=0.713 m, wet mass=164.48 kg, T_op: -2 to 70°C
%    MON-3 (oxidizer): D=0.715 m, wet mass=267.99 kg, T_op: -4 to 40°C
%  (MON-3 = NTO + 3% NO; freeze point ≈ -11.2°C, similar to pure NTO)
%
%  Treated as one thermal node (similar mass, similar mounting, close proximity).
%  Binding limits for single-node model:
%    T_op_min = -2°C  (MMH lower limit; stricter than MON-3's -4°C)
%    T_op_max = 40°C  (MON-3 upper limit; stricter than MMH's 70°C)
%    T_surv_min = -9°C (3°C margin above MON-3 freeze point of -11.2°C)
%    T_surv_max = 50°C (margin above op max)
%
%  NOTE: T_surv_min=-9°C is the binding survival constraint (MON-3 freeze).
%        If MMH/MON-3 ever drops to -9°C a heater must fire. Tank heaters
%        are almost certainly required — G=2 W/K gives essentially zero
%        self-heating, so all cold-case protection comes from heaters.
%
%  Source: Team Component Layout PDF (OD5, April 2026); Propellant freeze
%          points from USAF Propellant MSDS (NTO/MMH) and ECSS-E-HB-31A.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Chemical Tanks (RCS)';
comps(i).Q_hot_W      = 0;
comps(i).Q_cold_W     = 0;
comps(i).T_op_min_C   = -2;     % MMH lower limit (team PDF); binding over MON-3 -4°C
comps(i).T_op_max_C   = 40;     % MON-3 upper limit (team PDF); binding over MMH 70°C
comps(i).T_surv_min_C = -2;     % 3°C margin above MON-3 freeze (-11.2°C) [ECSS-E-HB-31A]
comps(i).T_surv_max_C = 40;     % Margin above op max; confirm with tank vendor
comps(i).G_W_K        = 2.0;    % Isolated strut mounts
                                  % TODO: calculate from strut geometry
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% =========================================================================
%  POWER SYSTEM
% =========================================================================

%% -------------------------------------------------------------------------
%  5. Solar Panels (2x)  [DECOUPLED — own radiation balance]
%
%  Temperature solved separately in main Section 5, NOT from T_int +/- Q/G.
%  Q_hot_W = 0 to prevent double-counting with the Section 5 balance.
%  3624 W is the team-provided electrical extraction value.
%
%  TODO: Clarify with EPS team: is 3624 W electrical output OR waste heat?
%        Interpretation changes the panel temperature calculation.
%  TODO: Confirm T_op_min (currently 18 C = 28 C optimal - 10 C margin).
%  TODO: Confirm T_surv_min with panel vendor (GaAs cells — typically -150 C).
%  TODO: Confirm G_boom with mechanical team (boom hinge design).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Solar Panels (2x)';
comps(i).Q_hot_W      = 0;      % DO NOT CHANGE — handled in main Section 5
comps(i).Q_cold_W     = 0;      % eclipse: no generation
comps(i).T_op_min_C   = 18;     % 28 C optimal - 10 C margin  % TODO: vendor
comps(i).T_op_max_C   = 140;    % team confirmed
comps(i).T_surv_min_C = -150;   % GaAs cold-tolerant  % TODO: confirm with vendor
comps(i).T_surv_max_C = 160;    % TODO: confirm with vendor
comps(i).G_W_K        = 5.0;    % Hinged deployment boom  % TODO: mech team
comps(i).always_on    = false;
comps(i).decoupled    = true;   % computed separately in main Section 5

%% -------------------------------------------------------------------------
%  6. Battery
%
%  Li-ion, interior-mounted. 130 W hot-case dissipation (team value).
%  Narrow temperature window (10–40 C) is the primary hot-case design driver.
%  40 W cold-case is eclipse discharge power to safe-mode loads.
%
%  WHY G = 15 W/K:
%    With G = 3 W/K (old value):  ΔT = 130/3 = 43.3 C  →  bus must reach -3 C
%    With G = 15 W/K (this value): ΔT = 130/15 = 8.7 C  →  bus needs to be ~31 C
%    15 W/K requires a dedicated Al or Cu cold plate with TIM.
%    This is a design REQUIREMENT, not just a passive interface estimate.
%
%  TODO: Confirm 130 W hot-case and 40 W cold-case with EPS team.
%  TODO: Design battery cold plate to achieve G = 15 W/K; verify with mech team.
%  TODO: Confirm survival limits from cell vendor (0 C / 60 C).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Battery';
comps(i).Q_hot_W      = 130;    % team value (charge + discharge losses)
comps(i).Q_cold_W     = 40;     % eclipse discharge to safe-mode loads  % TODO: EPS team
comps(i).T_op_min_C   = 10;     % team value (Li-ion lower charge limit)
comps(i).T_op_max_C   = 40;     % team value (Li-ion upper limit)
comps(i).T_surv_min_C = 10;      % Li-ion survival (no charging below 0 C)
comps(i).T_surv_max_C = 40;     % TODO: confirm with cell vendor
comps(i).G_W_K        = 15.0;   % REQUIRES dedicated cold plate with TIM
                                 % TODO: design cold plate; verify G achievability
comps(i).always_on    = true;
comps(i).decoupled    = false;

%% =========================================================================
%  PAYLOAD
% =========================================================================

%% -------------------------------------------------------------------------
%  7. Camera
%
%  Off in safe mode. G = 3 W/K is achievable with a bolted bracket and TIM.
%  TODO: Confirm 20 W is peak (worst case); use peak for worst-case hot.
%  TODO: Confirm survival temperature with payload team.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Camera';
comps(i).Q_hot_W      = 20;     % peak dissipation  % TODO: payload team
comps(i).Q_cold_W     = 0;
comps(i).T_op_min_C   = 0;
comps(i).T_op_max_C   = 40;
comps(i).T_surv_min_C = 0;      %
comps(i).T_surv_max_C = 40;     % 
comps(i).G_W_K        = 3.0;    % Bolted bracket with TIM
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  8. Robotic Arm
%
%  Off in safe mode (stowed). Actuator waste heat in hot case.
%  TODO: All values are placeholders — get from mechanisms/structures team.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Robotic Arm';
comps(i).Q_hot_W      = 10;     % TODO: actuator power — placeholder
comps(i).Q_cold_W     = 0;      % assumed stowed and off in safe mode
comps(i).T_op_min_C   = -50;    % TODO: mechanisms team
comps(i).T_op_max_C   = 100;     % TODO: mechanisms team
comps(i).T_surv_min_C = -50;    % TODO: mechanisms team
comps(i).T_surv_max_C = 100;    % TODO: mechanisms team
comps(i).G_W_K        = 2.0;    % Main body attachment point
                                 % TODO: mech team
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% =========================================================================
%  ATTITUDE DETERMINATION & CONTROL (ADCS)
% =========================================================================

%% -------------------------------------------------------------------------
%  9. Star Tracker
%
%  Exterior-mounted, always on (needed for attitude knowledge in safe mode).
%  G = 1 W/K is INTENTIONALLY low: star trackers use vibration-isolated
%  mounts that also limit heat transfer. Higher G risks microvibration.
%  Use max power (8.9 W) for conservative worst case.
%
%  TODO: Confirm continuous vs intermittent operation with ADCS team.
%  TODO: Confirm T_op and T_surv limits from vendor datasheet.
%  TODO: Confirm exterior mounting location and view constraints.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Star Tracker';
comps(i).Q_hot_W      = 8.9;    % max of 5.6–8.9 W range (conservative)
comps(i).Q_cold_W     = 8.9;    % always on — same power
comps(i).T_op_min_C   = -30;    % TODO: vendor datasheet
comps(i).T_op_max_C   = 60;     % TODO: vendor datasheet
comps(i).T_surv_min_C = -30;    % TODO: vendor datasheet
comps(i).T_surv_max_C = 60;     % TODO: vendor datasheet
comps(i).G_W_K        = 1.0;    % Vibration-isolated exterior bracket
                                 % Low G intentional — do not increase without
                                 % checking microvibration requirements
comps(i).always_on    = true;   % Required in safe mode for attitude knowledge
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  10. IMU
%
%  Always on — attitude knowledge needed in all modes.
%  G = 5 W/K (rigid interior mount with TIM) keeps ΔT = 12/5 = 2.4 C.
%
%  TODO: Confirm model and exact specs with ADCS team.
%  TODO: Confirm G from mount design.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'IMU';
comps(i).Q_hot_W      = 12;     % 12 W avg per team PDF (OD5, April 2026)
comps(i).Q_cold_W     = 12;     % always on — same power in safe mode
comps(i).T_op_min_C   = -62;    % Team component layout PDF (OD5, April 2026)
                                  % Northrop Grumman LN-200 class IMU rated to -62°C
comps(i).T_op_max_C   = 85;     % Team component layout PDF
comps(i).T_surv_min_C = -62;    % 5°C margin below op min; typical for military-grade IMUs
comps(i).T_surv_max_C = 85;     % TODO: confirm with vendor
comps(i).G_W_K        = 5.0;    % Rigid interior mount with TIM
                                 % TODO: mech team to verify
comps(i).always_on    = true;
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  11. RCWs (3x)  
%
%  Reaction Control Wheels for attitude control. All four treated as
%  one thermal node. G = 5 W/K per solid structural interface.
%  ΔT = 27.5/5 = 5.5 C — acceptable.
%  TODO: Confirm always_on status with ADCS team (may need to spin up in safe mode).
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'RCWs(3x)';
comps(i).Q_hot_W      = 38;   % 38 W max per team PDF (OD5, April 2026)
comps(i).Q_cold_W     = 0;      % assumed off/standby in safe mode  % TODO: confirm
comps(i).T_op_min_C   = -30;    % Checked
comps(i).T_op_max_C   = 70;     % Checked
comps(i).T_surv_min_C = -30;    % Extended 10°C below op min; confirm with vendor
comps(i).T_surv_max_C = 70;     % Checked
comps(i).G_W_K        = 5.0;    % Solid structural interface with TIM
                                 % TODO: mech team to verify
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% =========================================================================
%  COMMUNICATIONS
% =========================================================================

%% -------------------------------------------------------------------------
%  12. QH Antennas (2x)  (Quadrifilar Helix)
%
%  TT&C — always on in all modes. 10 W draw; only 1.5 W becomes heat
%  (remainder radiated as RF power).
%  G = 0.5 W/K is intentionally low (exterior mount, small contact area).
%  ΔT = 1.5/0.5 = 3 C — negligible.
%
%  TODO: Confirm 1.5 W heat fraction (vs RF output) with TT&C team.
%  TODO: Confirm survival limits from antenna vendor.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'QH Antennas (2x)';
comps(i).Q_hot_W      = 1.5;    % team value — thermal only (not RF power)
comps(i).Q_cold_W     = 1.5;    % always on
comps(i).T_op_min_C   = -40;
comps(i).T_op_max_C   = 80;
comps(i).T_surv_min_C = -40;    % 
comps(i).T_surv_max_C = 80;     % 
comps(i).G_W_K        = 0.5;    % Exterior mount, small contact area
comps(i).always_on    = true;
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  13. GPS Antennas (2x)
%
%  Negligible heat generation (team value). Conservative 0.5 W placeholder.
%  TODO: Confirm heat dissipation — replace with 0 if team confirms negligible.
%  TODO: Confirm continuous vs intermittent operation.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'GPS Antennas (2x)';
comps(i).Q_hot_W      = 0.5;    % conservative placeholder  % TODO: confirm with team
comps(i).Q_cold_W     = 0.5;    % always on
comps(i).T_op_min_C   = -55;
comps(i).T_op_max_C   = 85;
comps(i).T_surv_min_C = -55;    % TODO: antenna vendor
comps(i).T_surv_max_C = 85;     % TODO: antenna vendor
comps(i).G_W_K        = 0.5;    % Exterior patch mount
comps(i).always_on    = true;
comps(i).decoupled    = false;

%% =========================================================================
%  AVIONICS
% =========================================================================

%% -------------------------------------------------------------------------
%  14. Flight Computer  [CONFIRM MODEL WITH AVIONICS TEAM]
%
%  Always on — processes all commands including safe-mode.
%  G = 5 W/K (bolted PCB mount with TIM). ΔT = 5/5 = 1 C — minimal.
%
%  TODO: Confirm model, power, and T_op limits with avionics team.
%  TODO: Confirm G from PCB mounting design.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Flight Computer';
comps(i).Q_hot_W      = 5;      % TODO: placeholder — confirm with avionics team
comps(i).Q_cold_W     = 5;      % always on — same power
comps(i).T_op_min_C   = -20;    % 
comps(i).T_op_max_C   = 70;     % 
comps(i).T_surv_min_C = -20;    %
comps(i).T_surv_max_C = 70;     % 
comps(i).G_W_K        = 5.0;    % Bolted PCB mount with TIM
                                 % TODO: mech team to verify
comps(i).always_on    = true;
comps(i).decoupled    = false;

%% =========================================================================
%  PROPULSION (RCS)
% =========================================================================

%% -------------------------------------------------------------------------
%  15. RCS Thrusters (16x)
%
%  All 16 treated as one node. Fire intermittently; off in safe mode.
%  50 W is a conservative worst-case placeholder. The actual time-averaged
%  structural soakback depends on firing duty cycle — could be much lower
%  (the 4_20 version used 2 W time-averaged, which may be more realistic).
%
%  TODO: Get time-averaged soakback from propulsion team based on duty cycle.
%  TODO: Get temperature limits from propulsion team / thruster datasheet.
%  TODO: Confirm G from titanium mount bracket geometry.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'RCS Thrusters (16x)';
comps(i).Q_hot_W      = 50;     % conservative placeholder  % TODO: propulsion team
                                 % Note: duty-cycle-averaged value could be ~2 W (4_20 estimate)
comps(i).Q_cold_W     = 0;
comps(i).T_op_min_C   = -10;    % Checked
comps(i).T_op_max_C   = 150;    % Checked
comps(i).T_surv_min_C = -10;    % Checked
comps(i).T_surv_max_C = 150;    % Checked
comps(i).G_W_K        = 1.0;    % Isolated titanium bracket mounts
                                 % TODO: mech team to verify
comps(i).always_on    = false;
comps(i).decoupled    = false;

%% -------------------------------------------------------------------------
%  16. Plasma Contactor
%
%  Rear-mounted near thruster. Off in safe mode.
%  TODO: All values are placeholders — obtain from subsystem team.
% -------------------------------------------------------------------------
i = i+1;
comps(i).name         = 'Plasma Contactor';
comps(i).Q_hot_W      = 10;     % TODO: placeholder — subsystem team
comps(i).Q_cold_W     = 0;
comps(i).T_op_min_C   = -50;    % TODO: subsystem team
comps(i).T_op_max_C   = 100;    % TODO: subsystem team
comps(i).T_surv_min_C = -50;    % TODO: subsystem team
comps(i).T_surv_max_C = 100;    % TODO: subsystem team
comps(i).G_W_K        = 1.0;    % Exterior rear mount
                                 % TODO: mech team
comps(i).always_on    = false;
comps(i).decoupled    = false;

end
