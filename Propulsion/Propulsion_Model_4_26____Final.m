%% AAE 450- Propulsion Tank Sizing System
clc
clear all 
close all
%% Notes
% - This model correctly added extra Chemical ACS for debris avoidance manuevers:

% Parking(Straight off the launch vehicle)- 0.55 m/s
% 
% Transfer to Debris: 0.14 m/s
% Rendezvous= 0.5 m/s
% De-Orbit= 1.08 m/s
% Re-Orbit= 65.94 m/s

% Repeat last four for the next 2 debris(except for the re-orbit for the
% final debris, as the ship will go down with the final debris)


% - Need to add xenon for E passivation(Negligible)

%% Overall Mission Ops
% Stage 1: Electric prop Transfer to Debris
% Stage 2: Chemical Prop ACS + Transfer to Capture Debris
% Stage 3: Electric prop transfer + Chemical prop ACS to De-Orbit Debris

%% Constants
g = 9.81; % m/s^2

%% =========================================================
%% SHIP SELECTOR — change this single value to switch ships
active_ship = 3;  % 1, 2, or 3
%% =========================================================
%% Dry Mass (kg)
sc_dry_mass_no_prop_system= 1020; %kg - Whatever Dry Mass Gives ~1500 kg wet mass(970 kg), estimating about 1020 kg to account for thermal management system

%% Per-Ship Mission Profile Definitions
% --- Ship 1 ---
ship1.sc_dry_mass = sc_dry_mass_no_prop_system; 
ship1.dV_eprop_initial = 600;   
ship1.legs(1).dV_eprop_to_debris  = 0;       
ship1.legs(1).dV_chem_rendezvous  = 20;
ship1.legs(1).dV_eprop_to_dropoff = 194.3;
ship1.legs(1).dV_chem_rcs_dropoff = 15;
ship1.legs(1).dV_chem_rcs_boost   = 15;
ship1.legs(1).debris_mass_leg     = 2000;
ship1.legs(2).dV_eprop_to_debris  = 307.3;
ship1.legs(2).dV_chem_rendezvous  = 20;
ship1.legs(2).dV_eprop_to_dropoff = 196.4;
ship1.legs(2).dV_chem_rcs_dropoff = 15;
ship1.legs(2).dV_chem_rcs_boost   = 15;
ship1.legs(2).debris_mass_leg     = 4000;
ship1.legs(3).dV_eprop_to_debris  = 329.8;
ship1.legs(3).dV_chem_rendezvous  = 20;
ship1.legs(3).dV_eprop_to_dropoff = 178.2;
ship1.legs(3).dV_chem_rcs_dropoff = 15;
ship1.legs(3).dV_chem_rcs_boost   = 0;      
ship1.legs(3).debris_mass_leg     = 2000;

% --- Ship 2 ---
ship2.sc_dry_mass = sc_dry_mass_no_prop_system; 
ship2.dV_eprop_initial = 600;
ship2.legs(1).dV_eprop_to_debris  = 0;
ship2.legs(1).dV_chem_rendezvous  = 20;
ship2.legs(1).dV_eprop_to_dropoff = 231.8;
ship2.legs(1).dV_chem_rcs_dropoff = 15;
ship2.legs(1).dV_chem_rcs_boost   = 15;
ship2.legs(1).debris_mass_leg     = 2000;
ship2.legs(2).dV_eprop_to_debris  = 729.1;
ship2.legs(2).dV_chem_rendezvous  = 20;
ship2.legs(2).dV_eprop_to_dropoff = 177.7;
ship2.legs(2).dV_chem_rcs_dropoff = 15;
ship2.legs(2).dV_chem_rcs_boost   = 15;
ship2.legs(2).debris_mass_leg     = 4000;
ship2.legs(3).dV_eprop_to_debris  = 223.9;
ship2.legs(3).dV_chem_rendezvous  = 20;
ship2.legs(3).dV_eprop_to_dropoff = 157.1;
ship2.legs(3).dV_chem_rcs_dropoff = 15;
ship2.legs(3).dV_chem_rcs_boost   = 0;
ship2.legs(3).debris_mass_leg     = 4000;

% --- Ship 3 ---
ship3.sc_dry_mass = sc_dry_mass_no_prop_system; 
ship3.dV_eprop_initial = 600;
ship3.legs(1).dV_eprop_to_debris  = 0;
ship3.legs(1).dV_chem_rendezvous  = 20;
ship3.legs(1).dV_eprop_to_dropoff = 175.9;
ship3.legs(1).dV_chem_rcs_dropoff = 15;
ship3.legs(1).dV_chem_rcs_boost   = 15;
ship3.legs(1).debris_mass_leg     = 4000;
ship3.legs(2).dV_eprop_to_debris  = 571.8;
ship3.legs(2).dV_chem_rendezvous  = 20;
ship3.legs(2).dV_eprop_to_dropoff = 193.8;
ship3.legs(2).dV_chem_rcs_dropoff = 15;
ship3.legs(2).dV_chem_rcs_boost   = 15;
ship3.legs(2).debris_mass_leg     = 4000;
ship3.legs(3).dV_eprop_to_debris  = 288.3;
ship3.legs(3).dV_chem_rendezvous  = 20;
ship3.legs(3).dV_eprop_to_dropoff = 200.1;
ship3.legs(3).dV_chem_rcs_dropoff = 15;
ship3.legs(3).dV_chem_rcs_boost   = 0;
ship3.legs(3).debris_mass_leg     = 4000;

%% Select Active Ship
ships = {ship1, ship2, ship3};
ship = ships{active_ship};
sc_dry_mass = ship.sc_dry_mass;
n_legs = length(ship.legs);

fprintf('====================================================\n');
fprintf('  ACTIVE SHIP: Ship %d  (%d debris legs)\n', active_ship, n_legs);
fprintf('====================================================\n');

%% Important Inputs/ Data
% MF set to 1.0 — inert propulsion mass is now accounted for explicitly below
Isp_electric = 4155; 
Isp_chemical = 297;  
MF_electric = 1.0;   % inert mass handled explicitly via prop_inert_mass
MF_chemical = 1.0;   % inert mass handled explicitly via prop_inert_mass
rho_xenon = 2045.5; 
rho_NTO = 1440; 
rho_MMH = 880; 
OF_ratio = 1.61; 
ullage_xe = 0.15;
ullage_chem = 0.05;
% Margins from Alayna's SMAD (Table 17-10)
margin_off_nominal_performance= 0.01;
margin_off_operations_performance= 0.01;
margin_mission = 0.075;
margin_contingency= 0.075;
margin_residual= 0.02;
margin_loading_uncertainty = 0.005;
prop_margin = margin_off_nominal_performance + margin_off_operations_performance + margin_mission + margin_contingency + margin_residual + margin_loading_uncertainty;


%% ACS-specific Propellant Sizing
m_prop_plume_margin = 3.5*2*1.15;  
m_prop_ACS_s2 = 7;                
m_prop_ACS_s2_margin = m_prop_plume_margin + m_prop_ACS_s2;  % chemical adder per rendezvous leg
m_prop_ACS_s3_margin = 3;                                    % chemical adder per drop-off leg

%% Mini-Debris Avoidance ACS — delta-V values (m/s, chemical)
% Applied once at parking (post-launch), then per leg:
%   Transfer-to-debris phase, Rendezvous phase, De-Orbit phase, Re-Orbit phase.
% Re-Orbit is skipped on the final leg (ship de-orbits with last debris).
dV_minidebris_parking    = 0.55;   % m/s — once, at parking after launch vehicle separation
dV_minidebris_transfer   = 0.14;   % m/s — per leg, during eProp transfer-to-debris phase
dV_minidebris_rendezvous = 0.5;    % m/s — per leg, during rendezvous phase
dV_minidebris_deorbit    = 1.08;   % m/s — per leg, during de-orbit / dropoff phase
dV_minidebris_reorbit    = 65.94;  % m/s — per leg, EXCEPT final leg (ship stays down)

%% =========================================================
%% EXPLICIT PROPULSION HARDWARE MASS INPUTS
% These are fixed hardware masses independent of propellant load.
% Tank shell mass is computed each iteration and added here automatically.
% eProp system (NEXT-C per thruster string):
m_hw_thruster_eprop = 14;    % kg  NEXT-C thruster + harness (per unit)
m_hw_PPU            = 36;    % kg  Power Processing Unit (per unit)
m_hw_gimbal         = 6;     % kg  Gimbal (per unit)
m_hw_PMS            = 8;     % kg  Propellant Management System (HPA+LPA+plumbing)
n_eprop_thrusters   = 1;     % number of NEXT-C thruster strings
m_hw_eprop_fixed    = n_eprop_thrusters * (m_hw_thruster_eprop + m_hw_PPU + m_hw_gimbal) + m_hw_PMS;

% Chemical system (MOOG DST-12 per thruster):
m_hw_thruster_chem  = 0.64;  % kg  DST-12 per thruster
n_chem_thrusters    = 16;     % number of DST-12 thrusters (adjust as needed)
m_hw_feed_chem      = 30;    % kg  valves, manifolds, pressurant system, structure
m_hw_chem_fixed     = n_chem_thrusters * m_hw_thruster_chem + m_hw_feed_chem;
%% =========================================================

%% =========================================================
%% ITERATIVE SIZING LOOP
% Converges propellant and tank shell mass together, since tank mass
% depends on propellant, which depends on total spacecraft mass including tanks.
% =========================================================
conv_tol   = 0.01;   % kg — convergence tolerance
max_iter   = 20;
prop_inert_mass = 1; % initial guess: arbitrary

for iter = 1:max_iter
    prev_inert = prop_inert_mass;

    % Effective dry mass this iteration
    m_eff_dry = sc_dry_mass + prop_inert_mass;

    %% LEG-BY-LEG PROPELLANT SIZING (BACKWARDS INTEGRATION)
    leg_eprop_raw = zeros(1, n_legs);
    leg_chem_raw  = zeros(1, n_legs);
    m_running_wet = m_eff_dry;

    for k = n_legs:-1:1
        deb = ship.legs(k).debris_mass_leg;

        % ---- Mini-debris ACS: Re-Orbit (skipped on final leg) ----
        % Applied after dropoff boost, before moving to next leg.
        % In backwards integration this is the first burn we undo for leg k
        % (it occurs after boost in forward time, so it is the outermost step here).
        if k < n_legs   % not the final leg
            dV_reorbit_md = dV_minidebris_reorbit;
            MR = exp(dV_reorbit_md / (Isp_chemical * g));
            mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_chemical) / MF_chemical);
            leg_chem_raw(k) = leg_chem_raw(k) + mp;
            m_running_wet   = m_running_wet   + mp;
        end

        % 1. Post-Dropoff Boost (No Debris)
        dV_boost = ship.legs(k).dV_chem_rcs_boost;
        if dV_boost > 0
            MR = exp(dV_boost / (Isp_chemical * g));
            mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_chemical) / MF_chemical);
            leg_chem_raw(k) = leg_chem_raw(k) + mp;
            m_running_wet = m_running_wet + mp;
        end

        % 2. RCS at Dropoff (With Debris) + Mini-debris ACS: De-Orbit
        m_running_wet = m_running_wet + deb;
        % ---- Mini-debris ACS: De-Orbit burn (during dropoff / de-orbit phase) ----
        MR = exp(dV_minidebris_deorbit / (Isp_chemical * g));
        mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_chemical) / MF_chemical);
        leg_chem_raw(k) = leg_chem_raw(k) + mp;
        m_running_wet   = m_running_wet   + mp;
        % Existing RCS dropoff burn
        dV_rcs_drop = ship.legs(k).dV_chem_rcs_dropoff;
        MR = exp(dV_rcs_drop / (Isp_chemical * g));
        mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_chemical) / MF_chemical);
        leg_chem_raw(k) = leg_chem_raw(k) + mp;
        m_running_wet = m_running_wet + mp;

        % 3. eProp to Dropoff (With Debris) + Chemical ACS Adder (m_prop_ACS_s3_margin)
        m_running_wet = m_running_wet + m_prop_ACS_s3_margin;
        leg_chem_raw(k) = leg_chem_raw(k) + m_prop_ACS_s3_margin;  % chemical ACS → chem bucket
        dV_e2drop = ship.legs(k).dV_eprop_to_dropoff;
        MR = exp(dV_e2drop / (Isp_electric * g));
        mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_electric) / MF_electric);
        leg_eprop_raw(k) = leg_eprop_raw(k) + mp;
        m_running_wet = m_running_wet + mp;

        % 4. Rendezvous (No Debris) + Chemical ACS Adder (m_prop_ACS_s2_margin)
        %    + Mini-debris ACS: Rendezvous burn
        m_running_wet = m_running_wet - deb;
        m_running_wet = m_running_wet + m_prop_ACS_s2_margin;
        leg_chem_raw(k) = leg_chem_raw(k) + m_prop_ACS_s2_margin;  % Chem ACS adder → Chem bucket
        % ---- Mini-debris ACS: Rendezvous burn ----
        MR = exp(dV_minidebris_rendezvous / (Isp_chemical * g));
        mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_chemical) / MF_chemical);
        leg_chem_raw(k) = leg_chem_raw(k) + mp;
        m_running_wet   = m_running_wet   + mp;
        % Existing rendezvous burn
        dV_rend = ship.legs(k).dV_chem_rendezvous;
        MR = exp(dV_rend / (Isp_chemical * g));
        mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_chemical) / MF_chemical);
        leg_chem_raw(k) = leg_chem_raw(k) + mp;
        m_running_wet = m_running_wet + mp;

        % 5. eProp to Debris (No Debris) + Mini-debris ACS: Transfer burn
        % ---- Mini-debris ACS: Transfer-to-Debris burn ----
        MR = exp(dV_minidebris_transfer / (Isp_chemical * g));
        mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_chemical) / MF_chemical);
        leg_chem_raw(k) = leg_chem_raw(k) + mp;
        m_running_wet   = m_running_wet   + mp;
        % Existing eProp transfer burn
        dV_e2d = ship.legs(k).dV_eprop_to_debris;
        if dV_e2d > 0
            MR = exp(dV_e2d / (Isp_electric * g));
            mp = m_running_wet * (MR - 1) / (1 - (MR - 1) * (1 - MF_electric) / MF_electric);
            leg_eprop_raw(k) = leg_eprop_raw(k) + mp;
            m_running_wet = m_running_wet + mp;
        end
    end

    % ---- Mini-debris ACS: Parking burn (once, immediately after launch vehicle separation) ----
    MR_parking_md = exp(dV_minidebris_parking / (Isp_chemical * g));
    mp_parking_md = m_running_wet * (MR_parking_md - 1) / (1 - (MR_parking_md - 1) * (1 - MF_chemical) / MF_chemical);
    m_chem_parking_md = mp_parking_md;   % tracked separately for reporting
    m_running_wet = m_running_wet + mp_parking_md;

    % Initial Transfer (eProp)
    MR_init = exp(ship.dV_eprop_initial / (Isp_electric * g));
    m_eprop_initial_raw = m_running_wet * (MR_init - 1) / (1 - (MR_init - 1) * (1 - MF_electric) / MF_electric);

    %% Total raw propellant and Margins
    total_eprop_raw = m_eprop_initial_raw + sum(leg_eprop_raw);
    total_chem_raw  = sum(leg_chem_raw) + m_chem_parking_md;
    total_prop_electric = total_eprop_raw  * (1 + prop_margin);
    total_prop_chemical = total_chem_raw   * (1 + prop_margin);
    total_propellant    = total_prop_electric + total_prop_chemical;

    %% Tank Volume Sizing
    m_xe  = total_prop_electric;
    m_NTO = total_prop_chemical * (OF_ratio / (1 + OF_ratio));
    m_MMH = total_prop_chemical * (1 / (1 + OF_ratio));

    % Material Properties of Tanks
    sigma_y_Ti64 = 1100*10^6;
    rho_Ti64 = 4430;
    rho_CFRP = 1600;
    sigma_CFRP = 1200e6;
    SF_yield_chem = 1.5;
    SF_ultimate_chem = 2.25;   % optional (for burst check)
    SF_burst_copv = 2.25;
    knockdown_CFRP = 0.8;
    sigma_allow_ch = sigma_y_Ti64 / SF_yield_chem;
    t_min_m = 0.5*10^-3;

    % Xenon Tank
    P_xe = 150e5;
    V_xe_prop = m_xe / rho_xenon;
    V_xe_tank = V_xe_prop * (1 + ullage_xe);
    r_xe = ((3 * V_xe_tank) / (4 * pi))^(1/3);
    D_xe = 2 * r_xe;
    t_liner_xe = 0.8e-3;
    t_comp_xe = (SF_burst_copv * P_xe * r_xe) / (2 * sigma_CFRP * knockdown_CFRP);
    m_shell_xe_ideal = (rho_Ti64 * 4 * pi * r_xe^2 * t_liner_xe) + (rho_CFRP * 4 * pi * r_xe^2 * t_comp_xe);

    % Chemical Tanks
    P_chem = 22e5;
    V_NTO_tank = (m_NTO / rho_NTO) * (1 + ullage_chem);
    r_NTO = ((3 * V_NTO_tank) / (4 * pi))^(1/3);
    D_NTO = 2 * r_NTO;
    t_wall_NTO = max(P_chem * r_NTO / (2 * sigma_allow_ch), t_min_m);
    m_shell_NTO_ideal = rho_Ti64 * 4 * pi * r_NTO^2 * t_wall_NTO;
    V_MMH_tank = (m_MMH / rho_MMH) * (1 + ullage_chem);
    r_MMH = ((3 * V_MMH_tank) / (4 * pi))^(1/3);
    D_MMH = 2 * r_MMH;
    t_wall_MMH = max(P_chem * r_MMH / (2 * sigma_allow_ch), t_min_m);
    m_shell_MMH_ideal = rho_Ti64 * 4 * pi * r_MMH^2 * t_wall_MMH;

    % Hardware Factors
    HF_copv = 1.3;
    HF_chem = 1.4;
    m_shell_xe  = m_shell_xe_ideal  * HF_copv;
    m_shell_NTO = m_shell_NTO_ideal * HF_chem;
    m_shell_MMH = m_shell_MMH_ideal * HF_chem;

    % Total propulsion inert mass = fixed hardware + tank shells
    prop_inert_mass = m_hw_eprop_fixed + m_hw_chem_fixed + m_shell_xe + m_shell_NTO + m_shell_MMH;

    % Check convergence
    if abs(prop_inert_mass - prev_inert) < conv_tol
        fprintf('\nConverged in %d iterations. prop_inert_mass = %.2f kg\n', iter, prop_inert_mass);
        break;
    end
end

%% Print Outputs
fprintf('\n--- PER-LEG PROPELLANT BREAKDOWN (raw, before %.0f%% mission margin) ---\n', prop_margin*100);
fprintf('%-6s | %-22s | %-22s | %-22s\n', 'Leg', 'eProp Used raw (kg)', 'Chem Used raw (kg)', 'Debris Mass (kg)');
fprintf('-------+------------------------+------------------------+------------------------\n');
fprintf('%-6s | %-22.2f | %-22s | %-22s\n', 'Init', m_eprop_initial_raw, 'N/A', 'N/A');
fprintf('%-6s | %-22s | %-22.2f | %-22s\n', 'Park', 'N/A', m_chem_parking_md, 'N/A (mini-debris ACS)');
for k = 1:n_legs
    fprintf('%-6d | %-22.2f | %-22.2f | %-22.1f\n', k, leg_eprop_raw(k), leg_chem_raw(k), ship.legs(k).debris_mass_leg);
end

fprintf('\nSubtotal eProp (raw)               : %.2f kg\n', total_eprop_raw);
fprintf('Subtotal Chemical (raw)            : %.2f kg\n', total_chem_raw);
fprintf('\n--- PROPELLANT TOTALS AFTER SINGLE %.0f%% MISSION PROPELLANT MARGIN ---\n', prop_margin*100);
fprintf('Total Electric Propellant          : %.2f kg\n', total_prop_electric);
fprintf('Total Chemical Propellant          : %.2f kg\n', total_prop_chemical);
fprintf('\n--- PROPULSION INERT MASS BREAKDOWN ---\n');
fprintf('eProp Fixed Hardware (thrusters/PPU/gimbal/PMS) : %.2f kg\n', m_hw_eprop_fixed);
fprintf('Chem Fixed Hardware (thrusters/feed/pressurant) : %.2f kg\n', m_hw_chem_fixed);
fprintf('Xe Tank Shell (w/ HF)                           : %.2f kg\n', m_shell_xe);
fprintf('NTO Tank Shell (w/ HF)                          : %.2f kg\n', m_shell_NTO);
fprintf('MMH Tank Shell (w/ HF)                          : %.2f kg\n', m_shell_MMH);
fprintf('Total Propulsion Inert Mass                     : %.2f kg\n', prop_inert_mass);
fprintf('\nSpacecraft Dry Mass (non-prop)     : %.2f kg\n', sc_dry_mass);
fprintf('Total Launch Mass (dry+inert+prop) : %.2f kg\n', sc_dry_mass + prop_inert_mass + total_propellant);

% % Gas Comparison
% cost_xe_per_kg = 5000; cost_kr_per_kg = 500; cost_ar_per_kg = 10;
% rho_krypton = 2413; rho_argon = 1400;
% fprintf('--- EPROP PROPELLANT COMPARISON: Xe / Kr / Ar ---\n');
% m_ep = total_prop_electric;
% fprintf('  %-10s | Mass (kg) | Vol-Liq (L) | Cost (USD)\n', 'Propellant');
% fprintf('  -----------|-----------|-------------|------------\n');
% fprintf('  %-10s | %9.2f | %11.2f | $%11.0f\n', 'Xenon',   m_ep, (m_ep/rho_xenon)*1000,   m_ep*cost_xe_per_kg);
% fprintf('  %-10s | %9.2f | %11.2f | $%11.0f\n', 'Krypton', m_ep, (m_ep/rho_krypton)*1000, m_ep*cost_kr_per_kg);
% fprintf('  %-10s | %9.2f | %11.2f | $%11.0f\n', 'Argon',   m_ep, (m_ep/rho_argon)*1000,   m_ep*cost_ar_per_kg);

% Tank Table
fprintf('\n                             TABLE 1: PROPELLANT TANK INFORMATION FOR SHIP %d\n', active_ship);
fprintf('%-28s | %-12s | %-12s | %-30s | %-15s | %-12s | %-12s | %-10s\n', ...
    'Propellant (w/ margin)', 'Mass (kg)', 'Dens(kg/m3)', 'Vol (L) (w/ actual ullage)', 'Shell-ideal(kg)', 'Shell-HF(kg)', 'Wet (kg)', 'Diam (m)');
fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%-28s | %-12.1f | %-12.0f | %-30.1f | %-15.2f | %-12.2f | %-12.2f | %-10.3f\n', ...
    'Xenon (Xe) [15% ullage]', m_xe, rho_xenon, V_xe_tank*1000, m_shell_xe_ideal, m_shell_xe, m_xe+m_shell_xe, D_xe);
fprintf('%-28s | %-12.1f | %-12.0f | %-30.1f | %-15.2f | %-12.2f | %-12.2f | %-10.3f\n', ...
    'MON-3 (Ox) [5% ullage]', m_NTO, rho_NTO, V_NTO_tank*1000, m_shell_NTO_ideal, m_shell_NTO, m_NTO+m_shell_NTO, D_NTO);
fprintf('%-28s | %-12.1f | %-12.0f | %-30.1f | %-15.2f | %-12.2f | %-12.2f | %-10.3f\n', ...
    'MMH (Fuel) [5% ullage]', m_MMH, rho_MMH, V_MMH_tank*1000, m_shell_MMH_ideal, m_shell_MMH, m_MMH+m_shell_MMH, D_MMH);

%% FLEET TOTAL LAUNCH MASS
lm1 = calc_launch_mass(ship1, Isp_electric, Isp_chemical, MF_electric, MF_chemical, ...
    m_prop_ACS_s2_margin, m_prop_ACS_s3_margin, prop_margin, prop_inert_mass, g, ...
    dV_minidebris_parking, dV_minidebris_transfer, dV_minidebris_rendezvous, ...
    dV_minidebris_deorbit, dV_minidebris_reorbit);
lm2 = calc_launch_mass(ship2, Isp_electric, Isp_chemical, MF_electric, MF_chemical, ...
    m_prop_ACS_s2_margin, m_prop_ACS_s3_margin, prop_margin, prop_inert_mass, g, ...
    dV_minidebris_parking, dV_minidebris_transfer, dV_minidebris_rendezvous, ...
    dV_minidebris_deorbit, dV_minidebris_reorbit);
lm3 = calc_launch_mass(ship3, Isp_electric, Isp_chemical, MF_electric, MF_chemical, ...
    m_prop_ACS_s2_margin, m_prop_ACS_s3_margin, prop_margin, prop_inert_mass, g, ...
    dV_minidebris_parking, dV_minidebris_transfer, dV_minidebris_rendezvous, ...
    dV_minidebris_deorbit, dV_minidebris_reorbit);

fprintf('\n\n--- FLEET TOTAL LAUNCH MASS ---\n');
fprintf('Ship 1 Launch Mass : %.2f kg\n', lm1);
fprintf('Ship 2 Launch Mass : %.2f kg\n', lm2);
fprintf('Ship 3 Launch Mass : %.2f kg\n', lm3);
fprintf('Fleet Total        : %.2f kg\n', lm1 + lm2 + lm3);

fprintf('\nFleet Total With 3 Of The Heaviest Craft: %.2f kg', lm3*3);
fprintf('\nTotal Tank Diameter Combined: %.2f m\n', D_xe + D_NTO + D_MMH);

%% =========================================================
%% PLOTTING: MISSION SIMULATION & DEBRIS FOCUS
%% =========================================================
% Initialize tracking arrays
event_num = 1:43;
xe_rem = zeros(1, 43); 
acs_rem = zeros(1, 43);
m_curr_wet = sc_dry_mass + prop_inert_mass + total_propellant;

% Starting state (Event 1: Initial Launch)
xe_rem(1) = total_prop_electric;
acs_rem(1) = total_prop_chemical;

% --- 1. Pre-Mission Simulation (Events 1-3) ---
% Event 2: ACS Parking
mp = m_curr_wet * (1 - exp(-dV_minidebris_parking / (Isp_chemical * g)));
m_curr_wet = m_curr_wet - mp;
xe_rem(2) = xe_rem(1); acs_rem(2) = acs_rem(1) - mp;

% Event 3: Orbit Transfer (Ion)
mp = m_curr_wet * (1 - exp(-ship.dV_eprop_initial / (Isp_electric * g)));
m_curr_wet = m_curr_wet - mp;
xe_rem(3) = xe_rem(2) - mp; acs_rem(3) = acs_rem(2);

% --- 2. Leg-by-Leg Forward Simulation (Events 4-43) ---
curr_idx = 4;
for k = 1:n_legs
    deb = ship.legs(k).debris_mass_leg;

    % Event: Transfer (Ion)
    mp = m_curr_wet * (1 - exp(-ship.legs(k).dV_eprop_to_debris / (Isp_electric * g)));
    m_curr_wet = m_curr_wet - mp; 
    xe_rem(curr_idx) = xe_rem(curr_idx-1) - mp; acs_rem(curr_idx) = acs_rem(curr_idx-1);
    curr_idx = curr_idx + 1;

    % Event: ACS Transfer
    mp = m_curr_wet * (1 - exp(-dV_minidebris_transfer / (Isp_chemical * g)));
    m_curr_wet = m_curr_wet - mp; 
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1) - mp;
    curr_idx = curr_idx + 1;

    % Event: ACS Rendezvous
    mp = m_curr_wet * (1 - exp(-(dV_minidebris_rendezvous + ship.legs(k).dV_chem_rendezvous) / (Isp_chemical * g)));
    m_curr_wet = m_curr_wet - mp; 
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1) - mp;
    curr_idx = curr_idx + 1;

    % Event: Proximity Ops ACS
    m_curr_wet = m_curr_wet - m_prop_ACS_s2_margin;
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1) - m_prop_ACS_s2_margin;
    curr_idx = curr_idx + 1;

    % Event: Capture Burn / Capture Marker
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1);
    curr_idx = curr_idx + 1; 

    % Event: Capture (Mass Addition)
    m_curr_wet = m_curr_wet + deb;
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1);
    curr_idx = curr_idx + 1;

    % Event: De-Orbit ACS + Ion Drop-off
    m_curr_wet = m_curr_wet - m_prop_ACS_s3_margin;
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1) - m_prop_ACS_s3_margin;
    curr_idx = curr_idx + 1;

    mp = m_curr_wet * (1 - exp(-ship.legs(k).dV_eprop_to_dropoff / (Isp_electric * g)));
    m_curr_wet = m_curr_wet - mp;
    xe_rem(curr_idx) = xe_rem(curr_idx-1) - mp; acs_rem(curr_idx) = acs_rem(curr_idx-1);
    curr_idx = curr_idx + 1;

    % Event: ACS De-Orbit + RCS Drop-off
    mp = m_curr_wet * (1 - exp(-(dV_minidebris_deorbit + ship.legs(k).dV_chem_rcs_dropoff) / (Isp_chemical * g)));
    m_curr_wet = m_curr_wet - mp;
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1) - mp;
    curr_idx = curr_idx + 1;

    % Event: Release (Mass Removal)
    m_curr_wet = m_curr_wet - deb;
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1);
    curr_idx = curr_idx + 1;

    % Event: Re-Boost + Re-Orbit
    dV_reorbit = 0; if k < n_legs; dV_reorbit = dV_minidebris_reorbit; end
    mp = m_curr_wet * (1 - exp(-(ship.legs(k).dV_chem_rcs_boost + dV_reorbit) / (Isp_chemical * g)));
    m_curr_wet = m_curr_wet - mp;
    xe_rem(curr_idx) = xe_rem(curr_idx-1); acs_rem(curr_idx) = acs_rem(curr_idx-1) - mp;
    curr_idx = curr_idx + 1;
end

% Fill remaining indices with final values (residuals)
xe_rem(curr_idx:end) = xe_rem(curr_idx-1);
acs_rem(curr_idx:end) = acs_rem(curr_idx-1);

% --- 3. Refined Plotting ---
plot_range = 4:43;
yLimTop = max([max(xe_rem(plot_range)), max(acs_rem(plot_range))]) * 1.05;

cap_idx = [9, 22, 35]; 
rel_idx = [14, 27, 40];
bounds = [3.5, 16.5, 29.5, 43.5];
colors = [0.85 0.95 0.95; 0.98 0.95 0.85; 0.98 0.88 0.88]; 

figure('Color', 'w', 'Position', [100 100 1200 650]); hold on; grid on;

% Shading & Interior Labels
leg_titles = {'TARGET 1 LEG', 'TARGET 2 LEG', 'TARGET 3 LEG'};
for i = 1:3
    fill([bounds(i) bounds(i+1) bounds(i+1) bounds(i)], [0 0 yLimTop yLimTop], colors(i,:), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    text(cap_idx(i), yLimTop * 0.92, leg_titles{i}, ...
        'FontWeight', 'bold', 'Horiz', 'center', 'FontSize', 10, ...
        'BackgroundColor', 'w', 'EdgeColor', [0.7 0.7 0.7]);
end

% Plot Propellant Lines
h_xe = plot(event_num(plot_range), xe_rem(plot_range),'LineWidth', 2.5, 'Color', [0.1 0.4 0.6], 'DisplayName', 'Xenon Fuel (kg)');
h_acs = plot(event_num(plot_range), acs_rem(plot_range), 'LineWidth', 2, 'Color', [0.8 0.2 0.2], 'DisplayName', 'ACS Fuel (kg)');

% Event Markers
plot(cap_idx, xe_rem(cap_idx), 'og', 'MarkerFaceColor', 'g', 'MarkerSize', 10, 'DisplayName', 'Target Capture Event');
plot(rel_idx, xe_rem(rel_idx), 'dr', 'MarkerFaceColor', 'r', 'MarkerSize', 10, 'DisplayName', 'Target Release Event');

% Total Prop Mass at Beginning of Mission
text(4.0, xe_rem(1)*.93, sprintf('%.1f kg Xe loaded  ', xe_rem(1)),'Color',[0.1 0.4 0.6],'FontWeight','bold','FontSize',13,'HorizontalAlignment','left');
text(4.0, acs_rem(1)*1.02, sprintf('%.1f kg ACS loaded  ', acs_rem(1)),'Color',[0.8 0.2 0.2],'FontWeight','bold','FontSize',13,'HorizontalAlignment','left');

% Residual Labels at End of Mission
text(43.7, xe_rem(43), sprintf('  %.1f kg Extra Xe', xe_rem(43)), 'Color', [0.1 0.4 0.6], 'FontWeight', 'bold', 'FontSize', 13);
text(43.7, acs_rem(43), sprintf('  %.1f kg Extra ACS', acs_rem(43)), 'Color', [0.8 0.2 0.2], 'FontWeight', 'bold', 'FontSize', 13);

% X-Axis: Descriptive Stage Labels (Horizontal)
% We select key events to label to keep it readable
custom_ticks = [4, 9, 14, 18, 22, 27, 31, 35, 40, 43];
custom_labels = {'Transfer', 'Capture', 'Release', ...
                 'Transfer', 'Capture', 'Release', ...
                 'Transfer', 'Capture', 'Release', 'EOL'};
set(gca, 'XTick', custom_ticks, 'XTickLabel', custom_labels, ...
         'XTickLabelRotation', 0, 'FontWeight', 'bold', 'FontSize', 11);

% Formatting
xlim([3.5 50]); % Extended X-limit to fit residual text
ylim([0 yLimTop]);
xlabel('Mission Execution Stage', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Propellant Mass Remaining (kg)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Propellant Mass Mission Profile (w/ ~20%% total margin)', active_ship), 'FontSize', 14);
legend('Location', 'northeast', 'Box', 'on', 'FontSize', 13);

%% Helper Function
% Note: helper uses the converged prop_inert_mass from the active ship's iteration.
% For a fully rigorous fleet total each ship would converge independently,
% but since all three ships have the same hardware config this is consistent.
function launch_mass = calc_launch_mass(s, Isp_e, Isp_c, MF_e, MF_c, ACS_s2, ACS_s3, margin, prop_inert, g, ...
    dV_md_park, dV_md_xfer, dV_md_rndz, dV_md_deorb, dV_md_reorb)

    m_run = s.sc_dry_mass + prop_inert;
    ep_raw = 0;
    ch_raw = 0;
    n = length(s.legs);

    for k = n:-1:1
        deb = s.legs(k).debris_mass_leg;

        % Mini-debris ACS: Re-Orbit (skip final leg)
        if k < n
            mp = m_run * (exp(dV_md_reorb/(Isp_c*g))-1) / (1-(exp(dV_md_reorb/(Isp_c*g))-1)*(1-MF_c)/MF_c);
            ch_raw = ch_raw + mp; m_run = m_run + mp;
        end

        % Post-dropoff boost
        dV_boost = s.legs(k).dV_chem_rcs_boost;
        if dV_boost > 0
            mp = m_run * (exp(dV_boost/(Isp_c*g))-1) / (1-(exp(dV_boost/(Isp_c*g))-1)*(1-MF_c)/MF_c);
            ch_raw = ch_raw + mp; m_run = m_run + mp;
        end

        % RCS at dropoff + Mini-debris ACS: De-Orbit
        m_run = m_run + deb;
        mp = m_run * (exp(dV_md_deorb/(Isp_c*g))-1) / (1-(exp(dV_md_deorb/(Isp_c*g))-1)*(1-MF_c)/MF_c);
        ch_raw = ch_raw + mp; m_run = m_run + mp;
        mp = m_run * (exp(s.legs(k).dV_chem_rcs_dropoff/(Isp_c*g))-1) / (1-(exp(s.legs(k).dV_chem_rcs_dropoff/(Isp_c*g))-1)*(1-MF_c)/MF_c);
        ch_raw = ch_raw + mp; m_run = m_run + mp;

        % eProp to dropoff + ACS_s3 adder
        m_run = m_run + ACS_s3;
        ch_raw = ch_raw + ACS_s3;
        mp = m_run * (exp(s.legs(k).dV_eprop_to_dropoff/(Isp_e*g))-1) / (1-(exp(s.legs(k).dV_eprop_to_dropoff/(Isp_e*g))-1)*(1-MF_e)/MF_e);
        ep_raw = ep_raw + mp; m_run = m_run + mp;

        % Rendezvous + ACS_s2 adder + Mini-debris ACS: Rendezvous
        m_run = m_run - deb;
        m_run = m_run + ACS_s2;
        ch_raw = ch_raw + ACS_s2;
        mp = m_run * (exp(dV_md_rndz/(Isp_c*g))-1) / (1-(exp(dV_md_rndz/(Isp_c*g))-1)*(1-MF_c)/MF_c);
        ch_raw = ch_raw + mp; m_run = m_run + mp;
        mp = m_run * (exp(s.legs(k).dV_chem_rendezvous/(Isp_c*g))-1) / (1-(exp(s.legs(k).dV_chem_rendezvous/(Isp_c*g))-1)*(1-MF_c)/MF_c);
        ch_raw = ch_raw + mp; m_run = m_run + mp;

        % eProp to debris + Mini-debris ACS: Transfer
        mp = m_run * (exp(dV_md_xfer/(Isp_c*g))-1) / (1-(exp(dV_md_xfer/(Isp_c*g))-1)*(1-MF_c)/MF_c);
        ch_raw = ch_raw + mp; m_run = m_run + mp;
        if s.legs(k).dV_eprop_to_debris > 0
            mp = m_run * (exp(s.legs(k).dV_eprop_to_debris/(Isp_e*g))-1) / (1-(exp(s.legs(k).dV_eprop_to_debris/(Isp_e*g))-1)*(1-MF_e)/MF_e);
            ep_raw = ep_raw + mp; m_run = m_run + mp;
        end
    end

    % Mini-debris ACS: Parking burn (once)
    mp = m_run * (exp(dV_md_park/(Isp_c*g))-1) / (1-(exp(dV_md_park/(Isp_c*g))-1)*(1-MF_c)/MF_c);
    ch_raw = ch_raw + mp; m_run = m_run + mp;

    % Initial eProp transfer
    mp_init = m_run * (exp(s.dV_eprop_initial/(Isp_e*g))-1) / (1-(exp(s.dV_eprop_initial/(Isp_e*g))-1)*(1-MF_e)/MF_e);

    launch_mass = s.sc_dry_mass + prop_inert + (ep_raw + mp_init + ch_raw) * (1 + margin);
end


