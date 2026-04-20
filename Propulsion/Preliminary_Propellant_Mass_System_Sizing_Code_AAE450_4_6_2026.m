%% AAE 450- Propulsion Tank Sizing System
clc
clear all 
close all
%% Overall Mission Ops
% Stage 1: Electric prop Tansfer to Debris
% Stage 2: Chemical Prop ACS + Transfer to Capture Debris
% Stage 3: Electric prop transfer + Chemical prop ACS to De-Orbit Debris
%% Constants
g = 9.81; % m/s^2
%% Inputs
sc_dry_mass = 2300;  % kg
debris_mass = 4000;  % kg
payload_mass = sc_dry_mass + debris_mass;  
delta_V(1) = 50; 
delta_V(2) = 20;   
delta_V(3) = 600;  

%% Estimated Isp(s) for each prop type
Isp_electric = 4100; % Ion Thruster
Isp_chemical = 297;  % Hypergolic BiProp MMH/MON-3
MF_electric = 0.90;   
MF_chemical = 0.85;   

%% Estimated propellant densities
rho_xenon = 1900; % kg/m^3
rho_NTO = 1440; % kg/m^3
rho_MMH = 880; % kg/m^3
OF_ratio = 1.65; 
ullage = 0.10; 
prop_margin = 0.15; 

%% ACS-specific Propellant Sizing
m_prop_plume_margin = 3.5*2*1.15; 
m_prop_ACS_s2 = 70; 
m_prop_ACS_s2_margin = m_prop_plume_margin + m_prop_ACS_s2; 
m_prop_ACS_s3_margin = 30; 

%% Tsiolkovsky Calculation Loop
Isp_stage = [Isp_electric, Isp_chemical, Isp_electric];
n_stages  = 3;
m_propellant = zeros(1, n_stages);
m_propellant_margin = zeros(1, n_stages);
m_inert = zeros(1, n_stages);
m_wet = zeros(1, n_stages);
m_payload_effective = zeros(1, n_stages);
m_payload_effective(3) = payload_mass;

for i = n_stages:-1:1
    m_pl = m_payload_effective(i);
    MR = exp(delta_V(i) / (Isp_stage(i) * g)); 
    if i == 2
        extra_ACS_mass = m_prop_ACS_s2_margin;
        MF = MF_chemical;
    elseif i == 3
        extra_ACS_mass = m_prop_ACS_s3_margin;
        MF = MF_electric;
    else
        extra_ACS_mass = 0;
        MF = MF_electric;
    end
   
    m_propellant(i) = (m_pl + extra_ACS_mass) * (MR - 1) / (1 - (MR - 1) * (1 - MF) / MF);
    m_propellant_margin(i) = m_propellant(i) * (1 + prop_margin);
    m_inert(i) = m_propellant_margin(i) * (1 - MF) / MF;  
    m_wet(i) = m_pl + m_inert(i) + extra_ACS_mass + m_propellant_margin(i);

    if i > 1
        if i == 3
            m_payload_effective(i-1) = m_wet(i) - debris_mass; 
        else
            m_payload_effective(i-1) = m_wet(i);
        end
    end
end

%% Total propellant per stage
m_prop_stage2_total    = m_propellant_margin(2) + m_prop_ACS_s2_margin;
m_prop_stage3_electric = m_propellant_margin(3);
m_prop_stage3_chemical = m_prop_ACS_s3_margin;

m_final(1) = m_wet(1) - m_propellant_margin(1);
m_final(2) = m_wet(2) - m_prop_stage2_total;
m_final(3) = m_wet(3) - m_prop_stage3_electric - m_prop_stage3_chemical;

total_prop_electric = m_propellant_margin(1) + m_prop_stage3_electric;
total_prop_chemical = m_prop_stage2_total + m_prop_stage3_chemical;
total_propellant    = total_prop_electric + total_prop_chemical;

%% PROPELLANT TANK VOLUME SIZING
m_xe = total_prop_electric; % <--- FIXED: Defining m_xe here
m_NTO = total_prop_chemical * (OF_ratio / (1 + OF_ratio)); 
m_MMH = total_prop_chemical * (1 / (1 + OF_ratio)); 

%% Tankage Material Properties
sigma_y_Ti64 = 1100*10^6;   
rho_Ti64 = 4430;     
rho_CFRP = 1600;     
sigma_CFRP = 1200e6;   
SF_yield_chem = 2.0;      
SF_burst_copv = 1.5;      
sigma_allow_ch = sigma_y_Ti64 / SF_yield_chem;   
t_min_m = 0.5*10^-3;   
% hardware_factor_chem = 1.8;
% hardware_factor_copv = 1.5;

%% E-PROP XENON COPV TANK
P_xe = 150e5;         
V_xe_prop = m_xe / rho_xenon;
V_xe_tank = V_xe_prop * (1 + ullage);    
r_xe = ((3 * V_xe_tank) / (4 * pi))^(1/3);
D_xe = 2 * r_xe;
t_liner_xe = 0.8e-3;                         
t_comp_xe = (SF_burst_copv * P_xe * r_xe) / (2 * sigma_CFRP);
m_shell_xe_ideal = (rho_Ti64 * 4 * pi * r_xe^2 * t_liner_xe) + (rho_CFRP * 4 * pi * r_xe^2 * t_comp_xe);

% ACTUAL Mass (Heritage Scaling)
% V_heritage_xe = 50e-3;
% m_heritage_xe = 7.0;
% m_tank_xe_actual = m_heritage_xe * (V_xe_tank / V_heritage_xe);
%% CHEMICAL BIPROP TITANIUM TANKS
P_chem = 22e5;  
V_NTO_tank = (m_NTO / rho_NTO) * (1 + ullage);
r_NTO = ((3 * V_NTO_tank) / (4 * pi))^(1/3);
D_NTO = 2 * r_NTO;
t_wall_NTO = max(P_chem * r_NTO / (2 * sigma_allow_ch), t_min_m);
m_shell_NTO_ideal = rho_Ti64 * 4 * pi * r_NTO^2 * t_wall_NTO;
% m_tank_NTO_actual = m_shell_NTO_ideal * hardware_factor_chem;

V_MMH_tank = (m_MMH / rho_MMH) * (1 + ullage);
r_MMH = ((3 * V_MMH_tank) / (4 * pi))^(1/3);
D_MMH = 2 * r_MMH;
t_wall_MMH = max(P_chem * r_MMH / (2 * sigma_allow_ch), t_min_m);
m_shell_MMH_ideal = rho_Ti64 * 4 * pi * r_MMH^2 * t_wall_MMH;
% m_tank_MMH_actual = m_shell_MMH_ideal * hardware_factor_chem;

%% Print Stage Sizing Results
stage_names = {'Stage 1: Electric Ion Transfer', 'Stage 2: Chemical Hypergol BiProp Capture', 'Stage 3: Electric Ion De-orbit'};
for i = 1:n_stages
    fprintf('%s\n', stage_names{i});
    fprintf(' Effective Payload Mass : %.2f kg\n', m_payload_effective(i));
    if i == 1
        fprintf(' Electric Propellant (w/ margin): %.2f kg\n', m_propellant_margin(i));
        fprintf(' Electric Prop Inert Mass : %.2f kg\n', m_inert(i));
    elseif i == 2
        fprintf(' Chemical Translational Prop : %.2f kg\n', m_propellant_margin(i));
        fprintf(' Chemical ACS Prop : %.2f kg\n', m_prop_ACS_s2_margin);
        fprintf(' Total Stage 2 Propellant : %.2f kg\n', m_prop_stage2_total);
    elseif i == 3
        fprintf(' Electric Translational Prop : %.2f kg\n', m_prop_stage3_electric);
        fprintf(' Chemical ACS Prop : %.2f kg\n', m_prop_stage3_chemical);
        fprintf(' Total Stage 3 Propellant : %.2f kg\n', ...
        m_prop_stage3_electric + m_prop_stage3_chemical);
end
    fprintf(' Final Mass (after-burn) : %.2f kg\n\n', m_final(i));
end
fprintf('Total Electric Propellant : %.2f kg\n', total_prop_electric);
fprintf('Total Chemical Propellant : %.2f kg\n', total_prop_chemical);
fprintf('Total Propellant Mass : %.2f kg\n', total_propellant);
fprintf('Total Launch Mass (Stage 1 wet) : %.2f kg\n\n', m_wet(1));

%% FINAL SYSTEM OUTPUT TABLE
fprintf('                             TABLE 1: IDEAL TANK SHELL MODEL (EXCLUDING HARDWARE FACTORS)\n');
fprintf('%-28s | %-12s | %-12s | %-24s | %-12s | %-12s | %-10s\n', ...
    'Propellant (w/ 15% margin)', 'Mass (kg)', 'Dens(kg/m3)', 'Vol (L) (w/ 10%% ullage)', 'Shell (kg)', 'Wet (kg)', 'Diam (m)');
fprintf('-------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%-28s | %-12.1f | %-12.0f | %-24.1f | %-12.2f | %-12.2f | %-10.3f\n', 'Xenon (Xe)', m_xe, rho_xenon, V_xe_tank*1000, m_shell_xe_ideal, m_xe+m_shell_xe_ideal, D_xe);
fprintf('%-28s | %-12.1f | %-12.0f | %-24.1f | %-12.2f | %-12.2f | %-10.3f\n', 'MON-3 (Ox)', m_NTO, rho_NTO, V_NTO_tank*1000, m_shell_NTO_ideal, m_NTO+m_shell_NTO_ideal, D_NTO);
fprintf('%-28s | %-12.1f | %-12.0f | %-24.1f | %-12.2f | %-12.2f | %-10.3f\n', 'MMH (Fuel)', m_MMH, rho_MMH, V_MMH_tank*1000, m_shell_MMH_ideal, m_MMH+m_shell_MMH_ideal, D_MMH);
