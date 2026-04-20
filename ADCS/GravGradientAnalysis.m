%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Saturation Sims
%
% This is meant to simulate ACS Saturation
% This is V1. It's low fidelity with simple geometry and spin on one axis
% along with using only reaction wheels and CMGs
% Last modified:
%   - 4/9
%   - Atharva
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Orbital Parameters
mu = 3.986004418e14; % Earth gravitational parameter (m^3/s^2)
R_earth = 6371e3; % Earth radius (m)
alt = 800e3; % Orbital altitude (m)  assuming 800km debris orbit
R_orb = R_earth + alt; % Orbital radius (m)

%% Spacecraft MOI Calculation
% Chaser (Solid Cylinder)
m_c = 2000;  r_c = 2.0;  h_c = 10.0;
I_xx_c = (1/12) * m_c * (3*r_c^2 + h_c^2);
I_zz_c = (1/2) * m_c * r_c^2;

% Target Rocket Body 
m_t = 3000;  r_t = 2.0;  h_t = 10.0;
I_xx_t = (1/12) * m_t * (6*r_t^2 + h_t^2);
I_zz_t = m_t * r_t^2;

% Parallel Axis Theorem to find new combined CM
distance_between_CMs = (h_c/2) + (h_t/2); % 10 meters
z_cm_new = (m_c*0 + m_t*distance_between_CMs) / (m_c + m_t); % Relative to chaser CM

d_c = z_cm_new;                        
d_t = distance_between_CMs - z_cm_new; 

I_xx_stack = (I_xx_c + m_c*d_c^2) + (I_xx_t + m_t*d_t^2);
I_yy_stack = I_xx_stack; % Axisymmetric
I_zz_stack = I_zz_c + I_zz_t;

%% 3. Gravity Gradient Disturbance

theta_error_deg = 5.0; 
theta = deg2rad(theta_error_deg);

% Pitch axis 
T_d = (3 * mu / (2 * R_orb^3)) * abs(I_xx_stack - I_zz_stack) * sin(2 * theta);

fprintf('Combined I_xx: %.1f kg*m^2\n', I_xx_stack);
fprintf('Combined I_zz: %.1f kg*m^2\n', I_zz_stack);
fprintf('Calculated Gravity Gradient Torque: %.4f Nm\n\n', T_d);

%% Hardware Parameters & RCS
RW_H_max = 25;% Nms
CMG_H_max = 15; % Nms
tau_RCS = 22 * 1.0; % 22 N thrust at 1 m lever arm

sim_time_hours = 24;    
t_total_sec = sim_time_hours * 3600;
t = linspace(0, t_total_sec, 100000); 

%% Simulation & Saturation Metrics
RW_H = mod(T_d * t, RW_H_max);
CMG_H = mod(T_d * t, CMG_H_max);

RW_dumps = floor((T_d * t_total_sec) / RW_H_max);
CMG_dumps = floor((T_d * t_total_sec) / CMG_H_max);

% Time required to dump a fully saturated system
RW_burn_per_dump = RW_H_max / tau_RCS;
CMG_burn_per_dump = CMG_H_max / tau_RCS;

fprintf('--- Reaction Wheel (Capacity: %d Nms) ---\n', RW_H_max);
fprintf('Saturation Events: %d | Total RCS Burn: %.2f sec\n\n', RW_dumps, RW_dumps * RW_burn_per_dump);

fprintf('--- Control Moment Gyro (Capacity: %d Nms) ---\n', CMG_H_max);
fprintf('Saturation Events: %d | Total RCS Burn: %.2f sec\n', CMG_dumps, CMG_dumps * CMG_burn_per_dump);

%% Plotting
figure('Name', 'True Gravity Gradient Saturation');
subplot(2,1,1); plot(t/3600, RW_H, 'b'); yline(RW_H_max, 'r--'); 
title('Reaction Wheel Accumulation'); ylabel('Nms'); grid on;

subplot(2,1,2); plot(t/3600, CMG_H, 'g'); yline(CMG_H_max, 'r--');
title('Control Moment Gyro Accumulation'); ylabel('Nms'); xlabel('Hours'); grid on;

%% Propellant Mass Calculation
g0 = 9.80665; % (m/s^2)
I_sp = 297; % Specific Impulse in seconds (e.g., Hydrazine monopropellant)

% Calculate fuel mass required for desaturation (kg)
RW_fuel_mass_kg = RW_total_impulse / (I_sp * g0);
CMG_fuel_mass_kg = CMG_total_impulse / (I_sp * g0);

fprintf('--- Estimated Fuel Mass Required (Isp = %d s) ---\n', I_sp);
fprintf('Reaction Wheel:      %.3f kg of propellant\n', RW_fuel_mass_kg);
fprintf('Control Moment Gyro: %.3f kg of propellant\n\n', CMG_fuel_mass_kg);