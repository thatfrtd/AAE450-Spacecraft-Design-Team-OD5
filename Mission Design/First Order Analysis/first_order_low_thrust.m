%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% First Order Low Thrust Calculations
% Author: Saif Jalal, Travis Hastreiter 
% Created On: 13 February, 2026
% Description: Simple calculation of mission delta V using low thrust
% Last Modified On: 14 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Re = 6371;  %mean earth radius [km]
mu = 398600.435507;  %earth gravitational parameter [km^3/s^2]

F = 1;  %Electric Propulsion System Thrust [N] 
M_total = 1000;  %system mass [kg]
a_theta = F / M_total * 10 ^ -3;  %acceleration in the theta (tangential) direction [km / s ^ 2]

event_labels = string();

h_0 = 600; %initial parking orbit altitude [km]
h_f = 660; %debris orbit altitude [km]

event_labels(1) = "Transfer to Debris";
[V_1, t_1] = Spiral_Transfer_Altitude_Change(h_0, h_f, a_theta, Re, mu);  %delta v required to lift the spacecraft between parking and debris orbits using a spiral transfer [km/s]
P = Synodic_Period(h_0, h_f, Re, mu);  %syndic period of spacecraft and debris [s]

event_labels(2) = "Rendezvous and Docking";
[V_2, t_2] = Rendezvouz_and_Docking();  %total delta V used for rendezvous and docking

h_0 = 660; %debris orbit altitude [km]
h_f = 200; %Final debris disposal altitude [km]

event_labels(3) = "Deorbiting";
[V_3, t_3] = Spiral_Transfer_Altitude_Change(h_0, h_f, a_theta, Re, mu);  %delta v required to lower the spacecraft from debris orbit to disposal orbit using a spiral transfer [km/s]

h_drift = 2000; %Altitude to wait for RAAN drift

event_labels(4) = "Raise to Drift";
[V_4, t_4] = Spiral_Transfer_Altitude_Change(h_f, h_drift, a_theta, Re, mu);  %delta v required to lower the spacecraft from debris orbit to disposal orbit using a spiral transfer [km/s]

RAAN_diff = 20; % [deg]
event_labels(5) = sprintf("J2 Drift %.2g deg", RAAN_diff);
RAAN_drift = rad2deg(J2_RAAN_drift(Re + h_drift, 0, 71, mu, Re));
RAAN_drift_target = rad2deg(J2_RAAN_drift(Re + h_0, 0, 71, mu, Re));

RAAN_drift_wait = RAAN_diff / abs(RAAN_drift-RAAN_drift_target);
V_drift = 0;

event_labels(6) = "Lower from Drift";
[V_5, t_5] = Spiral_Transfer_Altitude_Change(h_drift, h_0, a_theta, Re, mu);  %delta v required to lower the spacecraft from debris orbit to disposal orbit using a spiral transfer [km/s]

Delta_V_total = V_1 + V_2 + V_3 + V_4 + V_5 %[km/s]  Total delta V for all the phases [km/s]

Mission_Time = (t_1 + t_2 + t_3 + t_4 + RAAN_drift_wait + t_5) / 3600 / 24  %mission time [days]
%%
figure
bar(event_labels, [V_1, V_2, V_3, V_4, V_drift, V_5])
title("Low Thrust Mission Segment Delta V")
subtitle(sprintf("Total dV: %.3f [km / s]", Delta_V_total))
ylabel("Delta V [km / s]")
grid on

figure 
bar(event_labels, [t_1, t_2, t_3, t_4, RAAN_drift_wait, t_5] / 3600 / 24)
title("Low Thrust Mission Segment Time")
subtitle(sprintf("Total time: %.3f [days]", Mission_Time))
ylabel("Time [days]")
grid on