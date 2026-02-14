clear
clc
close all
Re = 6371;  %mean earth radius [km]
mu = 398600.435507;  %earth gravitational parameter [km^3/s^2]

F = 0.1;  %Electric Propulsion System Thrust [N] 
M_total = 1000;  %system mass [kg]
a_theta = F / M_total * 10 ^ -3;  %acceleration in the theta (tangential) direction [km / s ^ 2]

h_0 = 600; %initial parking orbit altitude [km]
h_f = 660; %debris orbit altitude [km]

[V_1, t_1] = Spiral_Transfer_Altitude_Change(h_0, h_f, a_theta, Re, mu);  %delta v required to lift the spacecraft between parking and debris orbits using a spiral transfer [km/s]
P = Synodic_Period(h_0, h_f, Re, mu);  %syndic period of spacecraft and debris [s]

[V_2, t_2] = Rendezvouz_and_Docking();  %total delta V used for rendezvous and docking

h_0 = 660; %debris orbit altitude [km]
h_f = 200; %Final debris disposal altitude [km]

[V_3, t_3] = Spiral_Transfer_Altitude_Change(h_0, h_f, a_theta, Re, mu);  %delta v required to lower the spacecraft from debris orbit to disposal orbit using a spiral transfer [km/s]

Delta_V_total = V_1 + V_2 + V_3  %[km/s]  Total delta V for all the phases [km/s]

Mission_Time = (t_1 + t_2 + t_3) / 3600 / 24  %mission time [days]