function [delta_v, time] = Spiral_Transfer_Altitude_Change(h_0, h_f, a_theta, Re, mu)

    arguments (Input)
        h_0 = 600; %initial altitude [km]
        h_f = 2000; %final altitude [km]
        a_theta = 3 * 10 ^ -5 * 10 ^ -3;  %acceleration in the theta (tangential) direction [km / s ^ 2]
        Re = 6371;  %mean earth radius [km]
        mu = 398600.435507;  %earth gravitational parameter [km^3/s^2]
    end
    

    delta_v = abs(sqrt(mu / (h_0 + Re)) - sqrt(mu / (h_f + Re)));  %delta v required to lift the spacecraft between orbits using a spiral transfer [km/s]
    time = delta_v / a_theta;  %time required for the manuver [s]
end