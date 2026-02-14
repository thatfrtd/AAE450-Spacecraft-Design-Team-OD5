function [delta_v, time] = Hohmann_Transfer_Altitude_Change(h_0, h_f, Re, mu)

    arguments (Input)
        h_0 = 600; %initial altitude [km]
        h_f = 2000; %final altitude [km]
        Re = 6371;  %mean earth radius [km]
        mu = 398600.435507;  %earth gravitational parameter [km^3/s^2]
    end

    r_initial = Re + h_0; %radius of initial orbit [km]
    r_final = Re + h_f;  %radius of final orbit [km]

    a_transfer = (r_initial + r_final) / 2; %semi-major axis of the elliptical transfer orbit [km]

    v_initial = sqrt(mu / r_initial); %[km/s]

    v_final = sqrt(mu / r_final);  %[km/s]

    v_transfer_periapsis = sqrt(mu * (2 / r_initial - 1 / a_transfer)); %[km/s]

    v_transfer_apoapsis = sqrt(mu * (2 / r_final - 1 / a_transfer)); %[km/s]

    delta_v = abs(v_final - v_transfer_apoapsis) + abs(v_initial - v_transfer_periapsis);  %[km/s]

    time = pi * sqrt(a_transfer ^ 3 / mu);  %time required for the manuver [s]

end