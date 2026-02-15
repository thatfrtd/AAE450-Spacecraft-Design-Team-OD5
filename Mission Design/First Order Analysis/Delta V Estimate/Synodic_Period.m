function P = Synodic_Period(h_1, h_2, Re, mu)  %Synodic Period of objects 1 and 2 [s]

    arguments (Input)
        h_1 = 600; %average altitude of object 1 [km]
        h_2 = 660; %average altitude of object 2 [km]
        Re = 6371;  %mean earth radius [km]
        mu = 398600.435507;  %earth gravitational parameter [km^3/s^2]
    end
    r_1 = h_1 + Re;

    r_2 = h_2 + Re;
    
    P_1 = 2 * pi * sqrt(r_1 ^ 3 / mu);  %Object 1 Period [s]

    P_2 = 2 * pi * sqrt(r_2 ^ 3 / mu);  %Object 2 Period [s]

    P = (abs(1 / P_1 - 1 / P_2)) ^ -1;  %Synodic Period [s]
end