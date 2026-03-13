function [char_star] = load_charecteristic_values_Earth_meters()
    mu_E = 398600*(1e3)^3; % [m3 / s2] Earth gravitational parameter
    R_E = 6378.1e3; % [m] Earth radius

    char_star = load_charecteristic_values(mu_E, R_E);
    char_star.m = 1000; % [kg]
    char_star.F = char_star.m * char_star.a; % [N]
end