function [char_star] = load_charecteristic_values_Earth()
    mu_E = 398600; % [km3 / s2] Earth gravitational parameter
    R_E = 6378.1; % [km] Earth radius

    char_star = load_charecteristic_values(mu_E, R_E);
    char_star.m = 1000; % [kg]
    char_star.F = char_star.m * char_star.a; % [kN]
end