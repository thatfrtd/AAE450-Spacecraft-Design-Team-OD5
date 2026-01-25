function [char_star] = load_charecteristic_values()
    char_star.mu = 139348062043.343; %[ km3 / s2]
    char_star.l = 149597870.691; %[km]  
    char_star.t = sqrt(char_star.mu^-1 * char_star.l^3); %[s]
    char_star.v = char_star.l/char_star.t; %[km/s]
    char_star.a = char_star.v/char_star.t;%[km/s^2]
end