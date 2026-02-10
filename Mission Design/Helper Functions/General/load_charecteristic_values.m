function [char_star] = load_charecteristic_values(mu_star, l_star)
    char_star.mu = mu_star; % [LU3 / TU2]
    char_star.l = l_star; % [LU]  
    char_star.t = sqrt(char_star.mu^-1 * char_star.l^3); % [TU]
    char_star.v = char_star.l/char_star.t; % [LU / TU]
    char_star.a = char_star.v/char_star.t; % [LU / TU^2]
end