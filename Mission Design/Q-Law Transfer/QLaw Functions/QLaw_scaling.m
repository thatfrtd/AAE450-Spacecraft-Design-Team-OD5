function S_oe = QLaw_scaling(a, a_t, m, n, r)
arguments
    a % Current semimajor axis
    a_t % Desired semimajor axis
    m = 3 % Parameter (found by global optimization)
    n = 4 % Parameter (found by global optimization)
    r = 2 % Parameter (found by global optimization)
end

S_a = (1 + (abs(a - a_t) / (m * a_t)) ^ n) ^ (1 / r);
S_oe = [S_a; ones([4, 1])];

end