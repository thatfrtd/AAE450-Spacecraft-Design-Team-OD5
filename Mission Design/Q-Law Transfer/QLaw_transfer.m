function [transfer] = QLaw_transfer(x_keplerian_d, x_keplerian_c, mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, options)
%QLAW_TRANSFER Orbit transfer using Q-Law
%   Modified Equinoctial based Q-Law transfer. Uses Keplerian elements to
%   define spacecraft and target orbits. Q-Law uses a heuristic to
%   determine when to thrust and when to coast. There are >=10 parameters
%   that must be tuned to get the best results.
arguments
    x_keplerian_d % Spacecraft orbit
    x_keplerian_c % Target orbit
    mu % Gravitational parameter
    spacecraft_params % S/c Isp, max thrust (F_max), initial mass (m_0), dry mass (m_dry)
    Q_params % W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot (reference direction)
    penalty_params % Periapsis penalty function parmaeters: r_p_min, k, W_p
    Qdot_opt_params % Parameters for the optimization needed to determine efficiencies
    options.integration_tolerance = 1e-6 % Tolerance for integration of orbit
    options.angular_step = deg2rad(20) % Advised to be in range [5, 20] deg
    options.R_c = 1 % Quantity controlling penalization of remaining distance from target
    options.iter_max = 1e4 % Maximum number of timesteps
    options.a_disturbance = @(t,x) [0;0;0] % @(t, x) [3, 1] Disturbing accelerations
    options.return_dt_dm_only = false % Only return info to evaluate Q_params for global optimization of Qlaw parameters
end

% Constants
g_0 = 9.81e-3; % [km / s2]

% Adjust RAAN based on Theta_rot parameter which defines the reference
% direction (rotates x axis along z axis). Doesn't change actual orbit but
% changes the values of the modified equinoctial elements
x_keplerian_d(4) = x_keplerian_d(4) + Q_params.Theta_rot;
x_keplerian_c(4) = x_keplerian_c(4) + Q_params.Theta_rot;

% NONDIMENSIONALIZE (mu = 1)
char_star = load_charecteristic_values(mu, x_keplerian_c(1));
char_star.m = spacecraft_params.m_0;
char_star.F = char_star.m * char_star.a;
nd_scalar = [char_star.l; ones([5, 1]); char_star.m];
x_keplerian_d_nd = x_keplerian_d ./ nd_scalar(1:6);
x_keplerian_c_nd = x_keplerian_c ./ nd_scalar(1:6);
F_max_nd = spacecraft_params.F_max / 1000 / char_star.F; % F_max in N, char_star.F in kN
r_p_min_nd = penalty_params.r_p_min / char_star.l;
m_0_nd = spacecraft_params.m_0 / char_star.m;
m_dry_nd = spacecraft_params.m_dry / char_star.m;

% Get Modified Equinoctial elements
x_me_d_nd = keplerian_to_modified_equinoctial(x_keplerian_d_nd, []);
x_me_c_nd = keplerian_to_modified_equinoctial(x_keplerian_c_nd, []);

% Gather target elements (a, f, g, h, k)
oe_t = [x_keplerian_c(1); x_me_c_nd(2:5)];

% Get stopping condition
Q_stop = QLaw_stopping_condition(options.R_c, Q_params.W_oe);

% Prepare integration
tolerances = odeset(RelTol=options.integration_tolerance, AbsTol=options.integration_tolerance, Stats = "off");
x_me_mass_nd = [x_me_d_nd; m_0_nd]; % [p, f, g, h, k, L, mass]
u_nd = zeros([3, 1]); % [F_T, F_R, F_N]
alpha = 0;
beta = 0;
t_nd = 0;

% Simulate
iter = 1;
Q = 0;
while iter == 1 || Q(iter - 1) >= Q_stop && iter < options.iter_max && x_me_mass_nd(7, end) > m_dry_nd
    % Create slow QLaw orbital elements (a, f, g, h, k)
    e = sqrt(x_me_mass_nd(2, end) ^ 2 + x_me_mass_nd(3, end) ^ 2);
    a = x_me_mass_nd(1, end) / (1 - e ^ 2); % a = p / (1 - e ^ 2)
    w = 1 + x_me_mass_nd(2, end) .* cos(x_me_mass_nd(6, end)) + x_me_mass_nd(3, end) .* sin(x_me_mass_nd(6, end));
    r = x_me_mass_nd(1, end) / w; % r = p / (1 + f .* cos(L) + g .* sin(L))
    oe = [a; x_me_mass_nd(2:5, end)];
    L = x_me_mass_nd(6, end);

    % Calculate periapsis constraint
    r_p = a * (1 + e);
    P(iter) = periapsis_penalty(r_p, r_p_min_nd, penalty_params.k);
    
    % Calculate Q function
    S_oe = QLaw_scaling(oe(1), oe_t(1), Q_params.m, Q_params.n, Q_params.r);
    oedot_xx = QLaw_oe_max_rate(oe, 1, spacecraft_params.F_max);
    Q(iter) = Q_function(oe, oe_t, penalty_params.W_p, P(iter), S_oe, Q_params.W_oe, oedot_xx);
    
    % Calculate control
    [D, Qdot, partial_Q_partial_oe] = D_Qdot_partial_Q_partial_oe_func(oe, L, oe_t, Q_params.W_oe, Q_params.m, Q_params.n, Q_params.r, F_max_nd, penalty_params.W_p, r_p_min_nd, penalty_params.k);
    [u_nd(:, iter), alpha(iter), beta(iter)] = QLaw_thrust_mapping(D, F_max_nd);

    % Determine if thrusting should happen based on heuristic
    [Qdot_min, Qdot_max] = Qdot_extremize(oe, partial_Q_partial_oe, Qdot_opt_params);
    [eta_a, eta_r] = QLaw_efficiencies(Qdot, Qdot_min, Qdot_max);
    not_coast = (eta_a >= Q_params.eta_a_min ...
              && eta_r >= Q_params.eta_r_min);
    u_nd(:, iter) = u_nd(:, iter) * not_coast;
    a_control = @(x) u_nd(:, iter) / x_me_mass_nd(7, end); % a = F / m
    mdot = -F_max_nd / (spacecraft_params.Isp * g_0 / char_star.v) * not_coast;
    
    % Propagate orbit
    dt_step = options.angular_step * sqrt(r ^ 3 / 1) * sqrt(w) / (1 + e);
    tolerances.InitialStep = dt_step;
    tolerances.MaxStep = dt_step;

    [t_step, x_step] = ode45(@(t, x) [gauss_planetary_eqn(f0_modified_equinoctial(x, 1), B_modified_equinoctial(x, 1), a_control(x) + options.a_disturbance(t, x .* nd_scalar)); mdot], t_nd(end) + [0, dt_step], x_me_mass_nd(:, end), tolerances);

    t_nd = [t_nd; t_step(2:end)];
    x_me_mass_nd = [x_me_mass_nd, x_step(2:end, :)'];
    
    % Wrap up iteration
    iter = iter + 1;
    iter
end

% Package results
transfer = [];
transfer.converged = Q(end) < Q_stop;
if transfer.converged
    transfer.errors = "";
else % Tranfser failed...
    transfer.errors = 'Errors:';
    if iter == options.iter_max
        transfer.errors = strjoin({transfer.errors, 'MAX_ITER'});
    end
    if x_me_mass_nd(7, end) <= m_dry_nd
        transfer.errors = strjoin({transfer.errors, 'FUEL_DEPLETED'});
    end
end

% Calculate minimum amount of information to evaluate performance of
% Q_params for global optimization of the time-fuel pareto front
transfer.dt = t_nd(end) * char_star.t - 0;
transfer.delta_m = spacecraft_params.m_0 - x_me_mass_nd(7, end) * char_star.m;

transfer.Q = Q;
transfer.P = P;

% Dimensionalize results (and rotate back by Theta_rot)
if options.return_dt_dm_only == false
    x_me_mass_unrotated = x_me_mass_nd .* nd_scalar;
    transfer.x_keplerian_mass = [modified_equinoctial_to_keplerian_array(x_me_mass_unrotated(1:6, :)); x_me_mass_unrotated(7, :)];
    transfer.x_keplerian_mass(4, :) = transfer.x_keplerian_mass(4, :) - Q_params.Theta_rot;
    transfer.u = u_nd * char_star.F;
    transfer.t = t_nd * char_star.t;
    transfer.delta_V = spacecraft_params.Isp * g_0 * log(spacecraft_params.m_0 / x_keplerian_mass(7, end));
    
    transfer.alpha = alpha;
    transfer.beta = beta;
    
    transfer.Q = Q;
    transfer.P = P;
end
end
