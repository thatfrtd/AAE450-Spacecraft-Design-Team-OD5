function [transfer] = QLaw_transfer_fast(x_keplerian_d, x_keplerian_c, mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, options)
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
    options.thrust_during_eclipse = true % if an eclipse is detected at the start of the step, should thrust occur?
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
oe_t = [x_keplerian_c(1) / char_star.l; x_me_c_nd(2:5)];

% Get stopping condition
Q_stop = QLaw_stopping_condition(options.R_c, Q_params.W_oe);

% Prepare integration
x_me_mass_nd = [x_me_d_nd; m_0_nd]; % [p, f, g, h, k, L, mass]
u_nd = zeros([3, 1]); % [F_T, F_R, F_N]
alpha = [];
beta = [];
t_nd = 0;
u_cont_nd = zeros([3, 1]);

% Simulate
iter = 1;
Q = [];
P = [];
eclipsed = [];
not_coast = [];
while iter == 1 || Q(iter - 1) >= Q_stop && iter < options.iter_max && x_me_mass_nd(7, end) > m_dry_nd
    [t_step, x_step, Q(end + 1), P(end + 1), alpha(end + 1), beta(end + 1), u_nd(:, end + 1), not_coast(end + 1)] = QLaw_transfer_inner_mex(x_me_mass_nd(:, end), t_nd(end), oe_t, ...
                                                                                                                                            Q_params.W_oe, Q_params.m, Q_params.n, Q_params.r, Q_params.Theta_rot, Q_params.eta_a_min, Q_params.eta_r_min, ...
                                                                                                                                            penalty_params.W_p, penalty_params.k, ...
                                                                                                                                            F_max_nd, options.thrust_during_eclipse, r_p_min_nd, nd_scalar, char_star.t, char_star.v, ...
                                                                                                                                            spacecraft_params.Isp, Qdot_opt_params.num_start_points, options.integration_tolerance, options.angular_step);

    t_nd = [t_nd; t_step(2:end)];
    x_me_mass_nd = [x_me_mass_nd, x_step(2:end, :)'];
    u_cont_nd = [u_cont_nd, repmat(u_nd(:, iter), 1, numel(t_step(2:end)))];

    % fprintf("Iter %g\n", iter)
    % Wrap up iteration
    iter = iter + 1;
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
transfer.delta_V = spacecraft_params.Isp * g_0 * log(spacecraft_params.m_0 / (x_me_mass_nd(7, end) * char_star.m));

% Dimensionalize results (and rotate back by Theta_rot)
if options.return_dt_dm_only == false
    x_me_mass_unrotated = x_me_mass_nd .* nd_scalar;
    transfer.x_keplerian_mass = [modified_equinoctial_to_keplerian_array(x_me_mass_unrotated(1:6, :)); x_me_mass_unrotated(7, :)];
    transfer.x_keplerian_mass(4, :) = transfer.x_keplerian_mass(4, :) - Q_params.Theta_rot;
    transfer.u = u_nd * char_star.F;
    transfer.u_cont = u_cont_nd * char_star.F;
    transfer.t = t_nd * char_star.t;
    
    transfer.not_coast = not_coast;
    transfer.alpha = alpha;
    transfer.beta = beta;
    
    transfer.Q = Q;
    transfer.P = P;

    if ~options.thrust_during_eclipse
        transfer.eclipsed = eclipsed;
    end
end
end

function [x_cartesian] = x_me_nd_to_cartesian(x_me_nd, nd_scalar, Theta_rot, mu)
    x_me_unrotated = x_me_nd .* nd_scalar(1:6);
    x_keplerian = modified_equinoctial_to_keplerian(x_me_unrotated);
    x_keplerian(4) = x_keplerian(4) - Theta_rot;
    x_cartesian = keplerian_to_cartesian(x_keplerian, [], mu);
end