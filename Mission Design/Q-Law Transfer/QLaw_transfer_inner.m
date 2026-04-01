function [t_step, x_step, Q, P, alpha, beta, u_nd, not_coast] = QLaw_transfer_inner(x_me_mass_nd, t_nd, oe_t, Q_params_W_oe, Q_params_m, Q_params_n, Q_params_r, Q_params_Theta_rot, Q_params_eta_a_min, Q_params_eta_r_min, penalty_params_W_p, penalty_params_k, F_max_nd, thrust_during_eclipse, r_p_min_nd, nd_scalar, char_star_t, char_star_v, spacecraft_params_Isp, num_start_points, integration_tolerance, angular_step)
%QLAW_TRANSFER Orbit transfer using Q-Law
% Create slow QLaw orbital elements (a, f, g, h, k)
    e = sqrt(x_me_mass_nd(2) ^ 2 + x_me_mass_nd(3) ^ 2);
    a = x_me_mass_nd(1) / (1 - e ^ 2); % a = p / (1 - e ^ 2)
    w = 1 + x_me_mass_nd(2) .* cos(x_me_mass_nd(6)) + x_me_mass_nd(3) .* sin(x_me_mass_nd(6));
    r = x_me_mass_nd(1) / w; % r = p / (1 + f .* cos(L) + g .* sin(L))
    oe = [a; x_me_mass_nd(2:5)];
    L = x_me_mass_nd(6);

    % Calculate control
    [D, Q, Qdot, P, partial_Q_partial_oe] = D_Q_Qdot_P_partial_Q_partial_oe_func(oe, L, oe_t, Q_params_W_oe, Q_params_m, Q_params_n, Q_params_r, F_max_nd, penalty_params_W_p, r_p_min_nd, penalty_params_k);
    [u_nd, alpha, beta] = QLaw_thrust_mapping(D, F_max_nd);

    % Determine if thrusting should happen based on heuristic
    [Qdot_min, Qdot_max] = Qdot_extremize_fast(oe, partial_Q_partial_oe, num_start_points);
    [eta_a, eta_r] = QLaw_efficiencies(Qdot, Qdot_min, Qdot_max);
    % if ~thrust_during_eclipse % Should also add state of charge based thrusting
    %     % Eclipse
    %     R_E = 6378.1; % [km] Earth radius
    %     AU = 149597898; % [km] astronautical unit, Earth-Sun difference
    %     mu_sun = 132712440017.99; % [km3 / s2] Sun gravitational parameter
    %     n_E = sqrt(mu_sun / AU ^ 3); % [rad / s] Earth mean motion
    % 
    %     x_cartesian = x_me_nd_to_cartesian(x_me_mass_nd(1:6), nd_scalar, Q_params_Theta_rot, mu);
    %     rvec = x_cartesian(1:3);
    % 
    %     [~, eclipsed] = check_eclipse(t_nd * char_star_t, rvec, R_E, AU, n_E); % Need to properly get sun position and Earth tilt
    % 
    %     eclipse_no_thrust = eclipsed(iter);
    % else % Thrust during eclispe
    eclipse_no_thrust = false;
    % end
    not_coast = (eta_a >= Q_params_eta_a_min ...
              && eta_r >= Q_params_eta_r_min) ...
              && ~eclipse_no_thrust;
    u_nd = u_nd * not_coast;
    %a_control = @(x) u_nd(:, iter) / x_me_mass_nd(7, end); % a = F / m
    a_control = @(x) u_nd / x(end);
    g_0 = 9.81e-3; % [km / s2]
    mdot = -F_max_nd / (spacecraft_params_Isp * g_0 / char_star_v) * not_coast;

    % Propagate orbit
    dt_step = angular_step * sqrt(r ^ 3 / 1) * sqrt(w) / (1 + e);
    tolerances = odeset(RelTol=integration_tolerance, AbsTol=integration_tolerance, Stats = "off");
    tolerances.InitialStep = dt_step;
    tolerances.MaxStep = dt_step;

    [t_step, x_step] = ode45(@(t, x) [gauss_planetary_eqn(f0_modified_equinoctial(x, 1), B_modified_equinoctial(x, 1), a_control(x)); mdot], t_nd + [0, dt_step], x_me_mass_nd, tolerances);
end