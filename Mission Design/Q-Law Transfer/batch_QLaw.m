function [dV_ToF] = batch_QLaw(basic_vars, x0_d_keplerian, x0_c_keplerian, spacecraft_params, penalty_params, Qdot_opt_params, options)
%BATCH_QLAW Summary of this function goes here
%   Batch run of QLaw for optimization of parameters
arguments
    basic_vars % (N, 11) W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
    x0_d_keplerian % (N, 6)
    x0_c_keplerian % (N, 6)
    spacecraft_params
    penalty_params
    Qdot_opt_params
    options.mu = 398600 % [km3 / s2] Earth gravitational parameter
    options.iter_max = 50000
    options.angular_step = deg2rad(20)
end

N_i = size(basic_vars, 1);
dV_ToF = zeros([N_i, 2]);

parfor i = 1 : N_i
    % Define Q-Law feedback controller: W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
    Q_params = struct();
    Q_params.W_oe = basic_vars(i, 1:5)'; % Element weights 
    Q_params.eta_a_min = basic_vars(i, 6); % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = basic_vars(i, 7); % Minimum relative efficiency for thrusting instead of coasting
    Q_params.m = basic_vars(i, 8);
    Q_params.n = basic_vars(i, 9);
    Q_params.r = basic_vars(i, 10);
    Q_params.Theta_rot = basic_vars(i, 11);

    Qtransfer = QLaw_transfer(x0_d_keplerian(:, i), x0_c_keplerian(:, i), options.mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = true, iter_max = options.iter_max, angular_step = options.angular_step);
    
    if Qtransfer.converged
        dV_ToF(i, :) = [Qtransfer.delta_V, Qtransfer.dt / 60 / 60 /24];
    else
        dV_ToF(i, :) = [nan, nan]; %[Qtransfer.delta_V, Qtransfer.dt / 60 / 60 /24] * 2;
    end
end

end