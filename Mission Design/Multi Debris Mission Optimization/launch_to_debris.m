%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Transfer from a launch insertion orbit to first debris
% Author: Travis Hastreiter 
% Created On: 11 April, 2026
% Description: Orbit transfer using Q-Law from the orbit the launch vehicle 
% inserted us into to new debris not accounting for rendezvous
% Most Recent Change: 11 April, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main design parameters: transfer efficiency (assuming constant)
% Constraints: max delta V, max ToF
% Need to look at how changing initial time for routing changes results...
% Optimization inner loop
% 1. Calculate transfer from initial (a, i) to target debris (a, i)
% 2. Use transfer time to get drift time
% 3. Repeat 1-2 for other debris

%% Load Info
routing_solution = load("Multi Debris Mission Optimization\routing_solution_3x3_t0of1yr.mat").routing_solution;
N_ships = numel(routing_solution.t_per_sc);
dV_per_sc = routing_solution.dV_per_sc;
t_per_sc = routing_solution.t_per_sc;

transfer_dataset_inputs = load("Multi Debris Mission Optimization\transfer_dataset_inputs_fixedreorbit.mat").transfer_dataset_inputs;
debris_IDs = transfer_dataset_inputs.debris_ID;

x0_keplerian_debris = load("Debris Orbits\initial_debris_keplerian_orbits.mat").x0_keplerian_debris;
x_keplerian_targ = x0_keplerian_debris(:,  routing_solution.debris_i(:, 1));

%%

R_E = 6378.137; % [km] Earth radius
mu_E = 398600.4418; % [km3 / s2] Earth gravitational parameter
J_2_val = 1.08262668e-3; % [] Earth J2

char_star = load_charecteristic_values_Earth();

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 4155; % [s]
spacecraft_params.m_0 = 2000; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.235; % [N]

% Min Periapsis soft constraint
penalty_params = struct();
penalty_params.k = 100; % Smoothing parameter
penalty_params.W_p = 1; % Penalty weight
penalty_params.r_p_min = R_E + 120; % [km] min periapsis

% Define Q-Law feedback controller: W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
Q_params = struct();
Q_params.W_oe = 1 * ones([5, 1]); % Element weights 
Q_params.eta_a_min = 0.5; % Minimum absolute efficiency for thrusting instead of coasting
Q_params.eta_r_min = 0.5; % Minimum relative efficiency for thrusting instead of coasting
Q_params.m = 3;
Q_params.n = 4;
Q_params.r = 2;
Q_params.Theta_rot = 0;

% Parameters for the optimization needed to determine efficiencies
Qdot_opt_params = struct();
Qdot_opt_params.num_start_points = 10;
Qdot_opt_params.strategy = "Best Start Points";
Qdot_opt_params.plot_minQdot_vs_L = false;

% Integration error tolerance
default_tolerance = 1e-10;

%% Optimize Insertion Orbit with Transfers that Leverage J2 using Genetic Algorithm 
% Optimization variables
eta1_bounds = [0, 0.5]; % [] initial -> intermediate transfer min efficiency
eta2_bounds = [0, 0.5]; % [] intermediate -> target transfer min efficiencya_int_bounds = ([600, 1400] + R_E) / char_star.l; % [km] 
a_int_bounds = [600, 1400] + R_E; % [km] 
%e_int_bounds = [1e-5, 0.04]; % [] - assume 0 (1e-5)
i_int_bounds = [0.97, 1.03] * deg2rad(98.5); % [rad]
Omega_int_bounds = [min(x_keplerian_targ(4, :)) - 0.2, max(x_keplerian_targ(4, :)) + 0.2];
% omega_int - assume same as original - target almost circular
var_bounds = [a_int_bounds; i_int_bounds; Omega_int_bounds; repmat([eta1_bounds; eta2_bounds; a_int_bounds; i_int_bounds], N_ships, 1)];

max_total_t = 5; % [yr]
max_total_dV = 4; % [km / s]
max_ToF = (max_total_t - t_per_sc) * 0 + 1;
max_dV = max_total_dV - dV_per_sc;

%% Create Pool
p = gcp("nocreate"); % If no pool, do not create new one.
if isempty(p)
    p = parpool(8);
end

%% Optimization parameters
fminconOptions = optimoptions(@fmincon,'Display','iter','PlotFcn',{'optimplotfval','optimplotx'});
options = optimoptions('particleswarm','SwarmSize',50,'Display','iter','MaxStallIterations',6,'MaxIterations',30,'MaxTime',600,'UseParallel',true,'HybridFcn',{@fmincon, fminconOptions});

%% Solve
rng(0, 'twister');
[xbest, fbest, exitflag, output] = particleswarm(@(x) Qlaw_insertion_objective(x, dV_per_sc, x_keplerian_targ, mu_E, R_E, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV, max_ToF), ...
    size(var_bounds, 1), var_bounds(:, 1), var_bounds(:, 2), options);

%% Helper Functions
function [J, Qtransfers_to_int, Qtransfers_to_targ] = Qlaw_insertion_objective(x, dV_per_sc, x_keplerian_targ, mu, R, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV, max_ToF)
    % Objective: max total dV
    
    % Extract variables from x
    [x_keplerian_insertion, x_keplerian_int, eta] = extract_insertion_and_transfer_info(x);

    % Run transfers
    [dVs_ToFs, Qtransfers_to_int, Qtransfers_to_targ] = QLaw_insertion_transfers(x_keplerian_insertion, x_keplerian_int, x_keplerian_targ, eta, mu, R, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV, max_ToF);

    % Calculate objective
    J = max(dVs_ToFs(:, 1) + dV_per_sc);
end

function [x_keplerian_insertion, x_keplerian_int, eta] = extract_insertion_and_transfer_info(x)
    x_keplerian_insertion = [x(1); 1e-5; x(2); x(3); 0; 0];
    transfer_info = reshape(x(4:end), 4, 3);
    N_ships = size(transfer_info, 2);
    x_keplerian_int = [transfer_info(3, :); 1e-5 * ones([1, N_ships]); transfer_info(4, :); x(3) * ones([1, N_ships]); zeros([2, N_ships])];
    eta = transfer_info(1:2, :);
end

function [dVs_ToFs, Qtransfers_to_int, Qtransfers_to_targ] = QLaw_insertion_transfers(x_keplerian_insertion, x_keplerian_int, x_keplerian_targ, eta, mu, R, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV, max_ToF)
    N_ships = size(x_keplerian_int, 2);

    dVs_ToFs = zeros([N_ships, 2]);
    % Call QLaw_J2_drift_transfer for each ship with corresponding
    % eta1, eta2, a, i starting from insertion
    for i = 1 : N_ships
        [dVs_ToFs(i, :), Qtransfers_to_int(i), Qtransfers_to_targ(i)] = QLaw_J2_drift_transfer(x_keplerian_insertion, x_keplerian_int(:, i), x_keplerian_targ(:, i), eta(:, i), mu, R, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV(i), max_ToF(i));
    end
end

function [dV_ToF, Qtransfer_to_int, Qtransfer_to_targ] = QLaw_J2_drift_transfer(x_keplerian_0, x_keplerian_int, x_keplerian_targ, eta, mu, R, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV, max_ToF)
    arguments
        x_keplerian_0
        x_keplerian_int
        x_keplerian_targ
        eta
        mu
        R
        J_2_val
        spacecraft_params
        Q_params
        penalty_params
        Qdot_opt_params
        max_dV
        max_ToF
    end

    % Transfer to intermediate orbit
    Q_params.eta_a_min = eta(1); % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = eta(1); % Minimum relative efficiency for thrusting instead of coasting

    if x_keplerian_0(1) == x_keplerian_int(1)
        x_keplerian_int(1) = x_keplerian_int(1) + 1e-5;
    end

    [Qtransfer_to_int] = QLaw_transfer_fast(x_keplerian_0, x_keplerian_int, mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 50000, angular_step=deg2rad(20), max_t = max_ToF * 60 * 60 * 24 * 365.25, max_dV = max_dV);
    transfer_drift_1 = sum(J2_RAAN_drift(Qtransfer_to_int.x_keplerian_mass(1, :), Qtransfer_to_int.x_keplerian_mass(2, :), Qtransfer_to_int.x_keplerian_mass(3, :), mu, R, J_2_val) .* [diff(Qtransfer_to_int.t)', 0]);

    % Transfer to target orbit
    Q_params.eta_a_min = eta(2); % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = eta(2); % Minimum relative efficiency for thrusting instead of coasting
    
    spacecraft_params.m_0 = spacecraft_params.m_0 - Qtransfer_to_int.delta_m;
        
    [Qtransfer_to_targ] = QLaw_transfer_fast(x_keplerian_int, [x_keplerian_targ(1:3); x_keplerian_0(4); x_keplerian_targ(5:6)], mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 50000, angular_step=deg2rad(20), max_t = max_ToF * 60 * 60 * 24 * 365.25, max_dV = max_dV);
    transfer_drift_2 = sum(J2_RAAN_drift(Qtransfer_to_targ.x_keplerian_mass(1, :), Qtransfer_to_targ.x_keplerian_mass(2, :), Qtransfer_to_targ.x_keplerian_mass(3, :), mu, R, J_2_val) .* [diff(Qtransfer_to_targ.t)', 0]);

    % Calculate wait time for RAAN phasing accounting for drift during transfers
    targ_Omega_transfer_drift = J2_RAAN_drift(x_keplerian_targ(1, :), x_keplerian_targ(2, :), x_keplerian_targ(3, :), mu, R, J_2_val) * (Qtransfer_to_int.dt + Qtransfer_to_targ.dt);
    delta_Omega = wrapTo2Pi((x_keplerian_targ(4) + targ_Omega_transfer_drift) - (x_keplerian_0(4) + transfer_drift_1 + transfer_drift_2));
    rel_Omega_drift = J2_RAAN_drift(x_keplerian_int(1, :), x_keplerian_int(2, :), x_keplerian_int(3, :), mu, R, J_2_val) ...
                    - J2_RAAN_drift(x_keplerian_targ(1, :), x_keplerian_targ(2, :), x_keplerian_targ(3, :), mu, R, J_2_val);
    t_wait = delta_Omega ./ rel_Omega_drift .* (rel_Omega_drift > 0) ...
           + (delta_Omega - 2 * pi) ./ rel_Omega_drift .* (rel_Omega_drift < 0);

    % Package outputs
    dV_total = Qtransfer_to_int.delta_V + Qtransfer_to_targ.delta_V;
    ToF_total = (Qtransfer_to_int.dt + t_wait + Qtransfer_to_targ.dt) / 60 / 60 / 24 / 365.25;

    % Constraints
    n_constraints = 4; % max dV, min dV, transfer 1 converge, transfer 2 converge
    violated_constraints = (dV_total > max_dV) + (ToF_total > max_ToF) + ~Qtransfer_to_int.converged + ~Qtransfer_to_targ.converged;
    dV_ToF = [dV_total + 1e2 * (violated_constraints / n_constraints),... 
              ToF_total + 1e2 * (violated_constraints / n_constraints)];
end
