%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Targeted Debris Visualization
% Author: Travis Hastreiter 
% Created On: 1 April, 2026
% Description: Use Matlab Aerospace toolkit to visualize debris orbits from
% TLEs
% Most Recent Change: 1 April, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create Mission
mission.StartDate = datetime(2026, 1, 3, 12, 0, 0);
duration_yr = 6;
mission.Duration  = years(duration_yr);
mission.StopDate = mission.StartDate + mission.Duration;
sampleTime = 60000;
sc = satelliteScenario(mission.StartDate, mission.StopDate, sampleTime);

% Load Debris
debris_ID = [20791, 28480, 31114, 32063, 37766, 39203, 39261, 41858, 44548];
debris_mass = [2000, 4000, 4000, 2000, 4000, 4000, 2000, 4000, 4000];
N_debris = numel(debris_ID);
tleFiles = strings(size(debris_ID));
for i = 1 : N_debris
    tleFiles(i) = sprintf("%g.tle", debris_ID(i));
end

debris_1 = satellite(sc, tleFiles(1), ...
                           "Name",sprintf("debris_%g", debris_ID(1)), ...
                           "OrbitPropagator","sgp4");
debris = repmat(debris_1, 1, N_debris);
for i = 2 : N_debris 
    debris(i) = satellite(sc, tleFiles(i), ...
                           "Name",sprintf("debris_%g", debris_ID(i)), ...
                           "OrbitPropagator","sgp4");
end

% Visualize
%v = satelliteScenarioViewer(sc);
%play(sc)

%%
char_star = load_charecteristic_values_Earth;
mu_E = char_star.mu;
R_E = char_star.l;

% Calculate RAAN drift
drift = zeros(size(debris_ID));
a = zeros(size(debris_ID));
for i = 1 : N_debris
    d_oe = debris(i).orbitalElements;
    a(i) = ((d_oe.Period / (2 * pi)) ^ 2 * mu_E) ^ (1 / 3);
    drift(i) = J2_RAAN_drift(a(i), d_oe.Eccentricity, deg2rad(d_oe.Inclination), mu_E, R_E);
end

[r, v] = states(debris_1);
x = [r; v] / 1000;
x_keplerian = cartesian_to_keplerian_array(x, [0; 0; 1], [1; 0; 0], mu_E);
x_keplerian_array = zeros(6, size(x_keplerian, 2), N_debris);
x_keplerian_array(:, :, 1) = x_keplerian;
for i = 2 : N_debris
    [r, v] = states(debris(i));
    x = [r; v] / 1000;
    x_keplerian_array(:, :, i) = cartesian_to_keplerian_array(x, [0; 0; 1], [1; 0; 0], mu_E);
end

%% Compare relative RAAN drift
num_transfers_per_debris = zeros(size(debris_ID, 2), size(debris_ID, 2));
delta_RAAN_step = 30; % [deg] how big RAAN difference has to be before there needs to be another data point

for i = 1 : N_debris
    figure
    plot(linspace(0, duration_yr, size(x_keplerian, 2)), squeeze(rad2deg(unwrap(x_keplerian_array(4, :, :)) - unwrap(x_keplerian_array(4, :, i)))))
    relative_shifts = squeeze(rad2deg(unwrap(x_keplerian_array(4, :, :)) - unwrap(x_keplerian_array(4, :, i))));
    relative_RAAN_range = relative_shifts(1, :) - relative_shifts(end, :);
    num_transfers_per_debris(i, :) = ceil(abs(relative_RAAN_range) / delta_RAAN_step) + 1;
end
num_transfers_per_debris = num_transfers_per_debris - eye(N_debris);

total_transfers = sum(num_transfers_per_debris, "all");

%% Calculate deorbit transfers for each debris

char_star = load_charecteristic_values_Earth();

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 4155; % [s]
spacecraft_params.m_0 = 2000; % [kg] Of spacecraft
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.235; % [N]

% Integration error tolerance
default_tolerance = 1e-10;

a_d_0 = @(t, x) zeros([3, 1]); % Disturbance function

% Min Periapsis soft constraint
penalty_params = struct();
penalty_params.k = 100; % Smoothing parameter
penalty_params.W_p = 1; % Penalty weight
penalty_params.r_p_min = R_E + 0; % [km] min periapsis

% Parameters for the optimization needed to determine efficiencies
Qdot_opt_params = struct();
Qdot_opt_params.num_start_points = 10;
Qdot_opt_params.strategy = "Best Start Points";
Qdot_opt_params.plot_minQdot_vs_L = false;

deorbit_transfer_drifts = zeros([N_debris, 1]);
dVs_deorbit = zeros([N_debris, 1]);
ToFs_deorbit = zeros([N_debris, 1]);

clear Qtransfer_deorbit
parfor i = 1 : N_debris
    % Define Q-Law feedback controller: W_oe, eta_a_min, eta_r_min, m, n, r, Theta_rot
    Q_params = struct();
    Q_params.W_oe = 1 * ones([5, 1]); % Element weights 
    Q_params.eta_a_min = 0.7; % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = 0.7; % Minimum relative efficiency for thrusting instead of coasting
    Q_params.m = 3;
    Q_params.n = 4;
    Q_params.r = 2;
    Q_params.Theta_rot = 0;

    r_a_deorbit = x_keplerian_array(1, 1, i) * (1 + x_keplerian_array(2, 1, i)); % [km] apoapsis
    r_p_deorbit = R_E + 135; % [km] periapsis
    e_deorbit = (1 - r_p_deorbit / r_a_deorbit) / (1 + r_p_deorbit / r_a_deorbit); % [] eccentricity
    a_deorbit = r_p_deorbit / (1 - e_deorbit); % [km] semi-major axis

    x_deorbit_keplerian = x_keplerian_array(:, 1, i);
    x_deorbit_keplerian(1:2) = [a_deorbit; e_deorbit];
    spacecraft_params_i = spacecraft_params;
    spacecraft_params_i.m_0 = spacecraft_params.m_0 + debris_mass(i);
    [Qtransfer_deorbit(i)] = QLaw_transfer_fast(x_keplerian_array(:, 1, i), x_deorbit_keplerian, mu_E, spacecraft_params_i, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 1500000, angular_step=deg2rad(2), thrust_during_eclipse = true, integration_tolerance=default_tolerance);

    deorbit_initial_drifts(i) = sum(J2_RAAN_drift(Qtransfer_deorbit(i).x_keplerian_mass(1, 1), Qtransfer_deorbit(i).x_keplerian_mass(2, 1), Qtransfer_deorbit(i).x_keplerian_mass(3, 1), mu_E, R_E) .* [diff(Qtransfer_deorbit(i).t)', 0]);
    deorbit_transfer_drifts(i) = sum(J2_RAAN_drift(Qtransfer_deorbit(i).x_keplerian_mass(1, :), Qtransfer_deorbit(i).x_keplerian_mass(2, :), Qtransfer_deorbit(i).x_keplerian_mass(3, :), mu_E, R_E) .* [diff(Qtransfer_deorbit(i).t)', 0]);
    dVs_deorbit(i) = Qtransfer_deorbit(i).delta_V;
    ToFs_deorbit(i) = Qtransfer_deorbit(i).dt;

    if Qtransfer_deorbit(i).converged
        fprintf("Q-Law Transfer Converged! Took %.3f Days Using %.3f kg Propellant\n", Qtransfer_deorbit(i).dt / 60 / 60 / 24, Qtransfer_deorbit(i).delta_m)
    else
        fprintf("Q-Law Transfer Failed with %s\n", Qtransfer_deorbit(i).errors)
    end
end

% save("Multi Debris Mission Optimization\deorbit_transfers_info.mat", "ToFs_deorbit", "dVs_deorbit", "deorbit_transfer_drifts")

%% Construct Dataset Inputs
% For each delta RAAN point for two debris % save the starting terminator
% orbit and the next debris orbit to call deorbit_to_new_debris
t_transfer = zeros([1, total_transfers]);
x_initial_keplerian = zeros([6, total_transfers]);
x_final_keplerian = zeros([6, total_transfers]);
IDs_transfer = zeros([2, total_transfers]);

transfer_count = 1;
for i = 1 : N_debris
    for j = 1 : N_debris
        if i ~= j
            total_tspan = [ToFs_deorbit(i) / 60 / 60 / 24 / 365.25, duration_yr];
            t_k = linspace(total_tspan(1), total_tspan(2), num_transfers_per_debris(i, j));
            t_k_index = interp1(linspace(0, duration_yr, size(x_keplerian, 2)), 1 : size(x_keplerian, 2), t_k, "nearest");

            for k = 1 : num_transfers_per_debris(i, j)
                r_a_deorbit = R_E + 500; % [km] apoapsis - CHECK DRAG DURING DEORBIT FOR BETTER ESTIMATE FOR SOME RCS USED
                r_p_deorbit = R_E + 135; % [km] periapsis
                e_deorbit = (1 - r_p_deorbit / r_a_deorbit) / (1 + r_p_deorbit / r_a_deorbit); % [] eccentricity
                a_deorbit = r_p_deorbit / (1 - e_deorbit); % [km] semi-major axis
            
                x_deorbit_keplerian = x_keplerian_array(:, t_k_index(k), i);
                x_deorbit_keplerian(1:2) = [a_deorbit; e_deorbit];

                x_initial_keplerian(:, transfer_count) = x_deorbit_keplerian;
                x_final_keplerian(:, transfer_count) = x_keplerian_array(:, t_k_index(k), j);

                t_transfer(transfer_count) = t_k(k);
                IDs_transfer(:, transfer_count) = [i; j];
                transfer_count = transfer_count + 1;
            end
        end
    end
end

transfer_dataset_inputs = struct();
transfer_dataset_inputs.t0 = t_transfer;
transfer_dataset_inputs.IDs = IDs_transfer;
transfer_dataset_inputs.x_initial_keplerian = x_initial_keplerian;
transfer_dataset_inputs.x_final_keplerian = x_final_keplerian;
transfer_dataset_inputs.debris_mass = debris_mass;
transfer_dataset_inputs.debris_ID = debris_ID;
transfer_dataset_inputs.spacecraft_params = spacecraft_params;

save("Multi Debris Mission Optimization\transfer_dataset_inputs_fixedreorbit.mat", "transfer_dataset_inputs");
