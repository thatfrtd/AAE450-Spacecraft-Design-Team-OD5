%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Combined Transfer from Insertion Orbit, Target Rendezvous, and Deorbit
% for all 3 Ships
% Author: Travis Hastreiter 
% Created On: 12 April, 2026
% Description: 
% Most Recent Change: 12 April, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 4155; % [s]
spacecraft_params.m_0 = 2900; % [kg] Of spacecraft
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.235; % [N]

%% Load Info
% Load initial debris orbits
x0_keplerian_debris = load("Debris Orbits\initial_debris_keplerian_orbits.mat").x0_keplerian_debris;

% Load deorbit transfer info
deorbit_info = load("deorbit_transfers_info_fixedreorbit.mat");
dV_deorbit = deorbit_info.dVs_deorbit;

% Load terminator to target transfer info
transfer_dataset_inputs = load("Multi Debris Mission Optimization\transfer_dataset_inputs_fixedreorbit.mat").transfer_dataset_inputs;
debris_IDs = transfer_dataset_inputs.debris_ID;

% Load vehicle routing solution
routing_solution = load("Multi Debris Mission Optimization\routing_solution_3x3_t0of1yr.mat").routing_solution;
N_ships = numel(routing_solution.t_per_sc);
dV_per_sc = routing_solution.dV_per_sc;
t_per_sc = routing_solution.t_per_sc;
x_keplerian_targ = x0_keplerian_debris(:,  routing_solution.debris_i(:, 1));

% Load selected launch insertion orbit
launch_insertion_orbit_transfers = load("Multi Debris Mission Optimization\launch_insertion_orbit_transfers.mat").launch_insertion_orbit_transfers;

%% Get Initial Guess Parameters for Whole Mission
N_Qtransfers = 9;
dVs = zeros([N_ships, N_Qtransfers]); % [km / s]
ToFs = zeros([N_ships, N_Qtransfers]); % [years]
etas = zeros([N_ships, N_Qtransfers]);
x_keplerian_int = cell(3);

% Launch insertion -> 1st target
dVs(:, 1:2) = [[launch_insertion_orbit_transfers.Qtransfers_to_int.delta_V]', ...
               [launch_insertion_orbit_transfers.Qtransfers_to_targ.delta_V]'];
ToFs(:, 1:2) = [[launch_insertion_orbit_transfers.Qtransfers_to_int.dt]', ...
                [launch_insertion_orbit_transfers.Qtransfers_to_targ.dt]'] / 60 / 60 / 24 / 365.25;
etas(:, 1:2) = launch_insertion_orbit_transfers.eta';
x_keplerian_int{1} = launch_insertion_orbit_transfers.x_keplerian_int;

% Deorbit 1st target
dVs(:, 3) = deorbit_info.dVs_deorbit(routing_solution.debris_i(:, 1));
ToFs(:, 3) = deorbit_info.ToFs_deorbit(routing_solution.debris_i(:, 1)) / 60 / 60 / 24 / 365.25;
etas(:, 3) = 0.7 * ones([N_ships, 1]);

% Terminator orbit -> 2nd target
dVs(:, 4:5) = routing_solution.dVs(:, 1) * 0.5 * ones([1, 2]); % REPLACE WITH ACTUAL ESTIMATED VALUES
ToFs(:, 4:5) = routing_solution.ToFs(:, 1) * 0.5 * ones([1, 2]); % REPLACE WITH ACTUAL ESTIMATED VALUES
etas(:, 4:5) = nan([N_ships, 2]); % REPLACE WITH ACTUAL ESTIMATED VALUES
x_keplerian_int{2} = {};

% Deorbit 2nd target
dVs(:, 6) = deorbit_info.dVs_deorbit(routing_solution.debris_i(:, 2));
ToFs(:, 6) = deorbit_info.ToFs_deorbit(routing_solution.debris_i(:, 2)) / 60 / 60 / 24 / 365.25;
etas(:, 6) = 0.7 * ones([N_ships, 1]);

% Terminator orbit -> 3rd target
dVs(:, 7:8) = routing_solution.dVs(:, 2) * 0.5 * ones([1, 2]); % REPLACE WITH ACTUAL ESTIMATED VALUES
ToFs(:, 7:8) = routing_solution.ToFs(:, 2) * 0.5 * ones([1, 2]) / 60 / 60 / 24 / 365.25; % REPLACE WITH ACTUAL ESTIMATED VALUES
etas(:, 7:8) = nan([N_ships, 2]); % REPLACE WITH ACTUAL ESTIMATED VALUES
x_keplerian_int{3} = {};

% Deorbit 3rd target
dVs(:, 9) = deorbit_info.dVs_deorbit(routing_solution.debris_i(:, 3));
ToFs(:, 9) = deorbit_info.ToFs_deorbit(routing_solution.debris_i(:, 3)) / 60 / 60 / 24 / 365.25;
etas(:, 9) = 0.7 * ones([N_ships, 1]);














%% Initialize Satellite Scenario (for visualization)
% Specify scenario
startTime = datetime(2020,5,11,12,35,38);
stopTime = startTime + days(200);
sampleTime = 60*2;
sc = satelliteScenario(startTime,stopTime,sampleTime);

% Load Satellite info
tleFile = "eccentricOrbitSatellite.tle";

sat = satellite(sc, pos, vel, CoordinateFrame="inertial");

%% 
    