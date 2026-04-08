%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Generate dV-ToF Paretos for transfers from a deorbit orbit (after 
% deorbiting debris) to a set of new debris
% Author: Travis Hastreiter 
% Created On: 15 March, 2026
% Description: Orbit transfer using Q-Law from deorbit orbit (after drop 
% off) to new debris not accounting for rendezvous (assuming not much extra 
% delta V and time). Does it between a set of debris for purpose of using
% to optimize spacecraft routing
% Hybrid because fmincon is used to improve solution returned by GA
%
% Most Recent Change: 15 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Load all Paretos
% 2) Calculate representative deorbit transfer
% 3) Define routing problem as mixed integer problem for genetic algorithm
% 4) Solve
% 5) Save and analyze results

% Objective: 
%   * Maximize deorbited objects
%   * Minimize avg fuel
%   * Make fuel use not vary too much - maybe not that important
% Constraints:
%   * One ship per debris
%   * Maximum total time
%   * Maximum total Delta V
%   * ToF within the pareto bounds (mostly care about mininimum)
% Variables: 
%   - Have max number of debris per ship (N_debris_max) and dummy debris ID 
%     which takes 0 dV, ToF to get to to allow varying # of debris
%   - Total # variables = # ships * N_debris_max * (1 (IDs) + 1 (ToFs))
%   * Debris IDs: (discrete) one for each debris + one for not going to a debris
%   * ToFs: (continuous) use to interpolate Pareto to get Delta V

% Create function that adjusts results based on changing initial mass

%% Initialize Problem
% Load deorbit transfer info
deorbit_info = load("deorbit_transfers_info.mat");
dV_deorbit = deorbit_info.dVs_deorbit;

% Load Pareto creation info
transfer_dataset_inputs = load("Multi Debris Mission Optimization\transfer_dataset_inputs.mat").transfer_dataset_inputs;
debris_IDs = transfer_dataset_inputs.debris_ID;
debris_weights = ones(size(debris_IDs)); % Use McKnight top 50 list score? - all the same??

% Load debris paretos
% Pareto has .ToF and .dV like [N_pareto, N_debris, N_debris] where 3rd dim 
% is starting debris ID and 4rth dim is ending debris ID.
paretos = load("Mission Design\Multi Debris Mission Optimization\Deorbit to Debris Paretos\Low Thrust\low_thrust_paretos.mat").paretos;

%%
% Problem parameters
N_ships = 3;
N_debris = numel(debris_IDs);
N_debris_max = 3;
N_IDs = N_debris + 1;
N_vars = N_ships * (2 * N_debris_max - 1);
max_dV = 10; % [km / s]
max_t = 5.25; % [yr]
J_weights = [1, ... % Debris weight left
             0.1]; % Avg dV

% Variable layout
ship_i = 1:N_ships;
debris_i = 1:N_debris_max;
var_i = 1:2; % 1 is ID, 2 is ToF
var_layout_table = combinations(ship_i, debris_i, var_i); % Every row is a variable
% Get rid of last ToF because there is one ToF for every index pair (transfer)
var_layout_table = var_layout_table(~(var_layout_table.debris_i == N_debris_max & var_layout_table.var_i == 2), :);
integer_vars_i = find(var_layout_table.var_i == 1); % Every ID is integer valued
% Variable bounds
IDs = 0:N_debris;
ID_weights = [0, debris_weights];
ToF_bounds = [min([paretos.ToF{:}], [], "all"), max([paretos.ToF{:}], [], "all")]; % [yr]
var_bounds = [min(IDs), max(IDs);
              ToF_bounds];
lb = var_bounds(var_layout_table.var_i, 1);
ub = var_bounds(var_layout_table.var_i, 2);

% Optimization parameters
opts = optimoptions(@ga, ...
                    'PopulationSize', 5000, ...
                    'MaxGenerations', 500, ...
                    'EliteCount',1000, ...
                    'FunctionTolerance', 1e-12, ...
                    'PlotFcn', @gaplotbestf);

%% Solve
rng(0, 'twister');
[xbest, fbest, exitflag] = ga(@(x) spacecraft_routing_objective(x, var_layout_table, ID_weights, J_weights, paretos, N_ships), N_vars, [], [], [], [], ...
    lb, ub, @(x) spacecraft_routing_nonlconstraints(x, var_layout_table, paretos, max_t, max_dV, N_debris, N_ships, N_debris_max, dV_deorbit, transfer_dataset_inputs.spacecraft_params.m_0, transfer_dataset_inputs.debris_mass'), integer_vars_i, opts);

%%
[IDs_best, ToFs_best, dVs_best] = extract_routing_info(xbest, var_layout_table, paretos, N_ships);
[c, ceq, dV_per_sc, t_per_sc] = spacecraft_routing_nonlconstraints(xbest, var_layout_table, paretos, max_t, max_dV, N_debris, N_ships, N_debris_max, dV_deorbit, transfer_dataset_inputs.spacecraft_params.m_0, transfer_dataset_inputs.debris_mass');

%%
debris_IDs_best = debris_IDs(IDs_best);

%%
ToFs_best_2 = zeros(size(ToFs_best));
dVs_best_2 = zeros(size(dVs_best));
N_pareto = size(paretos.dV, 1);
for s = 1 : N_ships
    N_debris = nnz(IDs_best(s, :));
    N_transfers = N_debris - 1;
    ToF_deorbit = 0;
    for t = 1 : N_debris
        if IDs_best(s, t) ~= 1
            ToF_deorbit = ToF_deorbit + paretos.t{IDs_best(s, t), 1}(1);
        else
            ToF_deorbit = ToF_deorbit + paretos.t{IDs_best(s, t), 2}(1);
        end
    end
    [ToFs_best_2(s, 1:N_transfers), dVs_best_2(s, 1:N_transfers)] = optimize_transfer_ToFs_time_varying(paretos, IDs_best(s, :), max_t - ToF_deorbit * N_debris, ToFs_best(s, :)');
end
num_debris_per_sc = sum(IDs_best ~= 0, 2);
dV_per_sc_2 = sum(dVs_best_2, 2) + num_debris_per_sc * dV_deorbit;
t_per_sc_2 = sum(ToFs_best_2, 2) + num_debris_per_sc * ToF_deorbit;

%% Helper Functions
function [c, ceq, dV_per_sc, t_per_sc] = spacecraft_routing_nonlconstraints(x, var_layout_table, paretos, max_t, max_dV, N_debris, N_ships, N_debris_max, dV_deorbit, spacecraft_mass, debris_mass)
    % Constraints:
    %   * One ship per debris (ineq) - not equality because sum debris
    %     might not be visited so number is 0 or 1
    %   * Maximum total time (ineq)
    %   * Maximum total Delta V (ineq)
    %   * ToF within the Pareto bounds (mostly care about mininimum) (ineq)
    % All inequality constraints are written so satisfied when <= 0

    [IDs, ToFs, dVs] = extract_routing_info(x, var_layout_table, paretos, N_ships);
    num_debris_per_sc = sum(IDs ~= 0, 2);

    % Make debris matrix with #rows = #ships and #cols = #debris storing
    % which debris ships went to so ToF and deorbit dV per ship can be
    % found easily
    debris_matrix = zeros([N_ships, N_debris]);
    for s = 1 : N_ships
        for i = 1 : N_debris
            debris_matrix(s, i) = any(IDs(s, :) == i);
        end
    end

    % One ship per debris - N_debris
    % Count how many times each ID shows up
    c_one_ship_debris = arrayfun(@(x) sum(IDs == x, "all"), (1:N_debris)') - 1; % Skips IDs with 0

    % Max total time - N_ships
    ToF_deorbit = zeros([N_ships, 1]);
    for s = 1 : N_ships
        employed_i = find(IDs(s, :) ~= 0);  % Should always be in order
        for t = 1 : numel(employed_i)
            % First element of .t is the deorbit time which doesn't depend
            % on destination
            if IDs(s, t) ~= 1
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 1}(1);
            else
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 2}(1);
            end
        end
    end
    t_per_sc = sum(ToFs, 2) + ToF_deorbit;
    ToF_deorbit = zeros([N_ships, 1]);
    for s = 1 : N_ships
        employed_i = find(IDs(s, :) ~= 0);  % Should always be in order
        for t = 1 : numel(employed_i)
            % First element of .t is the deorbit time which doesn't depend
            % on destination
            if IDs(s, t) ~= 1
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 1}(1);
            else
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 2}(1);
            end
        end
    end
    t_per_sc = sum(ToFs, 2) + ToF_deorbit;
    c_max_t = t_per_sc - max_t;

    % Max total Delta V - N_ships (Get approximate adjusted delta V as if
    % the spacecraft's mass never changed from the debris, should be within a few percent)
    dV_per_sc = sum(dVs, 2) + debris_matrix * (dV_deorbit .* (spacecraft_mass + debris_mass) ./ spacecraft_mass);
    c_max_dV = dV_per_sc - max_dV;
    
    % ToF within Pareto bounds - 2 * N_ships * (N_debris_max - 1)
    % NEED TO INTERPOLATE
    c_pareto_min_ToF = zeros([N_ships, N_debris_max - 1]);
    c_pareto_max_ToF = zeros([N_ships, N_debris_max - 1]);
    for s = 1 : N_ships
        t0 = 0;
        employed_i = find(IDs(s, :) ~= 0);  % Should always be in order
        for t = 1 : numel(employed_i) - 1
            if IDs(s, t) ~= IDs(s, t + 1)
                ID1 = IDs(s, t); % From 
                ID2 = IDs(s, t + 1); % To
                t0 = t0 + paretos.t{ID1, ID2}(1);

                % Calculate starting time of transfer
                if t ~= 1
                    ToF_total_t = sum(ToFs(s, 1 : (t - 1)));
                else
                    ToF_total_t = 0;
                end
                t0_transfer = min(ToF_total_t + t0, paretos.t{ID1, ID2}(end));

                c_pareto_min_ToF(s, t) = interp1(paretos.t{ID1, IDs(s, t + 1)}, paretos.ToF{ID1, ID2}(1, :), t0_transfer) - ToFs(s, t);
                c_pareto_max_ToF(s, t) = ToFs(s, t) - interp1(paretos.t{ID1, ID2}, paretos.ToF{ID1, ID2}(end, :), t0_transfer);
            end
        end
    end
    c_pareto_ToF = [c_pareto_min_ToF(:); c_pareto_max_ToF(:)];
    
    % Package constraints
    c = [c_one_ship_debris;
         c_max_t; 
         c_max_dV;
         c_pareto_ToF];
    ceq = []; % No equality constraints
end

function [J] = spacecraft_routing_objective(x, var_layout_table, ID_weights, J_weights, paretos, N_ships)
    % Objective: 
    %   * Maximize deorbited objects
    %   * Minimize avg fuel - might be best as a postprocess so that it
    %     doesn't interfere with the optimizer squeezing more debris in
    %   * Make fuel use not vary too much - maybe not that important

    [IDs, ToFs, dVs] = extract_routing_info(x, var_layout_table, paretos, N_ships);

    num_debris_per_sc = sum(IDs ~= 0, 2);
    employed_sc = find(num_debris_per_sc > 0);
    dV_per_sc = sum(dVs, 2);

    % Count number of deorbited debris multiplied by their weights and get
    % what is left (minimize leftover weight)
    debris_deorbited = arrayfun(@(x) any(IDs == x, "all"), (0 : (numel(ID_weights)) - 1)');
    J_deorbit = sum(ID_weights) - sum(ID_weights * debris_deorbited);

    % Calculate average dV per debris
    if numel(employed_sc) > 0
        J_avgdv = max(sum(dV_per_sc(employed_sc) ./ num_debris_per_sc(employed_sc)) / numel(employed_sc), 0);
    else
        J_avgdv = 0;
    end

    % Calculate mission time for each spacecraft
    ToF_deorbit = zeros([N_ships, 1]);
    for s = 1 : N_ships
        employed_i = find(IDs(s, :) ~= 0);  % Should always be in order
        for t = 1 : numel(employed_i)
            % First element of .t is the deorbit time which doesn't depend
            % on destination
            if IDs(s, t) ~= 1
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 1}(1);
            else
                ToF_deorbit(s) = ToF_deorbit(s) + paretos.t{IDs(s, t), 2}(1);
            end
        end
    end
    t_per_sc = sum(ToFs, 2) + ToF_deorbit;

    % Combine objectives
    J = dot(J_weights, [J_deorbit; J_avgdv]) + sum(t_per_sc) * 0.01;
end

function [IDs, ToFs, dVs] = extract_routing_info(x, var_layout_table, paretos, N_ships)
    % Extract raw variables
    IDs_raw = reshape(x(var_layout_table.var_i == 1), N_ships, []);
    ToFs_raw = reshape(x(var_layout_table.var_i == 2), N_ships , []);

    % Push zero indices to last
    IDs = zeros(size(IDs_raw));
    ToFs = zeros(size(ToFs_raw));
    for s = 1 : size(ToFs, 1)
        reordered_i = [find(IDs_raw(s, :) ~= 0), find(IDs_raw(s, :) == 0)];
        IDs(s, :) = IDs_raw(s, reordered_i);
        ToFs(s, 1:min(numel(find(IDs_raw(s, :) ~= 0)), size(ToFs, 2))) = ToFs_raw(s, 1:min(numel(find(IDs_raw(s, :) ~= 0)), size(ToFs, 2)));
    end

    % Interploate Pareto curves to get dV
    dVs = interp_paretos_time_varying(IDs, ToFs, paretos);
end

function [dVs] = interp_paretos_time_varying(IDs, ToFs, paretos)
    % Need to flatten vars to make interp1 happy....
    % Use IDs into paretos...
    % IDs has shape [N_ships, N_debris_max]
    % ToFs has shape [N_ships, N_debris_max - 1]
    % dVs has shape [N_ships, N_debris_max - 1]
    % Pareto has .ToF and .dV like {N_debris, N_debris}[N_pareto, N_times] 
    % and .t like {N_debris, N_debris}[N_times]

    dVs = zeros(size(ToFs));
    for s = 1 : size(ToFs, 1)
        t0 = 0;
        for t = 1 : size(ToFs, 2)
            if IDs(s, t + 1) ~= 0 && IDs(s, t) ~= IDs(s, t + 1)
                ID1 = IDs(s, t); % From 
                ID2 = IDs(s, t + 1); % To
                t0 = t0 + paretos.t{ID1, ID2}(1);

                % Calculate starting time of transfer
                if t ~= 1
                    ToF_total_t = sum(ToFs(s, 1 : (t - 1)));
                else
                    ToF_total_t = 0;
                end
                t0_transfer = min(ToF_total_t + t0, paretos.t{ID1, ID2}(end));

                % Interplolate points
                interpolated_paretos_ToF = interp1(paretos.t{ID1, ID2}, paretos.ToF{ID1, ID2}', t0_transfer)';
                interpolated_paretos_dV = interp1(paretos.t{ID1, ID2}, paretos.dV{ID1, ID2}', t0_transfer)';

                % Bound ToFs
                ToF_bounded = max(min(ToFs(s, t), interpolated_paretos_ToF(end)), interpolated_paretos_ToF(1));
                
                % Interpolate Pareto
                dVs(s, t) = interp1(interpolated_paretos_ToF, interpolated_paretos_dV, ToF_bounded, "linear");
            else
                break;
            end
        end
    end
end