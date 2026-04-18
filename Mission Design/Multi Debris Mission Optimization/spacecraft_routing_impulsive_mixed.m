%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Generate dV-ToF Paretos for transfers from a deorbit orbit (after 
% deorbiting debris) to a set of new debris
% Author: Travis Hastreiter 
% Created On: 15 March, 2026
% Description: Orbit transfer using Q-Law from deorbit orbit (after drop 
% off) to new debris not accounting for rendezvous (assuming not much extra 
% delta V and time). Does it between a set of debris for purpose of using
% to optimize spacecraft routing. 
% Mixed because GA uses fmincon in every fitness evaluation (GA optimizes 
% debris IDs order, fmincon optimizes ToF for each segment)
% 
% Right now it uses impulsive Lambert transfers to make sure the procedure
% works.
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

% ADD DEORBIT ORBIT DV AND TOF AND DRIFT (drift has to be included in pareto)
dV_deorbit = 0.6 * 3; % Mass about 3 times larger with debris so multiply dV by 3 (uses 3 times for fuel)
ToF_deorbit = 240 / 365; % [yr]

%% Initialize Problem

% Load debris paretos
% Pareto has .ToF and .dV like [N_pareto, N_debris, N_debris] where 3rd dim 
% is starting debris ID and 4rth dim is ending debris ID.
paretos = load("C:\Users\thatf\OneDrive\Documents\Purdue Classes\AAE 450\AAE450-Spacecraft-Design-Team-OD5\Mission Design\Multi Debris Mission Optimization\Deorbit to Debris Paretos\Impulsive\fake_paretos.mat").pareto;
% debris_IDs = paretos.debris_IDs;
debris_IDs = 1 : size(paretos.dV, 2);
debris_weights = ones(size(debris_IDs)); % Use McKnight top 50 list score? - all the same??
%% Sort Paretos and Create Bounds (should be done already!!!)
paretos.bounds = zeros(2, numel(debris_IDs), numel(debris_IDs));
for d_1 = 1 : numel(debris_IDs)
    for d_2 = 1 : numel(debris_IDs)
        [ToF_sorted, ToF_sorted_i] = sort(paretos.ToF(:, d_1, d_2));
        for i = 1 : (numel(ToF_sorted) - 1) % Fix duplicate ToFs...
            if ToF_sorted(i) == ToF_sorted(i + 1)
                ToF_sorted(i + 1) = ToF_sorted(i + 1) + 1e-10;
            end
        end
        paretos.ToF(:, d_1, d_2) = ToF_sorted;
        paretos.dV(:, d_1, d_2) = paretos.dV(ToF_sorted_i, d_1, d_2);
        paretos.bounds(:, d_1, d_2) = [min(ToF_sorted); max(ToF_sorted)];
    end
end

%%
% Problem parameters
N_ships = 4;
N_debris = numel(debris_IDs);
N_debris_max = 3;
N_IDs = N_debris + 1;
N_vars = N_ships * N_debris_max;
max_dV = 6.2; % [km / s]
max_t = 8; % [yr]
J_weights = [1, ... % Debris weight left
             0.1]; % Avg dV

% Variable layout
ship_i = 1:N_ships;
debris_i = 1:N_debris_max;
integer_vars_i = 1 : N_vars; % Every ID is integer valued (so every variable)
% Variable bounds
IDs = 0:N_debris;
ID_weights = [0, debris_weights];
var_bounds = [min(IDs), max(IDs)];
lb = repmat(var_bounds(:, 1), 1, N_vars);
ub = repmat(var_bounds(:, 2), 1, N_vars);

% Optimization parameters
opts = optimoptions(@ga, ...
                    'PopulationSize', 200, ...
                    'MaxGenerations', 500, ...
                    'EliteCount',20, ...
                    'FunctionTolerance', 1e-12, ...
                    'PlotFcn', @gaplotbestf);

%% Solve
rng(0, 'twister');
[xbest, fbest, exitflag] = ga(@(x) spacecraft_routing_objective(x, ID_weights, J_weights, paretos, N_ships, max_t, ToF_deorbit), N_vars, [], [], [], [], ...
    lb, ub, @(x) spacecraft_routing_nonlconstraints(x, paretos, max_t, max_dV, N_debris, N_ships, N_debris_max, dV_deorbit, ToF_deorbit), integer_vars_i, opts);

%%
[IDs_best, ToFs_best, dVs_best] = extract_routing_info(xbest, paretos, N_ships, max_t, ToF_deorbit);
[c, ceq, dV_per_sc, t_per_sc] = spacecraft_routing_nonlconstraints(xbest, paretos, max_t, max_dV, N_debris, N_ships, N_debris_max, dV_deorbit, ToF_deorbit);

%% Helper Functions
function [c, ceq, dV_per_sc, t_per_sc] = spacecraft_routing_nonlconstraints(x, paretos, max_t, max_dV, N_debris, N_ships, N_debris_max, dV_deorbit, ToF_deorbit)
    % Constraints:
    %   * One ship per debris (ineq) - not equality because sum debris
    %     might not be visited so number is 0 or 1
    %   * Maximum total time (ineq)
    %   * Maximum total Delta V (ineq)
    %   * ToF within the Pareto bounds (mostly care about mininimum) (ineq)
    % All inequality constraints are written so satisfied when <= 0

    [IDs, ToFs, dVs] = extract_routing_info(x, paretos, N_ships, max_t, ToF_deorbit);
    num_debris_per_sc = sum(IDs ~= 0, 2);

    % One ship per debris - N_debris
    % Count how many times each ID shows up
    c_one_ship_debris = arrayfun(@(x) sum(IDs == x, "all"), (1:N_debris)') - 1; % Skips IDs with 0

    % Max total time - N_ships
    t_per_sc = sum(ToFs, 2) + num_debris_per_sc * ToF_deorbit;
    c_max_t = t_per_sc - max_t;

    % Max total Delta V - N_ships
    dV_per_sc = sum(dVs, 2) + num_debris_per_sc * dV_deorbit;
    c_max_dV = dV_per_sc - max_dV;
    
    % ToF within Pareto bounds - 2 * N_ships * (N_debris_max - 1)
    c_pareto_min_ToF = zeros([N_ships, N_debris_max - 1]);
    c_pareto_max_ToF = zeros([N_ships, N_debris_max - 1]);
    for s = 1 : N_ships
        employed_i = find(IDs(s, :) ~= 0);  % Should always be in order
        for t = 1 : numel(employed_i) - 1
            c_pareto_min_ToF(s, t) = paretos.bounds(1, IDs(s, t), IDs(s, t + 1)) - ToFs(s, t);
            c_pareto_max_ToF(s, t) = ToFs(s, t) - paretos.bounds(2, IDs(s, t), IDs(s, t + 1));
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

function [J] = spacecraft_routing_objective(x, ID_weights, J_weights, paretos, N_ships, max_t, ToF_deorbit)
    % Objective: 
    %   * Maximize deorbited objects
    %   * Minimize avg fuel - might be best as a postprocess so that it
    %     doesn't interfere with the optimizer squeezing more debris in
    %   * Make fuel use not vary too much - maybe not that important

    [IDs, ToFs, dVs] = extract_routing_info(x, paretos, N_ships, max_t, ToF_deorbit);

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

    % Combine objectives
    J = dot(J_weights, [J_deorbit; J_avgdv]);
end

function [IDs, ToFs, dVs] = extract_routing_info(x, paretos, N_ships, max_t, ToF_deorbit)
    % Extract raw variables
    IDs_raw = reshape(x, N_ships, []);

    % Push zero indices to last and get rid of duplicates
    IDs = zeros(size(IDs_raw));
    for s = 1 : size(IDs_raw, 1)
        reordered_i = [find(IDs_raw(s, :) ~= 0), find(IDs_raw(s, :) == 0)];
        IDs(s, :) = IDs_raw(s, reordered_i);
        for t = 1 : (size(IDs_raw, 2) - 1) % Get rid of duplicates
            if IDs(s, t) == IDs(s, t + 1)
                IDs(s, t + 1) = 0;
            end
        end
        % Should not need to repeat this... oh well
        reordered_i = [find(IDs(s, :) ~= 0), find(IDs(s, :) == 0)];
        IDs(s, :) = IDs(s, reordered_i);
        for t = 1 : (size(IDs_raw, 2) - 1) % Get rid of duplicates
            if IDs(s, t) == IDs(s, t + 1)
                IDs(s, t + 1) = 0;
            end
        end
    end

    % Interploate Pareto curves to get dV
    [ToFs, dVs] = interp_paretos(IDs, paretos, max_t, ToF_deorbit);
end

function [ToFs_best, dVs_best] = interp_paretos(IDs, paretos, max_t, ToF_deorbit)
    % Need to flatten vars to make interp1 happy....
    % Use IDs into paretos...
    % IDs has shape [N_ships, N_debris_max]
    % dVs has shape [N_ships, N_debris_max - 1]
    % Pareto has .ToF and .dV like [N_pareto, N_debris, N_debris] where 3rd dim 
    % is starting debris ID and 4rth dim is ending debris ID.

    ToFs_best = zeros(size(IDs, 1), size(IDs, 2) - 1);
    dVs_best = zeros(size(ToFs_best));
    N_pareto = size(paretos.dV, 1);
    for s = 1 : size(IDs, 1)
        N_debris = nnz(IDs(s, :));
        N_transfers = N_debris - 1;
        paretos_ToF = zeros([N_pareto, N_transfers]);
        paretos_dV = zeros([N_pareto, N_transfers]);
        ToF_bounds = zeros([2, N_transfers]);
        for t = 1 : N_transfers
            paretos_ToF(:, t) = paretos.ToF(:, IDs(s, t), IDs(s, t + 1));
            paretos_dV(:, t) = paretos.dV(:, IDs(s, t), IDs(s, t + 1));
            ToF_bounds(:, t) = paretos.bounds(:, IDs(s, t), IDs(s, t + 1));
        end
        if N_transfers > 0
            [ToFs_best(s, 1:N_transfers), dVs_best(s, 1:N_transfers)] = optimize_transfer_ToFs(paretos_ToF, paretos_dV, ToF_bounds, max_t - ToF_deorbit * N_debris);
        end
    end
end
