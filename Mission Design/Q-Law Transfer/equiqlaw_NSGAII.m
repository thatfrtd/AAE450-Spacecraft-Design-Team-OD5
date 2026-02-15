%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Q-Law Orbit Transfer Pareto Front Optimization 
% Author: Travis Hastreiter 
% Created On: 14 February, 2026
% Description: Orbit transfer using modified equinoctial elements Q-Law 
% optimized using NSGA-II genetic algorithm to produce optimal fuel-ToF 
% pareto front 
% Most Recent Change: 15 February, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optimization Variables
W_oe_bounds = [1, 10]; % [1, 10] ? Element scaling
eta_a_min_bounds = [0, 1]; % (0, 1) Minimum absolute efficiency
eta_r_min_bounds = 1; % (0, 1) Minimum relative efficiency
m = 3; % (1, 5) ? S_a scaling parameter
n = 4; % (1, 5) ? S_a scaling parameter
r = 2; % (1, 5) ? S_a scaling parameter
Theta_rot = 0; % [0, 2pi) 

%
MultiObj.fun = @(x) [f1(x(:,1),x(:,2)), f2(x(:,1),x(:,2))];
MultiObj.nVar = 2;
MultiObj.var_min = -pi.*ones(1,MultiObj.nVar);
MultiObj.var_max = pi.*ones(1,MultiObj.nVar);
MultiObj.truePF = PF;

% Parameters
params.Np = 50;        % Population size
params.pc = 0.9;        % Probability of crossover
params.pm = 0.5;        % Probability of mutation
params.maxgen = 50;    % Maximum number of generations
params.ms = 0.05;       % Mutation strength

% NSGA-II algorithm
NSGAII(params, MultiObj);   