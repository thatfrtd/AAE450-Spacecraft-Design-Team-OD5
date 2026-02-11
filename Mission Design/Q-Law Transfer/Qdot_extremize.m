function [Qdot_min, Qdot_max] = Qdot_extremize(oe, partial_Q_partial_oe, Qdot_opt_params)
%QDOT_EXTREMIZE Summary of this function goes here
%   The minimum and maximum w.r.t. L of the minimum Qdot w.r.t. thrust 
% direction needs to be found so that the current L can be compared to the
% best and worst Ls to compute effeciences used to determine if coasting
% should be done or not. There is usually two minimums and maximums so the
% right ones should try to be chosen. Care must be taken to make this
% process fast.
arguments
    oe
    partial_Q_partial_oe
    Qdot_opt_params % num_start_points, strategy: "Best Start Points", "Multistart", "fminbnd", plot_minQdot_vs_L
end

% Know min w.r.t. alpha, beta of Qdot in terms of L
min_ab_Qdot = @(L) -norm(partial_Q_partial_oe * B_oe_func(oe, L));

% Initialize search for L
L_start = linspace(0, 2 * pi, Qdot_opt_params.num_start_points + 1);
Qdot_start = zeros([Qdot_opt_params.num_start_points, 1]);
for s = 1 : Qdot_opt_params.num_start_points
    Qdot_start(s) = min_ab_Qdot(L_start(s));
end

if Qdot_opt_params.strategy == "Best Start Points"
    % Just look at the min and max of start point array, no iteration

    % Min w.r.t. L (min w.r.t. alpha, beta (Qdot))
    [Qdot_min, L_min_i] = min(Qdot_start);
    L_min = L_start(L_min_i);
    
    % Max w.r.t. L (min w.r.t. alpha, beta (Qdot))
    [Qdot_max, L_max_i] = max(Qdot_start);
    L_max = L_start(L_max_i);
elseif Qdot_opt_params.strategy == "Multistart"
    % Use built-in Multistart with start points to initialize(?)
    % Might be better to do multistart manually

    opts = optimoptions(@fmincon,'Algorithm','sqp');
    ms = MultiStart;

    start_set = CustomStartPointSet(L_start);

    % Min w.r.t. L (min w.r.t. alpha, beta (Qdot))
    problem = createOptimProblem('fmincon','objective',min_ab_Qdot,'x0',0,'lb',0,'ub',2 * pi,'options',opts);
    [L_min, Qdot_min] = run(ms, problem, start_set);

    % Max w.r.t. L (min w.r.t. alpha, beta (Qdot))
    problem.objective = @(L) -min_ab_Qdot(L);
    [L_max, Qdot_max] = run(ms, problem, start_set);
elseif Qdot_opt_params.strategy == "fminbnd"
    % Use fminbnd - will it find the correct one extrema automatically?? 
    % Doesn't use start point info.

    % Min w.r.t. L (min w.r.t. alpha, beta (Qdot))
    [L_min, Qdot_min] = fminbnd(min_ab_Qdot, 0, 2*pi);

    % Max w.r.t. L (min w.r.t. alpha, beta (Qdot))
    [L_max, Qdot_max] = fminbnd(@(L) -min_ab_Qdot(L), 0, 2*pi);
end

% Plot to test how close Qdot_min and Qdot_max are from true values
if Qdot_opt_params.plot_minQdot_vs_L
    L_test = linspace(0, 2 * pi, 360);
    Qdot_test = zeros([numel(L_test), 1]);
    for t = 1 : numel(L_test)
        Qdot_test(t) = min_ab_Qdot(L_test(t));
    end

    clf
    figure(1)
    plot(rad2deg(L_test), Qdot_test, HandleVisibility = "off"); hold on
    scatter(rad2deg(L_start(1:end - 1)), Qdot_start, 36, DisplayName = "Start Points")
    scatter(rad2deg(L_min), Qdot_min, 48, "red", "filled", "v", DisplayName = "Est Min")
    scatter(rad2deg(L_max), Qdot_max, 48, "red", "filled", "^", DisplayName = "Est Max")
    xlabel("True Longitude (L) [deg]")
    ylabel("min Qdot w.r.t. alpha beta")
    xlim([0, 360])
    title("min Qdot w.r.t. Thrusting Direction vs True Longitude")
    grid on
    legend()
end

end