%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Transfer from a deorbit orbit (after deorbiting debris) to new debris
% Author: Travis Hastreiter 
% Created On: 13 March, 2026
% Description: Orbit transfer using Q-Law from deorbit orbit (after drop 
% off) to new debris not accounting for rendezvous (assuming not much extra 
% delta V and time).
% Most Recent Change: 15 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Dataset Inputs
transfer_dataset_inputs = load("Multi Debris Mission Optimization\transfer_dataset_inputs_fixedreorbit.mat").transfer_dataset_inputs;

transfers_i = 1;

N_transfers = numel(transfers_i);

%% Create Pool
p = gcp("nocreate"); % If no pool, do not create new one.
if isempty(p)
    p = parpool(8);
end

for transfer_number = 1 : N_transfers
    %% Initialize Run
    run_name = sprintf("transfer_%g.mat", transfers_i(transfer_number));
    ID1 = transfer_dataset_inputs.debris_ID(transfer_dataset_inputs.IDs(1, transfers_i(transfer_number)));
    ID2 = transfer_dataset_inputs.debris_ID(transfer_dataset_inputs.IDs(2, transfers_i(transfer_number)));

    R_E = 6378.137; % [km] Earth radius
    mu_E = 398600.4418; % [km3 / s2] Earth gravitational parameter
    J_2_val = 1.08262668e-3; % [] Earth J2
    
    % Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
    x0_c_keplerian = transfer_dataset_inputs.x_final_keplerian(:, transfers_i(transfer_number));
    x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, [], mu_E);
    
    % Initial conditions for spacecraft
    x0_d_keplerian = transfer_dataset_inputs.x_initial_keplerian(:, transfers_i(transfer_number));
    x0_d_cartesian = keplerian_to_cartesian(x0_d_keplerian, [], mu_E);
    
    char_star = load_charecteristic_values_Earth();
    
    % Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
    spacecraft_params = transfer_dataset_inputs.spacecraft_params;
    
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
    
    %% Optimize Transfer that Leverages J2 using Genetic Algorithm 
    % Optimization variables
    eta1_bounds = [0, 0.9]; % [] initial -> intermediate transfer min efficiency
    eta2_bounds = [0, 0.9]; % [] intermediate -> target transfer min efficiency
    a_int_bounds = ([600, 2000] + R_E) / char_star.l; % [km] 
    %e_int_bounds = [1e-5, 0.04]; % []
    i_int_bounds = [0.97, 1.03] * x0_c_keplerian(3); % [rad]
    % Omega_int - assume same as original - all adjustments done by J2
    % omega_int - assume same as original - target almost circular
    var_bounds = [eta1_bounds; eta2_bounds; a_int_bounds; i_int_bounds];
    
    % Set up MultiObj
    MultiObj.fun = @(x) QLaw_J2_drift_transfer(x0_d_keplerian, [x(3) * char_star.l; [1e-5; x(4)]; x0_d_keplerian(4:6)], x0_c_keplerian, x(1:2), mu_E, R_E, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, 1.5, 365.25 * 2);
    MultiObj.nVar = size(var_bounds, 1);
    MultiObj.var_min = var_bounds(:, 1)';
    MultiObj.var_max = var_bounds(:, 2)';
    MultiObj.obj_names = ["ToF [days]", "Delta V [km / s]"];
    
    %% Optimize Pareto front
    options = optimoptions('paretosearch', 'ParetoSetSize', 60, 'UseParallel',true, 'MaxTime', 1200,'Display','iter',...
        'PlotFcn',{'psplotparetof','psplotparetox'}); % Could use custom plotting function that shows orbits
    fun = MultiObj.fun;
    lb = MultiObj.var_min;
    ub = MultiObj.var_max;
    rng shuffle % For reproducibility
    min_r_p = 400; % [km] minimum allowable periapsis (for drag reasons)
    tic;
    [x,fval,exitflag,output] = paretosearch(fun,MultiObj.nVar,[],[],[],[],lb,ub,@(x) min_periapsis_constraint(x(:, 3) * R_E, 1e-5, min_r_p, R_E),options);
    opt_time = toc;

    fprintf("Completed %g-%g with %g in %.1f seconds\n", ID1, ID2, exitflag, opt_time)

    %% Save results
    pareto = struct();
    pareto.x = x;
    pareto.fval = fval;
    pareto.exitflag = exitflag;
    pareto.output = output;
    %save("Mission Design\Multi Debris Mission Optimization\Deorbit to Debris Paretos\Low Thrust\" + run_name, "pareto");

    %% Analyze Results
    % figure
    % scatter3(x(:, 3) * R_E - R_E, rad2deg(x(:, 5)), fval(:, 1));
    % grid on
    % xlabel("Semi Major Alt [km]")
    % ylabel("Inclination []")
    % zlabel("dV [km / s]")
    
    %%
    x_intermediate = [x(:, 3)' * char_star.l; [1e-5 * ones([1, size(x, 1)]); x(:, 4)']; repmat(x0_d_keplerian(4:6), 1, size(x, 1))];
    
    figure
    earthy(R_E, "Earth", 0.6, [0;0;0]); hold on;
    C = cool();
    colormap(C)
    colors = interp1(linspace(min(fval(:, 1)), max(fval(:, 1)), size(C, 1)), C, fval(:, 1));
    for i = 1 : size(x_intermediate, 2)
        plotOrbit3(x_intermediate(4, i), x_intermediate(3, i), x_intermediate(5, i), x_intermediate(1, i) .* (1 - x_intermediate(2, i) .^ 2), x_intermediate(2, i), linspace(0, 2 * pi, 200), colors(i, :), 1, 0, [0; 0; 0], 0, 1); hold on
    end
    plotOrbit3(x0_d_keplerian(4), x0_d_keplerian(3), x0_d_keplerian(5), x0_d_keplerian(1) .* (1 - x0_d_keplerian(2) .^ 2), x0_d_keplerian(2), linspace(0, 2 * pi, 200), "g", 1, 0, [0; 0; 0], 1, 1); hold on
    plotOrbit3(x0_d_keplerian(4), x0_c_keplerian(3), x0_c_keplerian(5), x0_c_keplerian(1) .* (1 - x0_c_keplerian(2) .^ 2), x0_c_keplerian(2), linspace(0, 2 * pi, 200), "r", 1, 0, [0; 0; 0], 1, 1); hold on
    plotOrbit3(x0_c_keplerian(4), x0_c_keplerian(3), x0_c_keplerian(5), x0_c_keplerian(1) .* (1 - x0_c_keplerian(2) .^ 2), x0_c_keplerian(2), linspace(0, 2 * pi, 200), "r", 1, 0, [0; 0; 0], 1, 1); hold on
    colorbar()
    clim([min(fval(:, 1)), max(fval(:, 1))])
    grid on
    xlabel("X [km]")    
    ylabel("Y [km]")
    zlabel("Z [km]")
    title(sprintf("Pareto Front Transfers From Debris %g Deorbit to Debris %g", ID1, ID2))
    axis equal
end


%%
[ToF_sorted, ToF_sorted_i] = sort(fval(:, 2));
for k = 1 : (numel(ToF_sorted) - 1) % Fix duplicate ToFs...
    if ToF_sorted(k) == ToF_sorted(k + 1)
        ToF_sorted(k + 1) = ToF_sorted(k + 1) + 1e-10;
    end
end

ToF = ToF_sorted / 365.25;
dV = fval(ToF_sorted_i, 1);

figure
plot(ToF, dV);
xlabel("ToF [year]")
ylabel("Delta V [km / s]")
grid on

%%
[dV_ToF, Qtransfer_to_int, Qtransfer_to_targ] = QLaw_J2_drift_transfer(x0_d_keplerian, [x(end, 3) * char_star.l; [1e-5; x(end, 4)]; x0_d_keplerian(4:6)], x0_c_keplerian, x(end, 1:2), mu_E, R_E, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, 1.5, 365.25 * 2, true, 10);



%%
figure
plot_orbit_transfer_histories(Qtransfer_to_int.t / 60 / 60, x0_c_keplerian' ./ [R_E, ones([1, 5])], Qtransfer_to_int.x_keplerian_mass(1:6, :)' ./ [R_E, ones([1, 5])], Qtransfer_to_int.u_cont');
figure
plot_orbit_transfer_histories(Qtransfer_to_targ.t / 60 / 60, x0_c_keplerian' ./ [R_E, ones([1, 5])], Qtransfer_to_targ.x_keplerian_mass(1:6, :)' ./ [R_E, ones([1, 5])], Qtransfer_to_targ.u_cont');

save("terminator_to_new_debris_Qtransfer.mat", "Qtransfer_to_targ", "Qtransfer_to_int")

%% Helper Functions
function [c, ceq] = min_periapsis_constraint(a, e, min_r_p, R_E)
    ceq = [];
    r_p = a .* (1 - e) - R_E;
    c = min_r_p - r_p;
end

function [dV_ToF, Qtransfer_to_int, Qtransfer_to_targ] = QLaw_J2_drift_transfer(x_keplerian_0, x_keplerian_int, x_keplerian_targ, eta, mu, R, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV, max_ToF, plot_results, angular_step_deg)
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
        plot_results = false
        angular_step_deg = 20
    end

    % Transfer to intermediate orbit
    Q_params.eta_a_min = eta(1); % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = eta(1); % Minimum relative efficiency for thrusting instead of coasting

    [Qtransfer_to_int] = QLaw_transfer_fast(x_keplerian_0, x_keplerian_int, mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 50000, angular_step=deg2rad(angular_step_deg), max_dV = max_dV, max_t = max_ToF * 60 * 60 * 24);
    transfer_drift_1 = sum(J2_RAAN_drift(Qtransfer_to_int.x_keplerian_mass(1, :), Qtransfer_to_int.x_keplerian_mass(2, :), Qtransfer_to_int.x_keplerian_mass(3, :), mu, R, J_2_val) .* [diff(Qtransfer_to_int.t)', 0]);

    % Transfer to target orbit
    Q_params.eta_a_min = eta(2); % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = eta(2); % Minimum relative efficiency for thrusting instead of coasting
    
    spacecraft_params.m_0 = spacecraft_params.m_0 - Qtransfer_to_int.delta_m;
        
    [Qtransfer_to_targ] = QLaw_transfer_fast(x_keplerian_int, [x_keplerian_targ(1:3); x_keplerian_0(4); x_keplerian_targ(5:6)], mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 50000, angular_step=deg2rad(angular_step_deg), max_dV = max_dV, max_t = max_ToF * 60 * 60 * 24);
    transfer_drift_2 = sum(J2_RAAN_drift(Qtransfer_to_targ.x_keplerian_mass(1, :), Qtransfer_to_targ.x_keplerian_mass(2, :), Qtransfer_to_targ.x_keplerian_mass(3, :), mu, R, J_2_val) .* [diff(Qtransfer_to_targ.t)', 0]);

    % Calculate wait time for RAAN phasing accounting for drift during transfers
    targ_Omega_transfer_drift = J2_RAAN_drift(x_keplerian_targ(1), x_keplerian_targ(2), x_keplerian_targ(3), mu, R, J_2_val) * (Qtransfer_to_int.dt + Qtransfer_to_targ.dt);
    delta_Omega = wrapTo2Pi((x_keplerian_targ(4) + targ_Omega_transfer_drift) - (x_keplerian_0(4) + transfer_drift_1 + transfer_drift_2));
    rel_Omega_drift = J2_RAAN_drift(x_keplerian_int(1, :), x_keplerian_int(2, :), x_keplerian_int(3, :), mu, R, J_2_val) ...
                    - J2_RAAN_drift(x_keplerian_targ(1, :), x_keplerian_targ(2, :), x_keplerian_targ(3, :), mu, R, J_2_val);
    t_wait = delta_Omega ./ rel_Omega_drift .* (rel_Omega_drift > 0) ...
           + (delta_Omega - 2 * pi) ./ rel_Omega_drift .* (rel_Omega_drift < 0);

    % Package outputs
    dV_total = Qtransfer_to_int.delta_V + Qtransfer_to_targ.delta_V;
    ToF_total = (Qtransfer_to_int.dt + t_wait + Qtransfer_to_targ.dt) / 60 / 60 / 24;

    % Constraints
    n_constraints = 4; % max dV, min dV, transfer 1 converge, transfer 2 converge
    violated_constraints = (dV_total > max_dV) + (ToF_total > max_ToF) + ~Qtransfer_to_int.converged + ~Qtransfer_to_targ.converged;
    dV_ToF = [dV_total .* (violated_constraints == 0) + 1e5 * (violated_constraints / n_constraints),... 
              ToF_total .* (violated_constraints == 0) + 1e5 * (violated_constraints / n_constraints)];

    if plot_results
        char_star = load_charecteristic_values_Earth();

        Omegachange_to_int = cumsum(J2_RAAN_drift(Qtransfer_to_int.x_keplerian_mass(1, :), Qtransfer_to_int.x_keplerian_mass(2, :), Qtransfer_to_int.x_keplerian_mass(3, :), mu, R, J_2_val) .* [diff(Qtransfer_to_int.t)', 0]);
        Qtransfer_to_int.x_keplerian_mass(4, :) = x_keplerian_0(4) + Omegachange_to_int - (2 * pi / (365.25 * 60 * 60 * 24) * Qtransfer_to_int.t');

        Omegachange_int = J2_RAAN_drift(x_keplerian_int(1), x_keplerian_int(2), x_keplerian_int(3), mu, R, J_2_val) * t_wait;

        Omegachange_to_targ = cumsum(J2_RAAN_drift(Qtransfer_to_targ.x_keplerian_mass(1, :), Qtransfer_to_targ.x_keplerian_mass(2, :), Qtransfer_to_targ.x_keplerian_mass(3, :), mu, R, J_2_val) .* [diff(Qtransfer_to_targ.t)', 0]);
        Qtransfer_to_targ.x_keplerian_mass(4, :) = x_keplerian_0(4) + Omegachange_to_int(end) + Omegachange_int + Omegachange_to_targ - (2 * pi / (365.25 * 60 * 60 * 24) * (Qtransfer_to_int.t(end) + t_wait + Qtransfer_to_targ.t'));

        figure
        x_keplerian_cartesian_d = keplerian_to_cartesian_array(Qtransfer_to_int.x_keplerian_mass(1:6, :), [], char_star.mu);
        not_coast_colors = interp1(1:numel(Qtransfer_to_int.not_coast), double(Qtransfer_to_int.not_coast), linspace(1, numel(Qtransfer_to_int.not_coast), numel(Qtransfer_to_int.t)), 'nearest');
        plot_cartesian_orbit_color_varying(x_keplerian_cartesian_d(:, 1:1:end), not_coast_colors, 3); hold on

        x_keplerian_cartesian_d = keplerian_to_cartesian_array(Qtransfer_to_targ.x_keplerian_mass(1:6, :), [], char_star.mu);
        not_coast_colors = interp1(1:numel(Qtransfer_to_targ.not_coast), double(Qtransfer_to_targ.not_coast), linspace(1, numel(Qtransfer_to_targ.not_coast), numel(Qtransfer_to_targ.t)), 'nearest');
        plot_cartesian_orbit_color_varying(x_keplerian_cartesian_d(:, 1:1:end), not_coast_colors, 3); hold on

        c = colorbar;
        plotOrbit3(x_keplerian_0(4), x_keplerian_0(3), x_keplerian_0(5), x_keplerian_0(1) * (1 - x_keplerian_0(2) ^2), x_keplerian_0(2), linspace(0, 2 * pi, 1000), "g", 1, 1, [0, 0, 0], 0.1, 2); hold on
        plotOrbit3(x_keplerian_0(4) + Omegachange_to_int(end) - (2 * pi / (365.25 * 60 * 60 * 24) * Qtransfer_to_int.t(end)'), x_keplerian_int(3), x_keplerian_int(5), x_keplerian_int(1) * (1 - x_keplerian_int(2) ^2), x_keplerian_int(2), linspace(0, 2 * pi, 1000), "b", 1, 1, [0, 0, 0], 0.1, 2)
        plotOrbit3(x_keplerian_0(4) + Omegachange_to_int(end) + Omegachange_int - (2 * pi / (365.25 * 60 * 60 * 24) * (Qtransfer_to_int.t(end) + t_wait)'), x_keplerian_int(3), x_keplerian_int(5), x_keplerian_int(1) * (1 - x_keplerian_int(2) ^2), x_keplerian_int(2), linspace(0, 2 * pi, 1000), "b", 1, 1, [0, 0, 0], 0.1, 2)
        plotOrbit3(x_keplerian_targ(4) + J2_RAAN_drift(x_keplerian_targ(1), x_keplerian_targ(2), x_keplerian_targ(3), mu, R, J_2_val) * (Qtransfer_to_int.t(end) + t_wait + Qtransfer_to_targ.t(end)) - (2 * pi / (365.25 * 60 * 60 * 24) * (Qtransfer_to_int.t(end) + t_wait + Qtransfer_to_targ.t(end))'), x_keplerian_targ(3), x_keplerian_targ(5), x_keplerian_targ(1) * (1 - x_keplerian_targ(2) ^2), x_keplerian_targ(2), linspace(0, 2 * pi, 1000), "k", 1, 1, [0, 0, 0], 0.1, 2); hold on
        grid on,
        earthy(char_star.l, "Earth", 0.5, [0;0;0]); hold on; 
        clim([0, 1])
        c.Ticks = [0, 0.5, 1];
        c.TickLabels = {'Coast', 'Eclipse', 'Thrust'};
        axis equal
        title("Q-Law Orbit Transfer")
        %legend("Spacecraft", "Target", "Initial")
        xlabel("X [km]")
        ylabel("Y [km]")
        zlabel("Z [km]")
    end
end


function [] = plot_orbit_transfer_histories(t_hr, x_c, x_d, u)
    delta_x = x_c(:, 1:5) - x_d(:, 1:5);
        
    nexttile
    plot(t_hr, delta_x(:, 1), t_hr, delta_x(:, 2), t_hr, delta_x(:, 3), t_hr, delta_x(:, 4), t_hr, delta_x(:, 5)); hold off 
    title("Nondimensionalized Orbital Element Error")
    xlabel("Time [hr]")
    ylabel("\delta x")
    legend("a", "e", "i", "\omega", "\Omega", Location="northeast")
    grid on
    
    nexttile
    plot(t_hr, u(:, 1), t_hr, u(:, 2), t_hr, u(:, 3)); hold on 
    plot(t_hr, vecnorm(u, 2, 2), LineStyle= "--"); hold off
    title("Control")
    xlabel("Time [hr]")
    ylabel("u [km / s2]")
    legend("u_1", "u_2", "u_3", "||u||", Location="northeast")
    grid on
end
