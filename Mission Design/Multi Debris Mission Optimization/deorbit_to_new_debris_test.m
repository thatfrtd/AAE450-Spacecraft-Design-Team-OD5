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

R_E = 6378.137; % [km] Earth radius
mu_E = 398600.4418; % [km3 / s2] Earth gravitational parameter
J_2_val = 1.08262668e-3; % [] Earth J2

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = 7044.7634; % [km] semi-major axis
e_c = 0.003390; % [] eccentricity
i_c = deg2rad(98.1114 ); % [rad] inclination
Omega_c = deg2rad(320.5520 - 20 ); % [rad] right ascension of ascending node
omega_c = deg2rad(301.2069 - 0 ); % [rad] argument of periapsis
nu_c = deg2rad(58.6658 ); % [rad] true anomaly at epoch

M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);

% Initial conditions for spacecraft
r_a_0 = R_E + 600; % [km] periapsis
r_p_0 = R_E + 135; % [km] periapsis
e_d = (1 - r_p_0 / r_a_0) / (1 + r_p_0 / r_a_0); % [] eccentricity
a_d = r_p_0 / (1 - e_d); % [km] semi-major axis
i_d = deg2rad(98.1114); % [rad] inclination
Omega_d = deg2rad(320.5520); % [rad] right ascension of ascending node
omega_d = deg2rad(301.2069); % [rad] argument of periapsis
nu_d = deg2rad(58.6658 ); % [rad] true anomaly at epoch

M_d = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_d, e_d), e_d);
x0_d_keplerian = [a_d; e_d; i_d; Omega_d; omega_d; M_d];
x0_d_cartesian = keplerian_to_cartesian(x0_d_keplerian, nu_d, mu_E);

char_star = load_charecteristic_values_Earth();

% Spacecraft Parameters: Isp, max thrust, initial mass, fuel mass
spacecraft_params = struct();
spacecraft_params.Isp = 4100; % [s]
spacecraft_params.m_0 = 1500; % [kg]
spacecraft_params.m_dry = 600; % [kg]
spacecraft_params.F_max = 0.25; % [N]

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
eta1_bounds = [0, 0.85]; % [] initial -> intermediate transfer min efficiency
eta2_bounds = [0, 0.85]; % [] intermediate -> target transfer min efficiency
a_int_bounds = ([600, 1400] + R_E) / char_star.l; % [km] 
e_int_bounds = [1e-5, 0.04]; % []
i_int_bounds = [0.97, 1.03] * i_c; % [rad]
% Omega_int - assume same as original - all adjustments done by J2
% omega_int - assume same as original - target almost circular
var_bounds = [eta1_bounds; eta2_bounds; a_int_bounds; e_int_bounds; i_int_bounds];

% Set up MultiObj
MultiObj.fun = @(x) QLaw_J2_drift_transfer(x0_d_keplerian, [x(3) * char_star.l; x(4:5)'; x0_d_keplerian(4:6)], x0_c_keplerian, x(1:2), mu_E, R_E, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, 2, 365 * 2);
MultiObj.nVar = size(var_bounds, 1);
MultiObj.var_min = var_bounds(:, 1)';
MultiObj.var_max = var_bounds(:, 2)';
MultiObj.obj_names = ["ToF [days]", "Delta V [km / s]"];

%% Test
x_test = var_bounds(:, 1)';
MultiObj.fun(x_test)

%% Create Pool
p = gcp("nocreate"); % If no pool, do not create new one.
if isempty(p)
    p = parpool(8);
end

%% Optimize Pareto front
options = optimoptions('paretosearch', 'ParetoSetSize', 60, 'UseParallel',true, 'MaxTime', 1200,'Display','iter',...
    'PlotFcn',{'psplotparetof','psplotparetox'}); % Could use custom plotting function that shows orbits
fun = MultiObj.fun;
lb = MultiObj.var_min;
ub = MultiObj.var_max;
rng shuffle % For reproducibility
min_r_p = 400; % [km] minimum allowable periapsis (for drag reasons)
[x,fval,exitflag,output] = paretosearch(fun,MultiObj.nVar,[],[],[],[],lb,ub,@(x) min_periapsis_constraint(x(:, 3) * R_E, x(:, 4), min_r_p, R_E),options);

%% Analyze Results
figure
scatter3(x(:, 3) * R_E - R_E, rad2deg(x(:, 5)), fval(:, 1));
grid on
xlabel("Semi Major Alt [km]")
ylabel("Inclination []")
zlabel("dV [km / s]")

%%
x_intermediate = [x(:, 3)' * char_star.l; x(:, 4:5)'; repmat(x0_d_keplerian(4:6), 1, size(x, 1))];

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
axis equal

%% Helper Functions
function [c, ceq] = min_periapsis_constraint(a, e, min_r_p, R_E)
    ceq = [];
    r_p = a .* (1 - e) - R_E;
    c = min_r_p - r_p;
end

function [dV_ToF] = QLaw_J2_drift_transfer(x_keplerian_0, x_keplerian_int, x_keplerian_targ, eta, mu, R, J_2_val, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, max_dV, max_ToF)
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

    [Qtransfer_to_int] = QLaw_transfer_fast(x_keplerian_0, x_keplerian_int, mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 50000, angular_step=deg2rad(20));
    transfer_drift_1 = sum(J2_RAAN_drift(Qtransfer_to_int.x_keplerian_mass(1, :), Qtransfer_to_int.x_keplerian_mass(2, :), Qtransfer_to_int.x_keplerian_mass(3, :), mu, R, J_2_val) .* [diff(Qtransfer_to_int.t)', 0]);

    % Transfer to target orbit
    Q_params.eta_a_min = eta(2); % Minimum absolute efficiency for thrusting instead of coasting
    Q_params.eta_r_min = eta(2); % Minimum relative efficiency for thrusting instead of coasting
    
    spacecraft_params.m_0 = spacecraft_params.m_0 - Qtransfer_to_int.delta_m;
        
    [Qtransfer_to_targ] = QLaw_transfer_fast(x_keplerian_int, [x_keplerian_targ(1:3); x_keplerian_0(4); x_keplerian_targ(5:6)], mu, spacecraft_params, Q_params, penalty_params, Qdot_opt_params, return_dt_dm_only = false, iter_max = 50000, angular_step=deg2rad(20));
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
    ToF_total = (Qtransfer_to_int.dt + t_wait + Qtransfer_to_targ.dt) / 60 / 60 / 24;

    % Constraints
    n_constraints = 4; % max dV, min dV, transfer 1 converge, transfer 2 converge
    violated_constraints = (dV_total > max_dV) + (ToF_total > max_ToF) + ~Qtransfer_to_int.converged + ~Qtransfer_to_targ.converged;
    dV_ToF = [dV_total .* (violated_constraints == 0) + 1e5 * (violated_constraints / n_constraints),... 
              ToF_total .* (violated_constraints == 0) + 1e5 * (violated_constraints / n_constraints)];
end
