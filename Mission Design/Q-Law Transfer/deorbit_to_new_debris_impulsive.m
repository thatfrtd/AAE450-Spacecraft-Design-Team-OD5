%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 450 Team OD5
% Transfer from a deorbit orbit (after deorbiting debris) to new debris
% Author: Travis Hastreiter 
% Created On: 13 March, 2026
% Description: Orbit transfer using Q-Law from deorbit orbit (after drop 
% off) to new debris not accounting for rendezvous (assuming not much extra 
% delta V and time).
% 
% Right now it uses impulsive Lambert transfers to make sure the procedure
% works.
% Most Recent Change: 14 March, 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_E = 6378.137; % [km] Earth radius
mu_E = 398600.4418; % [km3 / s2] Earth gravitational parameter
J_2_val = 1.08262668e-3; % [] Earth J2

% Initial conditions for target Earth orbit (in Earth Centered Inertial (ECI) frame)
a_c = 7044.7634; % [km] semi-major axis
e_c = 0.003390; % [] eccentricity
i_c = deg2rad(98.1114 ); % [rad] inclination
Omega_c = deg2rad(320.5520 + 280 ); % [rad] right ascension of ascending node
omega_c = deg2rad(301.2069 - 0 ); % [rad] argument of periapsis
nu_c = deg2rad(58.6658 ); % [rad] true anomaly at epoch

M_c = eccentric_to_mean_anomaly(true_to_eccentric_anomaly(nu_c, e_c), e_c);
x0_c_keplerian = [a_c; e_c; i_c; Omega_c; omega_c; M_c];
x0_c_cartesian = keplerian_to_cartesian(x0_c_keplerian, nu_c, mu_E);

% Initial conditions for spacecraft
r_a_0 = R_E + 600; % [km] periapsis
r_p_0 = R_E + 120; % [km] periapsis
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
spacecraft_params.m_dry = 1000; % [kg]
spacecraft_params.F_max = 0.25; % [N]

% Integration error tolerance
default_tolerance = 1e-10;

%% Optimize Transfer that Leverages J2 using Genetic Algorithm 
% Optimization variables
ToF1_bounds = [0.01, 2]; % [hr] initial -> intermediate transfer time of flight
ToF2_bounds = [0.01, 2]; % [hr] intermediate to target transfer time of flight
a_int_bounds = ([300, 2000] + R_E) / char_star.l; % [km] 
e_int_bounds = [0, 0.2]; % []
i_int_bounds = [0.9, 1.1] * i_c; % [rad]
% Omega_int - assume same as original - all adjustments done by J2
% omega_int - assume same as original - target almost circular
var_bounds = [ToF1_bounds; ToF2_bounds; a_int_bounds; e_int_bounds; i_int_bounds];

% Set up MultiObj
MultiObj.fun = @(x) batch_Lambert_J2_drift_transfer(x0_d_keplerian, [x(:, 3)' * char_star.l; x(:, 4:5)'; repmat(x0_d_keplerian(4:6), 1, size(x, 1))], x0_c_keplerian, x(:, 1:2)' * 3600, 20, mu_E, R_E, J_2_val, char_star, 2, 365);
MultiObj.nVar = size(var_bounds, 1);
MultiObj.var_min = var_bounds(:, 1)';
MultiObj.var_max = var_bounds(:, 2)';
MultiObj.obj_names = ["ToF [days]", "Delta V [km / s]"];

load_lambert();

%% Test
x_test = [var_bounds(:, 1)'; var_bounds(:, 2)'];
MultiObj.fun(x_test)

%%
options = optimoptions('paretosearch', 'ParetoSetSize', 100,'Display','iter',...
    'PlotFcn',{'psplotparetof' 'psplotparetox'});
fun = MultiObj.fun;
lb = MultiObj.var_min;
ub = MultiObj.var_max;
rng shuffle % For reproducibility
min_r_p = 400; % [km] minimum allowable periapsis (for drag reasons)
[x,fval,exitflag,output] = paretosearch(fun,MultiObj.nVar,[],[],[],[],lb,ub,@(x) min_periapsis_constraint(x(:, 3) * R_E, x(:, 4), min_r_p, R_E),options);

%%
unload_lambert();

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
function [dV_ToF] = batch_Lambert_J2_drift_transfer(x_keplerian_0, x_keplerian_int, x_keplerian_targ, ToFs, N_thetastar, mu, R, J_2_val, char_star, max_dV, max_ToF)
    
    N_pop = size(ToFs, 2);

    % Transfer # (2) * Direction (2) * N_thetastar * N_thetastar * N_population
    transfer_num = 1 : 2;
    directions = [1; -1];
    thetastars_1 = linspace(0, 2 * pi, N_thetastar);
    thetastars_2 = linspace(0, 2 * pi, N_thetastar);
    pop_ID = 1 : N_pop;
    transfer_combos = combinations(transfer_num, directions, thetastars_1, thetastars_2, pop_ID);

    % Compute all transfer options
    ToF_array = zeros([size(transfer_combos, 1), 1]);
    for i = 1 : numel(ToF_array)
        ToF_array(i) = ToFs(transfer_combos.transfer_num(i), transfer_combos.pop_ID(i));
    end
    x_1_array = [keplerian_to_cartesian_array(repmat(x_keplerian_0, 1, size(transfer_combos, 1) / 2), transfer_combos.thetastars_1(transfer_combos.transfer_num == 1), mu), ...
                 keplerian_to_cartesian_array(x_keplerian_int(:, transfer_combos.pop_ID(transfer_combos.transfer_num == 2)), transfer_combos.thetastars_1(transfer_combos.transfer_num == 2), mu)];
    x_2_array = [keplerian_to_cartesian_array(x_keplerian_int(:, transfer_combos.pop_ID(transfer_combos.transfer_num == 1)), transfer_combos.thetastars_2(transfer_combos.transfer_num == 1), mu), ...
                 keplerian_to_cartesian_array(repmat(x_keplerian_targ, 1, size(transfer_combos, 1) / 2), transfer_combos.thetastars_2(transfer_combos.transfer_num == 2), mu)];

    nd_scalar = [ones([3, 1]) * char_star.l; ones([3, 1]) * char_star.v];
    [vel1_raw, vel2_raw, dV_raw] = best_lambert_zeroN(x_1_array ./ nd_scalar, x_2_array ./ nd_scalar, ToF_array / char_star.t, 0, 0, direction = transfer_combos.directions);
    
    % Extract the best transfer info for each individual's two transfers
    dV_reshaped = reshape(dV_raw, [2, 2, N_thetastar, N_thetastar, N_pop]) * char_star.v;

    [dV1, dV1_i] = min(dV_reshaped(1, :, :, :, :), [], 1:4, "linear");
    [dV2, dV2_i] = min(dV_reshaped(2, :, :, :, :), [], 1:4, "linear");

    % Calculate wait time for RAAN phasing
    delta_Omega = wrapTo2Pi(x_keplerian_targ(4) - x_keplerian_0(4));
    rel_Omega_drift = J2_RAAN_drift(x_keplerian_int(1, :), x_keplerian_int(2, :), x_keplerian_int(3, :), mu, R, J_2_val) ...
                    - J2_RAAN_drift(x_keplerian_targ(1, :), x_keplerian_targ(2, :), x_keplerian_targ(3, :), mu, R, J_2_val);
    t_wait = delta_Omega ./ rel_Omega_drift .* (rel_Omega_drift > 0) ...
           + (delta_Omega - 2 * pi) ./ rel_Omega_drift .* (rel_Omega_drift < 0);

    % Package outputs
    dV_total = squeeze(dV1 + dV2);
    ToF_total = (ToFs(1, :)' + t_wait' + ToFs(2, :)') / 60 / 60 / 24;
    violated_constraints = (dV_total > max_dV) + (ToF_total > max_ToF);
    dV_ToF = [dV_total .* (violated_constraints == 0) + 1e5 * (violated_constraints / 2),... 
              ToF_total .* (violated_constraints == 0) + 1e5 * (violated_constraints / 2)];
end

function [c, ceq] = min_periapsis_constraint(a, e, min_r_p, R_E)
    ceq = [];
    r_p = a .* (1 - e) - R_E;
    c = min_r_p - r_p;
end