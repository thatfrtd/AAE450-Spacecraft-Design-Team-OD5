%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% SIM parameters
tspan = [0 200]; % Simulation time [s]

% Common Control Input (u) & Parameters (p)
thrust = [0.05; 0; 0]; % [kN] Constant thrust for mass depletion
p = [];                % Empty parameter vector

% Constant: Isp
Isp = 300; % [s]

%% CHW model constants
r_c = 6778;         % [km] Radius of circular chief orbit
c_cwh = [r_c; Isp]; 

%% Linearized model constants
a_c     = 6778; % [km] Semi-major axis (match CWH r_c for comparison)
e_c     = 0.01; % Eccentricity (Set slightly > 0 to see differences from CWH)
i_c     = deg2rad(51.6); % [rad] Inclination
Omega_c = 0;    % [rad] RAAN
omega_c = 0;    % [rad] Argument of perigee
M0_c    = 0;    % [rad] Initial Mean Anomaly

x_keplerian_c = [a_c; e_c; i_c; Omega_c; omega_c; M0_c];
c_lin = [x_keplerian_c; Isp];

%% Define initial state
% x = [r (3x1); v (3x1); m (1x1)]
r0 = [0.1; 0.05; -0.02];   % [km]
v0 = [0.001; 0; 0.0005];   % [km/s]
m0 = 1000;                 % [kg]
x0 = [r0; v0; m0];

%% Propagate dynamics
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);

% Propagate CWH
ode_cwh = @(t, x) CWH_relative_orbit_EoM(t, x, thrust, p, c_cwh);
[t_cwh, x_cwh] = ode45(ode_cwh, tspan, x0, options);

% Propagate Linearized (Eccentric) Model
ode_lin = @(t, x) linearized_relative_orbit_EoM(t, x, thrust, p, c_lin);
[t_lin, x_lin] = ode45(ode_lin, tspan, x0, options);

%% Plotting
figure('Name', 'Orbit Model Comparison', 'Color', 'w', 'Position', [100 100 1000 600]);

% 3d comparison
subplot(2, 2, [1 3]);
plot3(x_cwh(:,1), x_cwh(:,2), x_cwh(:,3), 'b-', 'LineWidth', 2); 
hold on;
plot3(x_lin(:,1), x_lin(:,2), x_lin(:,3), 'r--', 'LineWidth', 2);
grid on;
axis equal;
xlabel('Radial (x) [km]'); 
ylabel('Along-Track (y) [km]'); 
zlabel('Cross-Track (z) [km]');
title('3D Relative Motion Trajectory');
legend('CWH (Circular)', 'Linearized (Eccentric)', 'Location', 'best');

% in plane motion
subplot(2, 2, 2);
plot(x_cwh(:,2), x_cwh(:,1), 'b-', 'LineWidth', 1.5); 
hold on;
plot(x_lin(:,2), x_lin(:,1), 'r--', 'LineWidth', 1.5);
grid on; axis equal;
xlabel('Along-Track (y) [km]'); 
ylabel('Radial (x) [km]');
title('In-Plane Motion (Top Down)');

% Model deviation over time
deviation = vecnorm(x_cwh(1:(length(x_cwh)+1)/2, 1:3) - interp1(t_lin, x_lin(:, 1:3), t_cwh(1:(length(x_cwh)+1)/2)), 2, 2);
subplot(2, 2, 4);
plot(t_cwh(1:(length(x_cwh)+1)/2), deviation, 'k', 'LineWidth', 1.5);
grid on;
xlabel('Time [s]'); 
ylabel('Position Difference [km]');
title(sprintf('Model Drift (e_c = %.3f)', e_c));

