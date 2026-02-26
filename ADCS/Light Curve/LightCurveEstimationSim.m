% =========================================================================
% Light Curve Attitude Estimation via Particle Swarm Optimization (PSO)
% Last Modified
%   - 2/22
%   - Atharva
% =========================================================================
clear; clc; close all;

%% Simulation Setup & "truth" Light Curve Generation
% Time parameters
dt = 1; % 1 second measurement cadence
t_end = 100; % 100 seconds total
t_span = 0:dt:t_end;

% Satellite Properties (Simplified as a single flat plate for demonstration)
I_tensor = diag([10, 15, 20]); % Moment of Inertia (kg*m^2)
plate_area = 2.0; % m^2
plate_normal_body = [0; 0; 1]; % Normal vector pointing out of Z-axis

% Sun and Observer Geometry (Assumed fixed in inertial frame for short pass)
sun_vec_inertial = normalize([1; 1; 0], 'norm');
obs_vec_inertial = normalize([1; 0; 1], 'norm');

% True Initial State: [q1, q2, q3, q4, wx, wy, wz]
% q4 is the scalar part.
q0_true = [0; 0; 0; 1]; 
w0_true = [0.1; 0.2; 0.05]; % Tumbling angular velocity (rad/s)
true_initial_state = [q0_true; w0_true];

% Generate True Light Curve
[~, state_true] = ode45(@(t, y) rigid_body_dynamics(t, y, I_tensor), t_span, true_initial_state);
true_light_curve = zeros(length(t_span), 1);

for k = 1:length(t_span)
    q_k = state_true(k, 1:4)';
    % Convert quaternion to rotation matrix (Body to Inertial)
    R_B2I = quat2rotm_custom(q_k); 
    
    % Rotate plate normal into inertial frame
    n_inertial = R_B2I * plate_normal_body;
    
    % Basic Lambertian diffuse reflection model
    % Intensity ~ max(0, N dot L) * max(0, N dot V)
    sun_term = max(0, dot(n_inertial, sun_vec_inertial));
    obs_term = max(0, dot(n_inertial, obs_vec_inertial));
    
    true_light_curve(k) = plate_area * sun_term * obs_term;
end

% Add some measurement noise (simulating CCD sensor noise)
measured_light_curve = true_light_curve + 0.05 * randn(size(true_light_curve));
measured_light_curve = max(0, measured_light_curve); % Intensity can't be negative

%% Attitude Estimation using Particle Swarm
disp('Starting Particle Swarm Optimization for Attitude Estimation...');

% State to estimate: 3 components of Rodrigues Parameters (or Euler axis/angle) 
% + 3 angular velocities. (Using 6 params avoids 4D quaternion constraint issues in PSO).
% Let's optimize: [wx, wy, wz, e1, e2, e3] (where e is Euler vector)

num_vars = 6;
lower_bounds = [-0.5, -0.5, -0.5, -pi, -pi, -pi];
upper_bounds = [ 0.5,  0.5,  0.5,  pi,  pi,  pi];

% PSO Options (Increase SwarmSize for better global search, takes longer)
options = optimoptions('particleswarm', 'SwarmSize', 100, 'Display', 'iter', 'MaxIterations', 50);

% Define anonymous objective function
cost_func = @(x) light_curve_cost(x, t_span, measured_light_curve, I_tensor, plate_normal_body, plate_area, sun_vec_inertial, obs_vec_inertial);

% Run PSO
[best_estimate, fval] = particleswarm(cost_func, num_vars, lower_bounds, upper_bounds, options);

disp('Optimization Complete.');
disp('True Initial Angular Velocity (rad/s):'); 
disp(w0_true');
disp('Estimated Initial Angular Velocity (rad/s):'); 
disp(best_estimate(1:3));

%% 3. Plotting the Results
% Re-simulate with the estimated parameters to see how well the curves match
w_est = best_estimate(1:3)';
e_est = best_estimate(4:6)';
angle_est = norm(e_est);
if angle_est == 0
    q_est = [0;0;0;1];
else
    q_est = [e_est/angle_est * sin(angle_est/2); cos(angle_est/2)];
end
est_initial_state = [q_est; w_est];

[~, state_est] = ode45(@(t, y) rigid_body_dynamics(t, y, I_tensor), t_span, est_initial_state);
est_light_curve = zeros(length(t_span), 1);
for k = 1:length(t_span)
    R_B2I = quat2rotm_custom(state_est(k, 1:4)'); 
    n_inertial = R_B2I * plate_normal_body;
    est_light_curve(k) = plate_area * max(0, dot(n_inertial, sun_vec_inertial)) * max(0, dot(n_inertial, obs_vec_inertial));
end

figure;
plot(t_span, measured_light_curve, 'k.', 'DisplayName', 'Noisy Measurements'); hold on;
plot(t_span, true_light_curve, 'b-', 'LineWidth', 2, 'DisplayName', 'True Light Curve');
plot(t_span, est_light_curve, 'r--', 'LineWidth', 2, 'DisplayName', 'Estimated Light Curve');
xlabel('Time (s)'); ylabel('Brightness / Intensity');
title('Light Curve Inversion via PSO'); legend; grid on;

% =========================================================================
% Helper Functions
% =========================================================================

function J = light_curve_cost(x, t_span, measured_lc, I_tensor, n_body, area, sun_vec, obs_vec)
    % Extract angular velocity and attitude vector from PSO particle
    w0 = x(1:3)';
    e0 = x(4:6)';
    
    % Convert Euler vector to quaternion
    angle = norm(e0);
    if angle == 0
        q0 = [0;0;0;1];
    else
        q0 = [e0/angle * sin(angle/2); cos(angle/2)];
    end
    
    initial_state = [q0; w0];
    
    % Propagate dynamics
    [~, state] = ode45(@(t, y) rigid_body_dynamics(t, y, I_tensor), t_span, initial_state);
    
    % Generate Simulated Light Curve
    sim_lc = zeros(length(t_span), 1);
    for k = 1:length(t_span)
        R = quat2rotm_custom(state(k, 1:4)'); 
        n_inertial = R * n_body;
        sim_lc(k) = area * max(0, dot(n_inertial, sun_vec)) * max(0, dot(n_inertial, obs_vec));
    end
    
    % Cost is the Sum of Squared Errors (SSE)
    J = sum((measured_lc - sim_lc).^2);
end

function dy = rigid_body_dynamics(~, y, I)
    q = y(1:4); % Quaternion [qx, qy, qz, qw] (qw is scalar)
    w = y(5:7); % Angular velocity
    
    % Normalize quaternion to prevent numerical drift
    q = q / norm(q); 
    
    % Quaternion kinematics matrix (using Hamilton convention)
    % dq/dt = 0.5 * Omega(w) * q
    wx = w(1); wy = w(2); wz = w(3);
    Omega = [  0,  wz, -wy,  wx;
             -wz,   0,  wx,  wy;
              wy, -wx,   0,  wz;
             -wx, -wy, -wz,   0 ];
         
    dq = 0.5 * Omega * q;
    
    % Euler's equations for torque-free motion
    % I * dw/dt + w x (I * w) = 0
    dw = I \ (-cross(w, I * w));
    
    dy = [dq; dw];
end

function R = quat2rotm_custom(q)
    % Converts a [qx; qy; qz; qw] quaternion to a rotation matrix
    q = q / norm(q);
    qx = q(1); qy = q(2); qz = q(3); qw = q(4);
    
    R = [1 - 2*(qy^2 + qz^2),     2*(qx*qy - qz*qw),     2*(qx*qz + qy*qw);
             2*(qx*qy + qz*qw), 1 - 2*(qx^2 + qz^2),     2*(qy*qz - qx*qw);
             2*(qx*qz - qy*qw),     2*(qy*qz + qx*qw), 1 - 2*(qx^2 + qy^2)];
end