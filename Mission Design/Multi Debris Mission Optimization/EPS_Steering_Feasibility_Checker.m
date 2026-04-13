%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continuous EP Steering Feasibility Checker
% Analyzes tracking rate and torque limits for a continuous trajectory
% AAE 450 OD5
% Last modified:
% - Atharva
% - 4/11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;

% Calling in Travis's Generated Trajectory
%orbit_trajectory = load('terminator_to_new_debris_Qtransfer.mat')
load('deorbit_Qtransfers.mat');
orbit_trajectory = Qtransfer_deorbit(3);
%% Spacecraft Parameters
m_r = 4000; % kg
m_s = 2000; % kg
I_sc_val = 18600; % kg*m^2 (Worst-case Transverse Axis)
I_sc = diag([I_sc_val, I_sc_val, 2000]); % Simplified MOI Matrix
s = 2.5; % m
I_sc = diag([1/6*m_s*s^2, 1/6*m_s*s^2, 1/6*m_s*s^2]);

% Inertia of Rocket about its own CM
L = 12; % m 
R = 1.85;
J_r = diag([1/12*m_r*(3*R^2 + L^2), 1/12*m_r*(3*R^2 + L^2), 1/2*m_r*R^2]);

% Combined Center of Mass
z_s = -(L/2 + s/2);
z_cm = (m_r*0 + m_s*z_s) / (m_r + m_s);

% Combined MOI
d_r = abs(0 - z_cm);
d_s = abs(z_s - z_cm);
J_total = (J_r + m_r*diag([d_r^2, d_r^2, 0])) + ...
          (I_sc + m_s*diag([d_s^2, d_s^2, 0]));
I_sc = J_total;
% Reaction Wheel
RW.tau = 0.2; % Nm
RW.H = 12; %25;  % Nms
RW.name = 'Reaction Wheel';

% Control Moment Gyro
CMG.tau = 45.0; % Nm
CMG.H = 48.0; % Nms (System envelope)
CMG.name = 'Control Moment Gyro';

% Extracting trajectory parameters
u = orbit_trajectory.u_cont;

t = orbit_trajectory.t;
dt = gradient(t);

num_pts = length(t);
u_dir = zeros(3, num_pts);
for i = 1:num_pts
    if (norm(u(:, i)) ~= 0)
        u_dir(:, i) = u(:, i) / norm(u(:,i));
    end
end

%% Kinematics Stuff
u_dot = zeros(3, num_pts);
u_gradient = zeros(3, num_pts);
for row = 1:3
    u_gradient(row, :) = gradient(u_dir(row, :));
end
for i = 1:num_pts
    u_dot(:, i) = u_gradient(:, i) ./ dt(i); % Gradient calcs rate between entries
end



% w = u x u_dot 
w_req = zeros(3, num_pts);

for i = 1:num_pts
    w_req(:, i) = cross(u_dir(:,i), u_dot(:,i), 1);
end

% Calculate the required angular acceleration (alpha = dw/dt)
alpha_req = zeros(3, num_pts);
w_req_gradient = zeros(3, num_pts);
for row = 1:3
    w_req_gradient(row, :) = gradient(w_req(row, :));
end
for i = 1:num_pts
    alpha_req(:, i) = w_req_gradient(:, i) / dt(i);
end

%% Calculate Required Continuous Torque
tau_req = zeros(3, num_pts);
for i = 1:num_pts
    w = w_req(:, i);
    alpha = alpha_req(:, i);
    
    % Tau = I*alpha + w x (I*w)
    gyroscopic_term = cross(w, I_sc * w);
    tau_req(:, i) = (I_sc * alpha) + gyroscopic_term;
end

% Extract magnitude limits for the whole trajectory
max_w_mag = max(vecnorm(w_req, 2, 1)); % rad/s
[max_tau_mag, max_tau_ind] = max(vecnorm(tau_req, 2, 1)); % Nm
max_tau_mag = max(tau_req, [], [1 2]);
max_w_mag = max(w_req, [], [1 2]);

fprintf('--- EP Continuous Steering Profile Analysis ---\n');
fprintf('Duration: %.1f minutes\n', max(t)/60);
fprintf('Peak Tracking Rate: %.4f deg/s\n', rad2deg(max_w_mag));
fprintf('Peak Dynamic Torque: %.4f Nm\n\n', max_tau_mag);

%% Evaluate Hardware
evaluate_tracking(RW, I_sc_val, max_w_mag, max_tau_mag);
evaluate_tracking(CMG, I_sc_val, max_w_mag, max_tau_mag);

function evaluate_tracking(sys, I, max_w, max_tau)
    % Check if the maneuver requires a slew rate faster than the H_max limit
    w_limit = sys.H / I;
    
    if max_w > w_limit
        w_status = 'FAIL';
    else
        w_status = 'PASS';
    end
    
    % Check if the maneuver requires instantaneous torque > tau_max
    if max_tau > sys.tau
        tau_status = 'FAIL';
    else
        tau_status = 'PASS';
    end
    
    fprintf('[%s] %s\n', sys.name, repmat('-', 1, 30));
    fprintf('Rate Check:   [%s] (Req: %.4f deg/s | Limit: %.4f deg/s)\n', ...
            w_status, rad2deg(max_w), rad2deg(w_limit));
    fprintf('Torque Check: [%s] (Req: %.4f Nm | Limit: %.4f Nm)\n\n', ...
            tau_status, max_tau, sys.tau);
end

%%% Plotting

figure
subplot(2,1,1)
plot(t, rad2deg(w_req(1,:)), "r")
hold on
plot(t, rad2deg(w_req(2,:)), "g")
plot(t, rad2deg(w_req(3,:)), "b")
grid on
xlabel("Time (s)")
ylabel("w_req")
legend(["w(1)", "w(2)", "w(3)"])

subplot(2,1,2)
plot(t, alpha_req(1,:), "r")
hold on
plot(t, alpha_req(2,:), "g")
plot(t, alpha_req(3,:), "b")
grid on
xlabel("Time (s)")
ylabel("alpha_req")
legend(["w(1)", "w(2)", "w(3)"])

figure
subplot(2,1,1)
plot(t, vecnorm(tau_req, 2, 1), "r")
grid on
xlabel("Time (s)")
ylabel("Torque Required (Nm)")

subplot(2, 1, 2)
plot(t, tau_req)