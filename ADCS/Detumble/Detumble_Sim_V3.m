%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Detumble Sim Analysis V3 
% This time it has some acc control law
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Init Variables
m_r = 4000; % kg
m_s = 1500; % kg
R = 1.85;   % m
L = 30;     % m
s = 1.0;    % m

% Inertia of Rocket about its own CM
J_r = diag([1/12*m_r*(3*R^2 + L^2), 1/12*m_r*(3*R^2 + L^2), 1/2*m_r*R^2]);
% Inertia of Satellite about its own CM
J_s = diag([1/6*m_s*s^2, 1/6*m_s*s^2, 1/6*m_s*s^2]);

% Combined Center of Mass (relative to rocket center)
z_s = -(L/2 + s/2);
z_cm = (m_r*0 + m_s*z_s) / (m_r + m_s);

% Parallel Axis Theorem to find J_total about combined CM
d_r = abs(0 - z_cm);
d_s = abs(z_s - z_cm);
J_total = (J_r + m_r*diag([d_r^2, d_r^2, 0])) + ...
          (J_s + m_s*diag([d_s^2, d_s^2, 0]));

%% Initial Conditions 
w0 = [0.005; -0.005; (2 * 2 * pi / 60)]; % rad/s (X, Y, and Z axes)

% Start at a random off-nominal orientation
q0 = [0.4; 0.3; -0.5; 0.7]; 
q0 = q0 / norm(q0); % Must be normalized!

x0 = [q0; w0];

%% Controller Parameters
tau_max = 0.4; % Nm (Max torque output of your reaction wheels)
Kp = 1.5;      % Proportional gain (Orientation error)
Kd = 50;       % Derivative gain (Rate damping)

%% ODE45
tspan = [0 12000]; 
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t, x] = ode45(@(t, x) plant_derivative(t, x, J_total, Kp, Kd, tau_max), tspan, x0, opts);

%% Extract states and Recalculate Control Torques for Plotting
q = x(:, 1:4);
w = x(:, 5:7);

% Pre-allocate torque array
tau_history = zeros(length(t), 3);

% Recalculate the torque at each time step to plot it
for i = 1:length(t)
    q_curr = q(i, :)';
    w_curr = w(i, :)';
    
    % PD Control Law (Aiming for q_target = [0;0;0;1] and w_target = [0;0;0])
    tau_cmd = -Kp * q_curr(1:3) - Kd * w_curr; 
    
    % Saturate torque at +/- 0.4 Nm
    tau_history(i, :) = max(min(tau_cmd, tau_max), -tau_max)';
end

%% Plotting Proof of Control Authority
figure('Color', 'w', 'Position', [100, 100, 800, 900]);

% 1. Angular Rates Plot (Proof of Detumble)
subplot(3,1,1);
plot(t, w, 'LineWidth', 1.5);
grid on; title('Angular Rates (\omega) - Detumbling to Zero');
ylabel('rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
yline(0, 'k--');

% 2. Quaternion Plot (Proof of Reorientation)
subplot(3,1,2);
plot(t, q, 'LineWidth', 1.5);
grid on; title('Quaternion Orientation - Reorienting to [0,0,0,1]');
ylabel('Components'); legend('\epsilon_1', '\epsilon_2', '\epsilon_3', '\eta');
yline(0, 'k--'); yline(1, 'k--');

% 3. Control Torque Plot (Proof of Actuator Effort)
subplot(3,1,3);
plot(t, tau_history, 'LineWidth', 1.5);
grid on; title('Reaction Wheel Control Torques (\tau)');
ylabel('Torque (Nm)'); legend('\tau_x', '\tau_y', '\tau_z');
xlabel('Time (s)');
ylim([-0.5 0.5]); % Show the 0.4 Nm physical saturation limit

figure
plot(t, w, 'LineWidth', 1.5);
grid on; title('Angular Rates (\omega) - Detumbling to Zero');
ylabel('rad/s'); legend('\omega_x', '\omega_y', '\omega_z');
yline(0, 'k--');

% 2. Quaternion Plot (Proof of Reorientation)
figure
plot(t, q, 'LineWidth', 1.5);
grid on; title('Quaternion Orientation - Reorienting to [0,0,0,1]');
ylabel('Components'); legend('\epsilon_1', '\epsilon_2', '\epsilon_3', '\eta');
yline(0, 'k--'); yline(1, 'k--');

% 3. Control Torque Plot (Proof of Actuator Effort)
figure
plot(t, tau_history, 'LineWidth', 1.5);
grid on; title('Reaction Wheel Control Torques (\tau)');
ylabel('Torque (Nm)'); legend('\tau_x', '\tau_y', '\tau_z');
xlabel('Time (s)');
ylim([-0.5 0.5]); % Show the 0.4 Nm physical saturation limit

%% Helper Function (Combined Kinematics and Dynamics)
function dxdt = plant_derivative(~, x, J, Kp, Kd, tau_max)
    % Extract states
    q = x(1:4);
    w = x(5:7);
    
    % Re-normalize quaternion
    q = q / norm(q); 
    
    % 1. Controller (PD Law)
    % Aiming for q_target = [0;0;0;1]. The error is essentially the vector part.
    tau_cmd = -Kp * q(1:3) - Kd * w;
    
    % Saturate torque to physical wheel limits
    tau = max(min(tau_cmd, tau_max), -tau_max);
    
    % 2. Kinematics (Quaternion Derivative)
    wx = w(1); wy = w(2); wz = w(3);
    Omega = [  0,  wz, -wy,  wx;
             -wz,   0,  wx,  wy;
              wy, -wx,   0,  wz;
             -wx, -wy, -wz,   0 ];
    dqdt = 0.5 * Omega * q;
    
    % 3. Dynamics (Euler's Equations)
    dwdt = J \ (tau - cross(w, J * w));
    
    % Pack state vector
    dxdt = [dqdt; dwdt];
end