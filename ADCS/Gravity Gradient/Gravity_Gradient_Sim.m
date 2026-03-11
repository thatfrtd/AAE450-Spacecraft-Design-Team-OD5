%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Detumble Sim Analysis V3 (Reaction Wheel Saturation)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Init Variables
m_r = 3000; % kg
m_s = 500;  % kg
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

%% Orbital Parameters (Assuming Circular LEO)
mu = 3.986e14;       % Earth's gravitational parameter (m^3/s^2)
R_earth = 6378e3;    % Earth radius (m)
altitude = 200e3;    % 200 km orbit
R_c = R_earth + altitude; 
n = sqrt(mu / R_c^3); % Mean motion (rad/s)

%% Reaction Wheel Parameters
H_max = 20;         % Max momentum storage capacity (Nms) - adjust as needed
tau_max = 0.4;      % Max torque output per axis (Nm)

%% Initial Conditions
w0 = [0; 0; 0]; 
q0 = [1; 1; 1; 0];  
q0 = q0 / norm(q0); % Ensure initial quaternion is normalized
h_rw0 = [0; 0; 0];  % Initial momentum stored in the reaction wheels (Nms)

x0 = [q0; w0; h_rw0]; % Expanded state vector (10 elements)

%% ODE45 Propagation
tspan = [0 50000]; % Extend time to ensure we hit saturation
opts = odeset('Events', @(t, x) stop_event(t, x, H_max), 'RelTol', 1e-8, 'AbsTol', 1e-8);
[t, x] = ode45(@(t, x) plant_derivative(t, x, J_total, n, tau_max), tspan, x0, opts);

%% Plotting
q = x(:, 1:4);
w = x(:, 5:7);
h_rw = x(:, 8:10);

% Calculate the magnitude of the momentum vector over time
h_mag = sqrt(h_rw(:,1).^2 + h_rw(:,2).^2 + h_rw(:,3).^2);

figure('Color', 'w', 'Position', [100, 100, 800, 800]);

subplot(3,1,1);
plot(t, w, 'LineWidth', 1.5);
grid on; title('Spacecraft Angular Rates (\omega)');
ylabel('rad/s'); legend('\omega_x', '\omega_y', '\omega_z');

subplot(3,1,2);
plot(t, q, 'LineWidth', 1.5);
grid on; title('Quaternion Orientation (\epsilon, \eta)');
ylabel('Components'); legend('\epsilon_1', '\epsilon_2', '\epsilon_3', '\eta');

subplot(3,1,3);
plot(t, h_rw, 'LineWidth', 1.5); hold on;
plot(t, h_mag, 'k--', 'LineWidth', 2);
yline(H_max, 'r:', 'LineWidth', 2); yline(-H_max, 'r:', 'LineWidth', 2);
grid on; title('Reaction Wheel Momentum Build-Up');
ylabel('Momentum (Nms)'); legend('h_x', 'h_y', 'h_z', '|h_{total}|', 'Saturation Limit');
xlabel('Time (s)');

figure
plot(t, h_rw, 'LineWidth', 1.5); hold on;
plot(t, h_mag, 'k--', 'LineWidth', 2);
yline(H_max, 'r:', 'LineWidth', 2); yline(-H_max, 'r:', 'LineWidth', 2);
grid on; title('Reaction Wheel Momentum Build-Up');
ylabel('Momentum (Nms)'); legend('h_x', 'h_y', 'h_z', '|h_{total}|', 'Saturation Limit');
xlabel('Time (s)');
title("Reaction Wheel Saturation Time")

fprintf('Simulation stopped at t = %.1f seconds.\n', t(end));
fprintf('Reason: Reaction wheels saturated (Total momentum hit %d Nms).\n', H_max);

%% Helper Functions
function dxdt = plant_derivative(t, x, J, n, tau_max)
    % Extract states
    q = x(1:4);
    w = x(5:7);
    h_rw = x(8:10);
    
    % Re-normalize to prevent numerical drift
    q = q / norm(q); 
    
    % 1. Kinematics
    wx = w(1); wy = w(2); wz = w(3);
    Omega = [  0,  wz, -wy,  wx;
             -wz,   0,  wx,  wy;
              wy, -wx,   0,  wz;
             -wx, -wy, -wz,   0 ];
    dqdt = 0.5 * Omega * q;
    
    % 2. Gravity Gradient Torque
    u_ECI = [-cos(n*t); -sin(n*t); 0];
    
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    R_ECI2Body = [1-2*(q2^2+q3^2),  2*(q1*q2+q3*q4),  2*(q1*q3-q2*q4);
                  2*(q1*q2-q3*q4),  1-2*(q1^2+q3^2),  2*(q2*q3+q1*q4);
                  2*(q1*q3+q2*q4),  2*(q2*q3-q1*q4),  1-2*(q1^2+q2^2)];
              
    u_Body = R_ECI2Body * u_ECI;
    tau_gg = 3 * n^2 * cross(u_Body, J * u_Body);
    
    % 3. Control Law (Active gravity gradient cancellation + rate damping)
    % The wheels try to perfectly oppose the GG torque and stop any spin
    tau_cmd = -tau_gg - 50 * w; 
    
    % Saturate the commanded torque to physical wheel limits
    tau_rw = max(min(tau_cmd, tau_max), -tau_max);
    
    % 4. Wheel Dynamics
    % Momentum change of the wheel is opposite to the torque it applies to the spacecraft
    dh_rw = -tau_rw; 
    
    % 5. Spacecraft Dynamics (Euler's Eq modified for internal momentum)
    % J*w_dot + w x (J*w + h_rw) = tau_ext + tau_rw
    dwdt = J \ (tau_gg + tau_rw - cross(w, J * w + h_rw));
    
    % Pack state vector
    dxdt = [dqdt; dwdt; dh_rw];
end

function [value, isterminal, direction] = stop_event(~, x, H_max)
    % Extract current RW momentum
    h_rw = x(8:10);
    current_momentum = norm(h_rw);
    
    % Trigger event when current_momentum reaches H_max
    value = current_momentum - H_max; 
    isterminal = 1; % Stop the integration
    direction = 1;  % Only trigger when crossing from negative to positive
end