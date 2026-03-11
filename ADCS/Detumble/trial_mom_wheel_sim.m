%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detumble Simulation V2
% Full reaction wheel internal coupling with saturation
% Momentum conserving model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ===================== Spacecraft Parameters ============================
m_r = 3000;      % kg (rocket)
m_s = 500;       % kg (satellite)
R   = 1.85;      % m
L   = 30;        % m
s   = 1.0;       % m

% Rocket inertia about own CM
J_r = diag([1/12*m_r*(3*R^2 + L^2), ...
            1/12*m_r*(3*R^2 + L^2), ...
            1/2*m_r*R^2]);

% Satellite inertia about own CM
J_s = diag([1/6*m_s*s^2, ...
            1/6*m_s*s^2, ...
            1/6*m_s*s^2]);

% Satellite position
z_s = -(L/2 + s/2);
z_cm = (m_r*0 + m_s*z_s) / (m_r + m_s);

d_r = abs(0 - z_cm);
d_s = abs(z_s - z_cm);

J_total = (J_r + m_r*diag([d_r^2 d_r^2 0])) + ...
          (J_s + m_s*diag([d_s^2 d_s^2 0]));

%% ===================== Reaction Wheel Parameters ========================
tau_max = 0.3;         % Nm (max torque per wheel)
h_max   = 1200;        % Nms (max stored momentum per wheel)
K       = 50;          % detumble gain

%% ===================== Initial Conditions ===============================
rpm = 2;
w0 = [0; 0; rpm * 2*pi/60];  % initial body spin

q0 = [0; 0; 0; 1];           % quaternion
h_rw0 = [0; 0; 0];           % initial wheel momentum

x0 = [q0; w0; h_rw0];

%% ===================== Simulation =======================================
tspan = [0 6000];

opts = odeset('RelTol',1e-9,'AbsTol',1e-9, ...
              'Events',@(t,x) stop_event(t,x));

[t,x] = ode45(@(t,x) plant_derivative(t,x,J_total,K,tau_max,h_max), ...
              tspan,x0,opts);

%% ===================== Extract States ===================================
q     = x(:,1:4);
omega = x(:,5:7);
h_rw  = x(:,8:10);

%% ===================== Plot Results =====================================
figure('Color','w')

subplot(3,1,1)
plot(t,omega,'LineWidth',1.5)
grid on
title('Angular Velocity')
ylabel('rad/s')
legend('\omega_x','\omega_y','\omega_z')

subplot(3,1,2)
plot(t,h_rw,'LineWidth',1.5)
grid on
title('Reaction Wheel Momentum')
ylabel('Nms')
legend('h_x','h_y','h_z')

subplot(3,1,3)
plot(t,q,'LineWidth',1.5)
grid on
title('Quaternion')
ylabel('Components')
xlabel('Time (s)')
legend('q_1','q_2','q_3','q_4')

%% Angular Momentum Calculations

H_body = zeros(length(t),3);
H_total = zeros(length(t),3);

for k = 1:length(t)
    omega_k = omega(k,:)';
    h_k     = h_rw(k,:)';

    H_body(k,:)  = (J_total * omega_k)';
    H_total(k,:) = (J_total * omega_k + h_k)';
end

%% Plot Angular Momentum

figure('Color','w')

subplot(2,1,1)
plot(t,H_body,'LineWidth',1.5)
grid on
title('Spacecraft Body Angular Momentum (J\omega)')
ylabel('N·m·s')
legend('H_x','H_y','H_z')

subplot(2,1,2)
plot(t,vecnorm(H_total,2,2),'LineWidth',1.5)
grid on
title('Total Angular Momentum Magnitude')
ylabel('|H_{total}|')
xlabel('Time (s)')

%% Angular Momentum Analysis

H_body  = zeros(length(t),3);
H_total = zeros(length(t),3);

for k = 1:length(t)
    omega_k = omega(k,:)';
    h_k     = h_rw(k,:)';

    H_body(k,:)  = (J_total * omega_k)';
    H_total(k,:) = (J_total * omega_k + h_k)';
end

figure('Color','w')

subplot(3,1,1)
plot(t,H_body,'LineWidth',1.5)
grid on
title('Body Angular Momentum (J\omega)')
ylabel('Nms')

subplot(3,1,2)
plot(t,h_rw,'LineWidth',1.5)
grid on
title('Reaction Wheel Momentum (h_{rw})')
ylabel('Nms')

subplot(3,1,3)
plot(t,vecnorm(H_total,2,2),'LineWidth',1.5)
grid on
title('Total Angular Momentum Magnitude (Should Be Constant)')
ylabel('|H|')
xlabel('Time (s)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdt = plant_derivative(~,x,J,K,tau_max,h_max)

    q     = x(1:4);
    omega = x(5:7);
    h_rw  = x(8:10);

    q = q / norm(q);

    %% ================= Controller (B-dot style detumble) ================
    tau_cmd = K * omega;

    % Torque saturation
    tau_rw = max(min(tau_cmd, tau_max), -tau_max);

    % Momentum saturation (per axis)
    for i = 1:3
        if abs(h_rw(i)) >= h_max && sign(tau_rw(i)) == sign(h_rw(i))
            tau_rw(i) = 0;
        end
    end

    %% ================= Dynamics =========================================
    % Spacecraft rigid body equation with internal wheel coupling:
    % J*w_dot = -tau_rw - w × (Jw + h_rw)

    dwdt = J \ (-tau_rw - cross(omega, J*omega + h_rw));

    % Wheel momentum dynamics
    dhdt = tau_rw;

    % Quaternion kinematics
    dqdt = quaternion_kde(q,omega);

    dxdt = [dqdt; dwdt; dhdt];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dqdt = quaternion_kde(q,omega)

    w1 = omega(1);
    w2 = omega(2);
    w3 = omega(3);

    Omega = [ 0   w3 -w2  w1;
             -w3  0   w1  w2;
              w2 -w1  0   w3;
             -w1 -w2 -w3  0];

    dqdt = 0.5 * Omega * q;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [value,isterminal,direction] = stop_event(~,x)

    omega = x(5:7);

    value = norm(omega) - 1e-4;  % stop when nearly zero
    isterminal = 1;
    direction = -1;
end