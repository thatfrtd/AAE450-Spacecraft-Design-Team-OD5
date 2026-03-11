%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Detumble Sim Analysis V1
%
% This is meant to simulate the spacecraft detumble process.
% This is V1. It's low fidelity with simple geometry and spin on one axis
% along with using only reaction wheels
% Last modified:
%   - 2/21
%   - Atharva
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init Variables
m_r = 3000; % kg
m_s = 1500; % kg
R = 1.85; % m
L = 30; % m
s = 1.0; % m

% Inertia of Rocket about its own CM
J_r = diag([1/12*m_r*(3*R^2 + L^2), 1/12*m_r*(3*R^2 + L^2), 1/2*m_r*R^2]);
% Inertia of Satellite about its own CM
J_s = diag([1/6*m_s*s^2, 1/6*m_s*s^2, 1/6*m_s*s^2]);

% Combined Center of Mass (relative to rocket center)
% Rocket center at 0, Sat center at -L/2 - s/2
% Case 1) Satellite at end of Rocket by nozzle:
z_s = -(L/2 + s/2);
z_cm = (m_r*0 + m_s*z_s) / (m_r + m_s);

% Parallel Axis Theorem to find J_total about combined CM
d_r = abs(0 - z_cm);
d_s = abs(z_s - z_cm);

J_total = (J_r + m_r*diag([d_r^2, d_r^2, 0])) + ...
          (J_s + m_s*diag([d_s^2, d_s^2, 0]));

% Initial Conditions
rpm = 2;
w0 = [0; 0; rpm * 2 * pi / 60]; % Spin along long axis (z)
q0 = [0; 0; 0; 1];             % Identity quaternion
x0 = [q0; w0];

%% Control
% Typical Reaction Wheel torque is ~0.1 to 0.4 Nm
u_const = [0; 0; -0.3]; % const torque on z axis

%% ODE45
tspan = [0 5000]; % seconds
opts = odeset('Events', @stop_event, 'RelTol', 1e-8, 'AbsTol', 1e-8);
[t, x] = ode45(@(t, x) plant_derivative(t, x, J_total, u_const), tspan, x0, opts);

%% Plotting
% Extract states
q = x(:, 1:4);
w = x(:, 5:7);

figure('Color', 'w');
subplot(2,1,1);
plot(t, w, 'LineWidth', 1.5);
grid on; title('Angular Rates (\omega)');
ylabel('rad/s'); legend('\omega_x', '\omega_y', '\omega_z');

subplot(2,1,2);
plot(t, q, 'LineWidth', 1.5);
grid on; title('Quaternion Orientation (\epsilon)');
ylabel('Components'); legend('\epsilon_1', '\epsilon_2', '\epsilon_3', '\eta');
xlabel('Time (s)');


%% Helper Functions
function dxdt = plant_derivative(t, x, J, u_const)
    epsilon = x(1:4);
    omega   = x(5:7);

    % Re-normalize to prevent numerical drift
    epsilon = epsilon / norm(epsilon);

    % Calculate Kinematics 
    deps = quatern  ion_kde(t, epsilon, omega);

    % Calculate Dynamics 
    dwdt = rigid_body_dynamics(omega, J, u_const);

    % Pack back into state vector
    dxdt = [deps; dwdt];
end

function [value, isterminal, direction] = stop_event(~, x)
    % value: the quantity we want to reach zero (omega_z)
    % isterminal: 1 = stop the integration
    % direction: -1 = only trigger if approaching from positive side
    value = x(7); 
    isterminal = 1; 
    direction = -1; 
end