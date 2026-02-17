% SpaceX Falcon Heavy Analysis
% Assume full prop weight and full payload capacity
%

clear; clc; close all;

%% Vehicle Parameters
g = 9.81;                     % gravity (m/s^2)
L = 70;                       % rocket height (m)
m0 = 1.42e6;                 % initial mass (kg)
mp = 1.0e6;                  % propellant mass approx (kg)
burn_time = 150;             % first stage burn duration (s)
thrust = 22.8e6;             % thrust (N)

%% Structural properties (equivalent beam)
E = 70e9;                    % aluminum-lithium approx (Pa)
I = 18;                      % estimated second moment (m^4)
K = 3*E*I/L^3;               % effective bending stiffness

%% Damping ratio (typical rockets 1-5%)
zeta = 0.03;

%% Payload mass
mpayload = 63800;

%% Time Vector
dt = 0.01;
t = 0:dt:200;

%% Mass deplation
m = m0 - mp*(t/burn_time);
m(t>burn_time) = m0 - mp;

%% Axial Accleration
a = thrust./m - g;
a(a<0)=0;

%% MAx_q (aero load)
rho = 1.2*exp(-t/8);          % simple atmosphere model
V = cumtrapz(t,a);            % velocity
q = 0.5*rho.*V.^2;

Cd = 0.4;
A = pi*(3.7/2)^2;
Faero = q*Cd*A;

%% Bending mode model
% Single degree of freedom: Mx'' + Cx' + Kx = F(t)

M = m0/3;                    % modal mass approximation
C = 2*zeta*sqrt(K*M);

x = zeros(size(t));
xd = zeros(size(t));

for i=2:length(t)
    F = Faero(i);
    xdd = (F - C*xd(i-1) - K*x(i-1))/M;
    xd(i) = xd(i-1) + xdd*dt;
    x(i) = x(i-1) + xd(i)*dt;
end

%% Payload accleration
payload_accel = gradient(gradient(x,dt),dt);

%% PLots
figure;
plot(t,a/g,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Axial Acceleration (g)')
title('Rocket Axial g-Load During Launch')
grid on

figure;
plot(t,x,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Tip Deflection (m)')
title('Rocket Bending Deflection')
grid on

figure;
plot(t,payload_accel/g,'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Payload Acceleration (g)')
title('Payload Dynamic Load Environment')
grid on
