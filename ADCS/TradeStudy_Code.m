%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Trade Study Code
%
% This is meant to simulate ACS System Capabilities.
% This is V1. It's low fidelity with simple geometry and spin on one axis
% along with using only reaction wheels and CMGs
% Last modified:
%   - 4/9
%   - Atharva
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_s = 2000; % kg
r = 2.0; % m radius of the cylinder
h = 10.0; % m height/length of the cylinder

I_xx = (1/12) * m_s * (3*r^2 + h^2);
I_yy = I_xx; 
I_zz = (1/2) * m_s * r^2;

I_sc = diag([I_xx, I_yy, I_zz]);
I_sc = I_xx;
% Mission Parameters
slew_angle_deg = 90; % Maneuver angle (degrees)
theta = deg2rad(slew_angle_deg);

% Target Max Slew Rate for Part 3
desired_omega_deg = 2.0; % Target slew rate (deg/s)
desired_omega = deg2rad(desired_omega_deg);

% Typical Hardware Specs (Torque in Nm, Momentum in Nms)
RW_tau = 0.2; % Reaction Wheel 0.1 to 1.0 Nm Honeywell Aerospace "Space Momentum Products" Catalog
RW_H_max = 12; % 10 to 150 Nms     

CMG_tau = 15; % CMG Airbus Defence and Space Agile Actuators Datasheets (specifically the CMG 15-45S and CMG 75-75S hardware specifications).   
CMG_H_max = 75; % 15 to 75 Nms

RCS_tau = 44; % Our Current Thruster Abilities
RCS_H_max = 10000; % Max Momentum for RCS Thrusters (Limited not by value, but by total fuel weight. this increases based off how much we want to use)

fprintf('--- 1. Max Slew Rate Possible Given Torque & MOI ---\n');
[rw_max_w, rw_lim] = calc_max_slew(I_sc, RW_tau, theta, RW_H_max);
fprintf('Reaction Wheel: %.4f deg/s (Constrained by %s)\n', rad2deg(rw_max_w), rw_lim);

[cmg_max_w, cmg_lim] = calc_max_slew(I_sc, CMG_tau, theta, CMG_H_max);
fprintf('Control Moment Gyro:  %.4f deg/s (Constrained by %s)\n', rad2deg(cmg_max_w), cmg_lim);

[rcs_max_w, rcs_lim] = calc_max_slew(I_sc, RCS_tau, theta, RCS_H_max);
fprintf('Reaction Control System:  %.4f deg/s (Constrained by %s)\n\n', rad2deg(rcs_max_w), rcs_lim);

fprintf('--- 2. Torque Needed for Desired Max Slew Rate (%.1f deg/s) ---\n', desired_omega_deg);
rw_req_tau = calc_req_torque(I_sc, desired_omega, theta, RW_H_max);
fprintf('Reaction Wheel: %s\n', rw_req_tau);

cmg_req_tau = calc_req_torque(I_sc, desired_omega, theta, CMG_H_max);
fprintf('Control Moment Gyro:  %s\n', cmg_req_tau);

rcs_req_tau = calc_req_torque(I_sc, desired_omega, theta, RCS_H_max);
fprintf('Reaction Control System:  %s\n\n', rcs_req_tau);






function [actual_max_w, limiting_factor] = calc_max_slew(I, tau, theta, H_max)
    % Peak angular velocity based purely on kinematics (Bang-Bang)
    w_kinematic = sqrt((tau * theta) / I);
    
    % Maximum angular velocity based on momentum capacity
    w_momentum = H_max / I;
    
    if w_kinematic < w_momentum
        actual_max_w = w_kinematic;
        limiting_factor = 'Kinematics/Torque';
    else
        actual_max_w = w_momentum;
        limiting_factor = 'Momentum Saturation (H_max)';
    end
end

function result_str = calc_req_torque(I, w_max, theta, H_max)
    % First, check if the system can even hold enough momentum to reach w_max
    if (H_max / I) < w_max
        result_str = sprintf('FAILED - Cannot reach rate. Exceeds momentum capacity limit (Max capable: %.2f deg/s)', rad2deg(H_max/I));
        return;
    end
    
    % If momentum is sufficient, calculate required torque for a bang-bang maneuver
    tau_req = (I * w_max^2) / theta;
    result_str = sprintf('Requires %.2f Nm of torque.', tau_req);
end