function [passive_safety_constraint, min_safety, x_prop] = construct_passive_safety_constraint(x, x_0, f, A, safety_constraint, safety_constraint_linearized, T, N_safe, integration_tolerance)
%CONSTRUCT_PASSIVE_SAFETY_CONSTRAINT Constrains state at least safe time in safety time horizon to be safe
%   Uses state transition matrix to approximate mapping changes in initial
%   condition to state at worst safety
arguments
    x % State to evaluate (linearized) passive safety constraint 
    x_0 % Initial state (from reference trajectory)
    f % (t, x) Dynamics function
    A % (t, x) Jacobian of dynamics function
    safety_constraint % (t, x) <= 0 defines safe set
    safety_constraint_linearized % (t, x, x_ref_safe) 
    T % [s] Safety horizon - how long to propagate into the future to check safety
    N_safe % Number of evenly spaced time nodes to look at safety
    integration_tolerance = 1e-10
end

% Propagate for entire safety horizon
z_k = linspace(0, T, N_safe); % Time points to evaluate safety at
tolerances = odeset(RelTol=integration_tolerance, AbsTol=integration_tolerance);
[~, x_prop] = ode45(@(t_prop, x_prop) f(t_prop, x_prop), z_k, x_0, tolerances);
x_prop = x_prop';

% Find least safe time in horizon
safety = zeros([N_safe, 1]);
for k = 1 : N_safe  
    safety(k) = safety_constraint(z_k(k), x_prop(:, k));
end
[min_safety, min_safety_i] = max(safety);
z_min_safety = z_k(min_safety_i);
x_min_safety = x_prop(:, min_safety_i);
fprintf("min safety i: %g , min safety: %.3f\n", min_safety_i, min_safety)

% Calculate STM from initial to least safe time
A_min_safety = integrate_STM(x_0, f, A, [0, z_min_safety], tolerances);

% Define (approximate) safety constraint using STM to map perturbed initial
% to perturbed state at min safety
% Approximate because only doing the worst timestep
passive_safety_constraint = safety_constraint_linearized(z_min_safety, real(A_min_safety) * x, x_min_safety);

end

function [A_k] = integrate_STM(x0, f, A, tspan, tolerances)
    % Integrates STM and state
    %   Uses ODE45 to integrate the State Transition Matrix and the state using
    %   the given A matrix and dynamics f over the time period in tspan using
    %   the specified error tolerance.

    % Create initial condition
    nx = numel(x0);

    STM0 = eye(nx);
  
    y0 = [x0; STM0(:)];

    % Simulate    
    if tspan(1) ~= tspan(end)
        [~, y] = ode45(@(t, y) STM_diff_eq(t, y, A, f, nx), tspan, y0, tolerances);
    else
        y = y0';
    end

    y_f = y(end, :); 

    % Unpack solution
    A_k = reshape(y_f(:, (nx + 1) : (nx * (nx + 1))), nx, nx);
end

function [ydot] = STM_diff_eq(t, y, A, f, n)
    x = y(1:n);
    STM = reshape(y((n + 1) : (n * (n + 1))), n, n);

    xdot = f(t, x);
    A_kdot = A(t, x) * STM;

    ydot = [xdot; A_kdot(:)];
end