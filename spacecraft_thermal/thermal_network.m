function [nodes, C_vec, G_matrix] = thermal_network(geom, components, surf)
% THERMAL_NETWORK  Build the lumped-parameter thermal network
%
%  DESCRIPTION:
%    Constructs the thermal network for transient analysis.
%    Each spacecraft component + the shell is a thermal "node."
%    Nodes are coupled by:
%      - Conduction (solid connections, thermal straps/heat pipes)
%      - Radiation (internal cavity radiation between nodes)
%
%    The governing equation for each node i is:
%
%      C_i * dT_i/dt = Q_ext_i + Q_int_i
%                    + sum_j [ G_cond_ij * (T_j - T_i) ]
%                    + sum_j [ G_rad_ij  * (T_j^4 - T_i^4) ]
%                    - Q_rad_to_space_i
%
%    where:
%      C_i       = thermal capacitance [J/K]
%      G_cond_ij = conductive conductance [W/K]
%      G_rad_ij  = radiative coupling [W/K^4] (= eps_eff * sigma * A * F_ij)
%      Q_rad_to_space = radiation from node to deep space [W]
%
%  NODE MAP:
%    Nodes 1..N_comp = components (from component_heat_map.m)
%    Node  N_comp+1  = spacecraft structure/shell
%
%  OUTPUTS:
%    nodes    — struct array, one per node (name, T_init_K, C_J_K, etc.)
%    C_vec    — thermal capacitance vector [J/K], length N_nodes
%    G_matrix — conductive conductance matrix [W/K], N_nodes x N_nodes
%                (symmetric; diagonal = sum of off-diagonal for energy balance)
%
%  TODO: Replace contact conductance values with measured interface
%        conductances once mechanical drawings are finalized.
%  TODO: Add heat pipe / thermal strap conductances when routing is done.
%  TODO: Internal radiative couplings are set to 0 for now (conservative)
%        — add when interior surface map is available.
% =========================================================================

sigma = 5.6704e-8;  % Stefan-Boltzmann constant [W/m^2/K^4]
N_comp = length(components);
N_nodes = N_comp + 1;   % +1 for shell node
i_shell = N_nodes;      % index of shell node

%% ---- Build node array --------------------------------------------------
for i = 1:N_comp
    nodes(i).name        = components(i).name;
    nodes(i).location    = components(i).location;
    nodes(i).T_init_K    = 293.15;   % 20°C — initial guess; converges in transient
    nodes(i).C_J_K       = components(i).C_thermal_J_K;
    nodes(i).T_min_op_K  = components(i).T_op_min_C  + 273.15;
    nodes(i).T_max_op_K  = components(i).T_op_max_C  + 273.15;
    nodes(i).T_min_surv_K = components(i).T_surv_min_C + 273.15;
    nodes(i).T_max_surv_K = components(i).T_surv_max_C + 273.15;
    nodes(i).is_exterior = strcmp(components(i).location(end-7:end), 'exterior') | ...
                           contains(components(i).location, 'exterior');
end

% Shell node
nodes(i_shell).name       = 'Structure_Shell';
nodes(i_shell).location   = 'exterior';
nodes(i_shell).T_init_K   = 293.15;
nodes(i_shell).C_J_K      = geom.m_shell_kg * geom.material.Cp;
nodes(i_shell).T_min_op_K = (-60 + 273.15);   % °C
nodes(i_shell).T_max_op_K = (120 + 273.15);
nodes(i_shell).T_min_surv_K = (-80 + 273.15);
nodes(i_shell).T_max_surv_K = (150 + 273.15);
nodes(i_shell).is_exterior = true;

%% ---- Capacitance vector ------------------------------------------------
C_vec = zeros(N_nodes, 1);
for i = 1:N_nodes
    C_vec(i) = nodes(i).C_J_K;
end

%% ---- Conductance matrix ------------------------------------------------
%  G_matrix(i,j) = conductive coupling between nodes i and j [W/K]
%  G_matrix is symmetric; diagonal elements not used (computed in ODE solver)
G_matrix = zeros(N_nodes, N_nodes);

%% ---- Interface conductances --------------------------------------------
%  Contact conductance formula: G = k * A / L  (parallel paths sum)
%  OR:  G = h_contact * A  where h_contact is interface conductance [W/m^2/K]
%
%  Typical values for bolted aluminum interfaces:
%    h_contact_bolted = 1000–5000 W/m^2/K (with thermal interface material)
%    h_contact_bolted =  500–1000 W/m^2/K (bare metal, low contact pressure)
%  TODO: Replace with actual fastener patterns and TIM specs

h_bolted_TIM = 2000;   % W/m^2/K — with Bergquist GP3000 or similar  % TODO: update
h_bolted_dry = 500;    % W/m^2/K — dry bolted  % TODO: update

% Contact area estimates per component mount
% TODO: Get actual mounting footprints from mechanical team

% Component index shortcuts for readability
idx_thruster  = 1;
idx_panels    = 2;
idx_battery   = 3;
idx_camera    = 4;
idx_st        = 5;   % star tracker
idx_ant       = 6;   % antennas
idx_cmg       = 7;
idx_fcomp     = 8;   % flight computer
idx_imu       = 9;
idx_rcs       = 10;

% Ion thruster — bolted to bottom plate of shell
A_contact_thruster = 0.05;   % m^2  % TODO: from thruster mounting flange
G_matrix(idx_thruster, i_shell) = h_bolted_TIM * A_contact_thruster;
G_matrix(i_shell, idx_thruster) = G_matrix(idx_thruster, i_shell);

% Battery — mounted on internal shelf (good thermal contact desired)
A_contact_battery = 0.015;   % m^2  % TODO: battery mounting plate area
G_matrix(idx_battery, i_shell) = h_bolted_TIM * A_contact_battery;
G_matrix(i_shell, idx_battery) = G_matrix(idx_battery, i_shell);

% Camera — externally mounted
A_contact_camera = 0.003;    % m^2  % TODO
G_matrix(idx_camera, i_shell) = h_bolted_dry * A_contact_camera;
G_matrix(i_shell, idx_camera) = G_matrix(idx_camera, i_shell);

% Star tracker — externally mounted, low-conductance (for stability)
A_contact_st = 0.002;   % m^2  % TODO
G_matrix(idx_st, i_shell) = h_bolted_dry * A_contact_st;
G_matrix(i_shell, idx_st) = G_matrix(idx_st, i_shell);

% Antennas
A_contact_ant = 0.003;  % m^2  % TODO
G_matrix(idx_ant, i_shell) = h_bolted_dry * A_contact_ant;
G_matrix(i_shell, idx_ant) = G_matrix(idx_ant, i_shell);

% CMGs — internally mounted, good contact
A_contact_cmg = 0.02;   % m^2  % TODO: CMG flange area × 4 units
G_matrix(idx_cmg, i_shell) = h_bolted_TIM * A_contact_cmg;
G_matrix(i_shell, idx_cmg) = G_matrix(idx_cmg, i_shell);

% Flight computer
A_contact_fcomp = 0.005;  % m^2
G_matrix(idx_fcomp, i_shell) = h_bolted_TIM * A_contact_fcomp;
G_matrix(i_shell, idx_fcomp) = G_matrix(idx_fcomp, i_shell);

% IMU
A_contact_imu = 0.003;  % m^2
G_matrix(idx_imu, i_shell) = h_bolted_dry * A_contact_imu;
G_matrix(i_shell, idx_imu) = G_matrix(idx_imu, i_shell);

% RCS — attached to exterior brackets
A_contact_rcs = 0.004 * 16;  % m^2 per thruster × 16  % TODO
G_matrix(idx_rcs, i_shell) = h_bolted_dry * A_contact_rcs;
G_matrix(i_shell, idx_rcs) = G_matrix(idx_rcs, i_shell);

% Solar panels — connected through hinged boom (low conductance)
% TODO: Add actual panel-to-bus hinge/boom thermal conductance
G_panel_boom = 0.5;   % W/K — very rough estimate for hinged deployment  % TODO
G_matrix(idx_panels, i_shell) = G_panel_boom;
G_matrix(i_shell, idx_panels) = G_panel_boom;

%% ---- Thermal Straps / Heat Pipes (placeholder) -------------------------
%  TODO: Once routing decisions are made, add explicit strap conductances
%  Typical Al thermal strap: G ~ 2–10 W/K
%  Typical heat pipe (1/4" Al): G ~ 50–500 W/K (effective)
%
%  For now: no straps. Conduction through shell is the only path.
%  If analysis shows components overheating, add straps here.

%  Example (commented out — uncomment when strap routing is defined):
%  G_strap_battery_to_shell = 5.0;  % W/K  % TODO: design and specify
%  G_matrix(idx_battery, i_shell) += G_strap_battery_to_shell;
%  G_matrix(i_shell, idx_battery) += G_strap_battery_to_shell;

fprintf('  [Network] Thermal network built: %d nodes, %d couplings\n', ...
    N_nodes, nnz(triu(G_matrix,1)));
fprintf('  [Network] NOTE: Thermal straps/heat pipes NOT included yet.\n');
fprintf('  [Network] TODO: Add straps/pipes once routing is decided.\n');

end
