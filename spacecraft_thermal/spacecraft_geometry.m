function geom = spacecraft_geometry()
% SPACECRAFT_GEOMETRY  Define spacecraft shape, surface areas, and view factors
%
%  DESCRIPTION:
%    Models the spacecraft as a cylinder (main bus) + two rectangular
%    solar panels + a bottom thruster plate.
%    All areas are computed analytically. View factors to Earth, Sun,
%    and deep space are estimated per face based on attitude assumptions.
%
%  GEOMETRY (as given):
%    Cylinder: D = 3.0 m, H = 3.2 m
%    Solar panels (x2): 3.2 m x 9.0 m each
%    Top face: arm/docking — not available for radiators % TODO: confirm arm geometry
%    Bottom face: ion thruster — not available for radiators
%
%  OUTPUT:
%    geom — struct with fields:
%      .D, .H, .R           — cylinder dimensions (m)
%      .A_side              — lateral cylindrical area (m^2)
%      .A_top, .A_bottom    — end cap areas (m^2)
%      .A_panel_one_side    — area of one solar panel face (m^2)
%      .A_panels_total      — both panels, both sides total (m^2)
%      .faces               — struct array of named surface faces
%      .VF                  — view factor estimates per face
%
%  NOTE ON ATTITUDE:
%    The spacecraft attitude is variable. For worst-case analysis:
%      HOT  — face with largest area pointed at Sun
%      COLD — Sun-facing area minimized; Earth-facing maximized
%    TODO: Once ADCS team confirms pointing modes, update VF_sun per face.
% =========================================================================

%% ---- Main bus (cylinder) ------------------------------------------------
geom.D = 3.0;           % m — diameter
geom.H = 3.2;           % m — height (excl. panels)
geom.R = geom.D / 2;    % m — radius

geom.A_top    = pi * geom.R^2;               % m^2
geom.A_bottom = pi * geom.R^2;               % m^2
geom.A_side   = pi * geom.D * geom.H;        % m^2 — full lateral surface

%% ---- Solar panels -------------------------------------------------------
panel_width  = 3.2;   % m (same as body height, attached along length)
panel_length = 9.0;   % m
geom.A_panel_one_face  = panel_width * panel_length;   % m^2, one side, one panel
geom.A_panels_total    = 4 * geom.A_panel_one_face;    % both panels, both sides

%% ---- Available radiator surfaces ----------------------------------------
%  TOP    — arm present, NOT available for radiators  % TODO: confirm with mechanical
%  BOTTOM — ion thruster, NOT available for radiators % TODO: confirm plume heating
%  SIDES  — partially blocked by RCS thrusters
%    Assume 70% of side area is available for passive radiation
%    TODO: Get actual RCS blockage fraction from propulsion team

geom.radiator_side_fraction = 0.70;  % fraction of side area usable  % TODO: update
geom.A_radiator_available   = geom.A_side * geom.radiator_side_fraction;

fprintf('  [Geometry] Cylinder side area:       %.2f m^2\n', geom.A_side);
fprintf('  [Geometry] Available radiator area:  %.2f m^2\n', geom.A_radiator_available);
fprintf('  [Geometry] Each solar panel face:    %.2f m^2\n', geom.A_panel_one_face);

%% ---- Face definitions (for view factor assignment) ----------------------
% Each face has: name, area, normal vector (body frame), coating type
% Body frame: +Z = top (arm), -Z = bottom (thruster), +X/-X = panel side

faces(1).name     = 'side_cylindrical';
faces(1).area_m2  = geom.A_side;
faces(1).normal   = [0; 1; 0];   % average — cylinder integrates over phi
faces(1).coating  = 'body';
faces(1).available_radiator = true;

faces(2).name     = 'top_cap';
faces(2).area_m2  = geom.A_top;
faces(2).normal   = [0; 0; 1];   % +Z
faces(2).coating  = 'body';
faces(2).available_radiator = false;   % arm present

faces(3).name     = 'bottom_cap';
faces(3).area_m2  = geom.A_bottom;
faces(3).normal   = [0; 0; -1];  % -Z
faces(3).coating  = 'body';
faces(3).available_radiator = false;   % thruster

faces(4).name     = 'solar_panel_front';
faces(4).area_m2  = 2 * geom.A_panel_one_face;  % both panels, sun-facing side
faces(4).normal   = [1; 0; 0];   % panels track sun — normal toward Sun
faces(4).coating  = 'panel';
faces(4).available_radiator = false;   % solar cells on front

faces(5).name     = 'solar_panel_back';
faces(5).area_m2  = 2 * geom.A_panel_one_face;  % back side
faces(5).normal   = [-1; 0; 0];
faces(5).coating  = 'panel_back';     % TODO: confirm back coating (bare Al? OSR?)
faces(5).available_radiator = true;   % backs could serve as passive radiator

geom.faces = faces;
geom.n_faces = length(faces);

%% ---- View Factors to Sun / Earth / Space --------------------------------
%  VF_sun(i)   = fraction of face i's hemisphere subtended by Sun direction
%  VF_earth(i) = fraction subtended by Earth disk
%  VF_space(i) = remainder (deep space)
%
%  Approximate Earth view factor from LEO using:
%    F_earth = sin^2(rho)  where rho = asin(R_earth / R_orbit)
%
R_earth  = 6371e3;  % m
R_orbit  = R_earth + 1100e3;  % m — using 1100 km nominal
rho_rad  = asin(R_earth / R_orbit);
VF_earth_hemisphere = sin(rho_rad)^2;   % ~0.34 for 1100 km

fprintf('  [Geometry] Earth view factor (hemisphere): %.3f\n', VF_earth_hemisphere);

%  TODO: For a proper model, replace these with ESARAD/Thermal Desktop output
%        or a Monte Carlo ray-tracing result for your actual pointing modes.

%  For now: use simple bounding estimates per face
%  Attitude assumption: nadir-pointing (bottom toward Earth, top toward space)

VF(1).face = 'side_cylindrical';
VF(1).sun   = 1/pi;    % cylinder sees ~1/pi of solar hemisphere on average
VF(1).earth = VF_earth_hemisphere * 0.5;
VF(1).space = 1 - VF(1).sun - VF(1).earth;

VF(2).face = 'top_cap';
VF(2).sun   = 0.0;     % top faces zenith — back to Sun in nadir-point
VF(2).earth = 0.0;     % top faces away from Earth
VF(2).space = 1.0;

VF(3).face = 'bottom_cap';
VF(3).sun   = 0.0;
VF(3).earth = VF_earth_hemisphere;
VF(3).space = 1 - VF(3).earth;

VF(4).face = 'solar_panel_front';
VF(4).sun   = 1.0;     % panels track Sun — front always faces Sun
VF(4).earth = 0.0;
VF(4).space = 0.0;

VF(5).face = 'solar_panel_back';
VF(5).sun   = 0.0;
VF(5).earth = 0.0;
VF(5).space = 1.0;     % backs face deep space in sun-tracking config

geom.VF = VF;
fprintf('  [Geometry] View factors assigned.\n');

%% ---- Thermal Mass (rough estimates for transient) -----------------------
%  TODO: Replace with actual mass budget from structures team
%  Material: Aluminum 6061-T6: rho=2700 kg/m^3, Cp=896 J/kg/K, k=167 W/m/K
geom.material.name  = 'Al_6061';
geom.material.rho   = 2700;    % kg/m^3    % TODO: confirm alloy
geom.material.Cp    = 896;     % J/(kg·K)  % TODO: confirm
geom.material.k     = 167;     % W/(m·K)   % TODO: confirm
geom.material.wall_thickness = 0.003;  % m — 3 mm wall  % TODO: from structures

% Shell mass (cylinder side + caps)
A_shell   = geom.A_side + geom.A_top + geom.A_bottom;
geom.m_shell_kg = A_shell * geom.material.wall_thickness * geom.material.rho;
fprintf('  [Geometry] Estimated shell mass: %.1f kg\n', geom.m_shell_kg);

end
