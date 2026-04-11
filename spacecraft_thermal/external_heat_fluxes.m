function Q = external_heat_fluxes(geom, surf, env, scenario, t_s, orbit)
% EXTERNAL_HEAT_FLUXES  Compute absorbed external heat at time t
%
%  DESCRIPTION:
%    For a given time t in the orbit, computes the total absorbed external
%    heat [W] on each spacecraft face from:
%      1. Direct solar radiation
%      2. Earth albedo (reflected solar)
%      3. Earth infrared (outgoing longwave radiation)
%
%    Also determines eclipse state (sun/shadow).
%
%  INPUTS:
%    geom     — geometry struct (from spacecraft_geometry.m)
%    surf     — surface optical properties struct
%    env      — environment struct (solar flux, albedo, Earth IR)
%    scenario — scenario struct (hot/cold case parameters)
%    t_s      — current time in orbit [s]
%    orbit    — orbit struct
%
%  OUTPUTS:
%    Q        — struct with fields:
%      .solar_W      — absorbed solar per face [W], vector (n_faces x 1)
%      .albedo_W     — absorbed albedo per face [W]
%      .earth_IR_W   — absorbed Earth IR per face [W]
%      .total_W      — total absorbed per face [W]
%      .total_sum_W  — scalar total [W]
%      .in_eclipse   — boolean
%      .sun_fraction — 0 (eclipse) to 1 (full sun), for partial cases
%
%  KEY EQUATIONS:
%    Q_solar  = alpha * S * A * cos(theta) * F_sun * (1 - eclipse)
%    Q_albedo = alpha * a * S * A * F_earth * f_albedo
%    Q_IR     = eps  * q_IR * A * F_earth
%
%    where F_sun/F_earth are view factors, a is Earth albedo fraction
%
%  TODO: Replace simplified view factor model with actual attitude-dependent
%        factors once ADCS pointing profile is known.
%  TODO: For high-fidelity model, integrate with SPICE or MATLAB Aerospace
%        Toolbox for true Sun/Earth geometry.
% =========================================================================

n_faces = geom.n_faces;

%% ---- Eclipse determination ---------------------------------------------
%  Simple model: fraction of orbit in shadow based on eclipse_frac
%  t=0 is start of sunlit period; eclipse occurs at end of orbit
%
%  TODO: For true geometry, use cylindrical shadow model:
%        eclipse when Earth-spacecraft-Sun angle places s/c in shadow cone

T_orbit = orbit.period;   % s
t_mod   = mod(t_s, T_orbit);
t_sun   = (1 - scenario.eclipse_frac) * T_orbit;  % sunlit duration

if t_mod < t_sun
    in_eclipse   = false;
    sun_fraction = 1.0;
else
    in_eclipse   = true;
    sun_fraction = 0.0;
end

Q.in_eclipse   = in_eclipse;
Q.sun_fraction = sun_fraction;

%% ---- Solar flux for this scenario ---------------------------------------
S = scenario.solar_flux;    % W/m^2

%% ---- Initialize output arrays ------------------------------------------
Q.solar_W    = zeros(n_faces, 1);
Q.albedo_W   = zeros(n_faces, 1);
Q.earth_IR_W = zeros(n_faces, 1);

%% ---- Albedo correction factor for altitude / Earth geometry ------------
%  Albedo peaks when face is near sub-solar point
%  Simple model: use mean albedo flux; scale by Earth view factor
%
%  More precise formula:
%    q_albedo = a * S * F_albedo(altitude, face_orientation)
%  TODO: Use ESARAD tables or Gauss-quadrature integration over Earth disk
%        for exact albedo distribution at given altitude and beta angle.

a = scenario.albedo;           % Earth albedo fraction (dimensionless)
q_earth_IR = env.q_earth_IR;   % W/m^2 — Earth outgoing LW

%% ---- Loop over faces ---------------------------------------------------
for i = 1:n_faces
    VFi   = geom.VF(i);
    face  = geom.faces(i);
    A     = face.area_m2;       % m^2

    % Get optical properties for this face
    switch face.coating
        case 'body'
            alpha  = surf.alpha_body;
            eps    = surf.eps_body;
        case 'panel'
            alpha  = surf.alpha_panel;
            eps    = surf.eps_panel;
        case 'panel_back'
            alpha  = 0.15;    % bare Al back  % TODO: confirm back coating
            eps    = 0.80;
        otherwise
            alpha  = surf.alpha_body;
            eps    = surf.eps_body;
    end

    %% Direct solar
    %  Q_solar = alpha * S * A * VF_sun * sun_fraction
    %  VF_sun accounts for cos(angle of incidence) already in view factor
    %  TODO: For variable attitude, VF_sun should be computed as
    %        cos(theta_sun) where theta_sun is angle between face normal
    %        and Sun vector — update when ADCS provides pointing schedule
    Q.solar_W(i) = alpha * S * A * VFi.sun * sun_fraction;

    %% Earth albedo
    %  Q_albedo = alpha * a * S * A * VF_earth * f_geometry
    %  f_geometry: albedo is highest for Earth-facing surfaces near sub-solar
    %  TODO: use proper albedo view factor from Eckart or Siegel & Howell
    Q.albedo_W(i) = alpha * a * S * A * VFi.earth * sun_fraction;

    %% Earth IR (present even in eclipse)
    %  Q_IR = eps * q_IR * A * VF_earth
    %  Earth IR does NOT depend on sun/eclipse for short-term analysis
    %  (Earth IR varies by ±10% season, ~±3% lat — use mean value)
    Q.earth_IR_W(i) = eps * q_earth_IR * A * VFi.earth;
end

%% ---- Totals ------------------------------------------------------------
Q.total_W    = Q.solar_W + Q.albedo_W + Q.earth_IR_W;
Q.total_sum_W = sum(Q.total_W);

% Breakdown for logging
Q.solar_total_W    = sum(Q.solar_W);
Q.albedo_total_W   = sum(Q.albedo_W);
Q.earth_IR_total_W = sum(Q.earth_IR_W);

end
