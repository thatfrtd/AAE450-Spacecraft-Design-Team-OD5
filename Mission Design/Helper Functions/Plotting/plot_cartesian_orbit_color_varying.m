function [] = plot_cartesian_orbit_color_varying(x_cartesian, color, n_colors)
%PLOT_CARTESIAN_ORBIT_COLOR_VARYING Plot 3D orbit with color varying along it
%   Reference: https://www.mathworks.com/matlabcentral/answers/254696-how-to-assign-gradual-color-to-a-3d-line-based-on-values-of-another-vector
% Used for plotting orbits with thrusting arcs, coasting arcs, eclipsing, 
% etc
arguments
    x_cartesian
    color
    n_colors = []
end

if isempty(n_colors)
    C = cool();
else
    C = cool(n_colors);
end
colormap(C);
patch([x_cartesian(1, :) nan],[x_cartesian(2, :) nan],[x_cartesian(3, :) nan],[color nan],'FaceColor','none','EdgeColor','flat')
% How to label???????????? Make custom color map??
view(3)

end