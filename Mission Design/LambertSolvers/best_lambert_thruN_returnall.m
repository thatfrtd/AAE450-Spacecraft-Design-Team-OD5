function [vel1_reshaped, vel2_reshaped, v_filter] = best_lambert_thruN_returnall(x_1, x_2, ToF, N_max)
arguments
    x_1
    x_2
    ToF
    N_max
end
    % Solve Lambert up to a number of revolutions for positive and negative 
    % directions and evaluate all solutions with some arbitrary function

    % Solve Lambertus Maximus
    Q = size(x_1, 2);

    direction = [ones([Q, 1]); -ones([Q, 1])];

    % First find N (best number of revolutions)
    [v1vec, v2vec, ~, ~, ~] = ivLam_thruN_multipleInputDLL(2 * Q, repmat(x_1(1:3, :), 1, 2), repmat(x_2(1:3, :), 1, 2), repmat(ToF, 2, 1), direction, N_max);
    v1vec = cat(3, v1vec(:, 1 : size(v2vec, 2) / 2), v1vec(:, (size(v2vec, 2) / 2 + 1) : size(v2vec, 2)));
    v2vec = cat(3, v2vec(:, 1 : size(v2vec, 2) / 2), v2vec(:, (size(v2vec, 2) / 2 + 1) : size(v2vec, 2)));

    % Retrieve solutions
    [Ns, Qs] = meshgrid(0 : N_max, 1 : Q);
    jcolumn = Ni2col(Ns, Qs, N_max);

    % Filter out NaN and 0 solutions
    vel1_unfiltered = v1vec(1:3,jcolumn(:), :);
    vel2_unfiltered = v2vec(1:3,jcolumn(:), :);

    vel1_reshaped = reshape(vel1_unfiltered, 3, Q, N_max + 1, 2); % Double check this is correct
    vel2_reshaped = reshape(vel2_unfiltered, 3, Q, N_max + 1, 2);

    v_filter = all(~isnan(vel1_reshaped), 1) & any(vel1_reshaped ~= 0, 1) ...
               & all(~isnan(vel2_reshaped), 1) & any(vel2_reshaped ~= 0, 1);
end