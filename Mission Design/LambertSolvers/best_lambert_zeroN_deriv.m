function [vel1, vel2, dV, dpdz, d2pdz] = best_lambert_zeroN_deriv(x_1, x_2, ToF, v1_assist, v2_assist, options)
    arguments
        x_1
        x_2
        ToF
        v1_assist
        v2_assist
        options.include_second_order = false
        options.direction = 1
    end
    % z is the lambert inputs: [r1; r2; ToF]
    % p is the lambert outputs: [v1; v2]

    Q = size(x_1, 2);
    direction = options.direction * ones([Q, 1]);
    N_rev = zeros([Q, 1]);

    % Solve Lambertus Maximus
    [vel1, vel2, ~, ~, dpdz_transposed, d2pdz] = ivLam_NtildeWithDerivs_multipleInputDLL(Q, x_1(1:3, :), x_2(1:3, :), ToF, direction, N_rev, options.include_second_order);

    dpdz = pagetranspose(dpdz_transposed);

    % Calculate delta V
    dV = max(vecnorm(x_1(4: 6, :) - vel1) - v1_assist, 0) + max(vecnorm(x_2(4:6, :) - vel2) - v2_assist, 0);
end