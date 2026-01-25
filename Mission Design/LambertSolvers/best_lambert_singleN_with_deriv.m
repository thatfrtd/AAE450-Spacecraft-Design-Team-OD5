function [vel1, vel2, dV, dv1v2dr1r2ToF, d2zdyT] = best_lambert_singleN_with_deriv(x_1, x_2, ToF, v1_assist, v2_assist, options)
arguments
    x_1 
    x_2 
    ToF 
    v1_assist 
    v2_assist 
    options.N = zeros([size(x_1, 2), 1])
    options.include_second_order = false 
end

    Q = size(x_1, 2);
    direction = ones([Q, 1]);
    
    % Solve Lambertus Maximus
    [vel1, vel2, infoReturnStatus, infoHalfRevStatus, dzdyT, d2zdyT] = ivLam_NtildeWithDerivs_multipleInputDLL(Q, x_1(1:3, :), x_2(1:3, :), ToF, direction, options.N, options.include_second_order);

    dv1v2dr1r2ToF = pagetranspose(dzdyT);

    % Calculate delta V
    dV = max(vecnorm(x_1(4: 6, :) - vel1) - v1_assist, 0) + max(vecnorm(x_2(4:6, :) - vel2) - v2_assist, 0);
end