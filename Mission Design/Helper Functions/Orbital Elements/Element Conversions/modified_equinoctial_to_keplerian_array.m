function [x_keplerian, nu] = modified_equinoctial_to_keplerian_array(x_modified_equinoctial_array)
    x_keplerian = zeros(6, size(x_modified_equinoctial_array, 2));
    nu = zeros(size(x_modified_equinoctial_array, 2), 1);
    
    for ind = 1:size(x_modified_equinoctial_array, 2)
        x_modified_equinoctial = x_modified_equinoctial_array(:, ind);

        [x_keplerian(:, ind), nu(ind)] = modified_equinoctial_to_keplerian(x_modified_equinoctial);
    end
end