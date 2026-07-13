function u_center = evaluate_DG_on_centers(coeff, num_cells, p_order)
%EVALUATE_DG_ON_CENTERS Recover DG solution at cell centers.
%   coeff    : global DG coefficient vector, length = num_cells*(p_order+1)
%   num_cells: number of cells
%   p_order  : polynomial degree
%
%   u_center : value on each cell center, size = [num_cells,1]

    nloc = p_order + 1;
    coeff_mat = reshape(coeff, nloc, num_cells).';

    u_center = zeros(num_cells,1);

    % 参考区间中心 xi = 0 处的基函数值
    phi0 = zeros(nloc,1);
    for k = 0:p_order
        phi0(k+1) = reference_basis(0, k);
    end

    for elem = 1:num_cells
        u_center(elem) = coeff_mat(elem,:) * phi0;
    end
end