function [phi_coeff, E_coeff] = solve_1D_LDG(mesh, basis_data, carrier_coeff, time_now, params, poisson_data)
%SOLVE_1D_LDG Solve the Poisson subproblem with cached matrices.

    cell_sizes = mesh.cell_sizes;
    cell_centers = mesh.cell_centers;
    num_cells = mesh.num_cells;
    problem = params.problem;

    p_order = basis_data.p_order;
    inv_mass = basis_data.inv_mass(:);
    Gp = basis_data.gauss_points;
    Gw = basis_data.gauss_weights;
    basis_values = basis_data.basis_values;
    face_values = basis_data.face_values;

    M1 = poisson_data.main;
    M2 = poisson_data.aux;
    C_p = p_order / max(cell_sizes);

    carrier_coeff = reshape(carrier_coeff, p_order+1, num_cells);

    poisson_scale = problem.physics.poisson_scale;
    left_boundary_value = problem.boundary.phi_left(time_now);
    right_boundary_value = problem.boundary.phi_right(time_now);

    b1 = zeros(num_cells * (p_order+1), 1);
    left_block = -inv_mass .* face_values(1,:).' / cell_sizes(1) * left_boundary_value;
    right_block =  inv_mass .* face_values(2,:).' / cell_sizes(num_cells) * right_boundary_value;
    b1(1:p_order+1) = left_block;
    b1(end-p_order:end) = right_block;

    b2 = zeros(num_cells * (p_order+1), 1);
    b2(end-p_order:end) = C_p * right_boundary_value * face_values(2,:).';

    x = Gp * cell_sizes + cell_centers;
    carrier_values = basis_values * carrier_coeff;
    rhs_values = carrier_values - problem.source.poisson(x, time_now);

    f_terms = basis_values' * bsxfun(@times, Gw, rhs_values);
    f_terms = bsxfun(@times, f_terms, cell_sizes);
    f_terms = poisson_scale * f_terms(:) - M2 * b1 - b2;

    phi_coeff_vec = solve_with_linear_solver(poisson_data.solver, f_terms);
    E_coeff_vec = M1 * phi_coeff_vec + b1;

    phi_coeff = reshape(phi_coeff_vec, p_order+1, num_cells);
    E_coeff = reshape(E_coeff_vec, p_order+1, num_cells);
end
