function transport_rhs = assemble_transport_rhs_single(mesh, basis_data, carrier_coeff, E_coeff, time_now, params)
%ASSEMBLE_TRANSPORT_RHS_SINGLE Assemble the explicit transport RHS.

    num_cells = mesh.num_cells;
    cell_sizes = mesh.cell_sizes;
    cell_centers = mesh.cell_centers;
    problem = params.problem;

    p_order = basis_data.p_order;
    Gp = basis_data.gauss_points;
    Gw = basis_data.gauss_weights;
    basis_values = basis_data.basis_values;
    basis_grad_values = basis_data.basis_grad_values;
    face_values = basis_data.face_values;

    mobility = problem.physics.mobility;
    emf = mobility;

    carrier_coeff = reshape(carrier_coeff, p_order+1, num_cells);
    E_coeff = reshape(E_coeff, p_order+1, num_cells);

    carrier_values = basis_values * carrier_coeff;
    field_values = basis_values * E_coeff;
    flux_values = -(carrier_values .* field_values);
    volume_term = -(basis_grad_values' * bsxfun(@times, Gw, flux_values));

    carrier_left = face_values(1,:) * carrier_coeff;
    carrier_right = face_values(2,:) * carrier_coeff;
    field_left = face_values(1,:) * E_coeff;
    field_right = face_values(2,:) * E_coeff;

    flux_left = -(carrier_left .* field_left);
    flux_right = -(carrier_right .* field_right);

    right_term = zeros(1, num_cells);
    left_term = zeros(1, num_cells);

    if isfield(problem, 'transport_bc_type') && strcmpi(problem.transport_bc_type, 'dirichlet')
        % interior interfaces
        for cell_id = 1:num_cells-1
            right_term(cell_id) = 0.5 * (flux_right(cell_id) + flux_left(cell_id+1)) ...
                - 0.5 * emf * (carrier_left(cell_id+1) - carrier_right(cell_id));
            left_term(cell_id+1) = -0.5 * (flux_right(cell_id) + flux_left(cell_id+1)) ...
                + 0.5 * emf * (carrier_left(cell_id+1) - carrier_right(cell_id));
        end

        % left boundary
        carrier_bc_left = problem.boundary.transport_left(time_now);
        flux_bc_left = -(carrier_bc_left .* field_left(1));
        left_term(1) = -0.5 * (flux_bc_left + flux_left(1)) ...
            + 0.5 * emf * (carrier_bc_left - carrier_left(1));

        % right boundary
        carrier_bc_right = problem.boundary.transport_right(time_now);
        flux_bc_right = -(carrier_bc_right .* field_right(end));
        right_term(end) = 0.5 * (flux_right(end) + flux_bc_right) ...
            - 0.5 * emf * (carrier_bc_right - carrier_right(end));
    else
        next_idx = [2:num_cells, 1];
        prev_idx = [num_cells, 1:num_cells-1];

        right_term = 0.5 * (flux_right + flux_left(next_idx)) ...
            - 0.5 * emf * (carrier_left(next_idx) - carrier_right);
        left_term = -0.5 * (flux_right(prev_idx) + flux_left) ...
            + 0.5 * emf * (carrier_left - carrier_right(prev_idx));
    end

    interface_term = face_values(2,:).' * right_term + face_values(1,:).' * left_term;

    x = Gp * cell_sizes + cell_centers;
    source_values = problem.source.transport(x, time_now);
    source_term = basis_values' * bsxfun(@times, Gw, source_values);
    source_term = bsxfun(@times, source_term, cell_sizes);

    transport_rhs = emf * (volume_term + interface_term) + source_term;
    transport_rhs = transport_rhs(:);
end
