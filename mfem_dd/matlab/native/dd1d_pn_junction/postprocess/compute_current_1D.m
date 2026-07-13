function current_data = compute_current_1D(n_coeff, E_coeff, mesh, basis_data, params)
%COMPUTE_CURRENT_1D Compute J = -D n_x + mu n E cellwise and at contacts.

    n_coeff = reshape(n_coeff, basis_data.nloc, mesh.num_cells);
    E_coeff = reshape(E_coeff, basis_data.nloc, mesh.num_cells);

    D = params.problem.physics.diffusion;
    mu = params.problem.physics.mobility;

    basis_values = basis_data.basis_values;
    basis_grad_values = basis_data.basis_grad_values;
    face_values = basis_data.face_values;
    face_grad_values = basis_data.face_grad_values;

    n_val = basis_values * n_coeff;
    E_val = basis_values * E_coeff;
    n_x = (basis_grad_values * n_coeff) ./ mesh.cell_sizes;

    J = -D * n_x + mu * (n_val .* E_val);

    current_data.gauss_current = J;
    current_data.cell_mean = mean(J, 1).';

    n_left  = face_values(1,:) * n_coeff(:,1);
    n_right = face_values(2,:) * n_coeff(:,end);

    E_left  = face_values(1,:) * E_coeff(:,1);
    E_right = face_values(2,:) * E_coeff(:,end);

    n_x_left  = (face_grad_values(1,:) * n_coeff(:,1)) / mesh.cell_sizes(1);
    n_x_right = (face_grad_values(2,:) * n_coeff(:,end)) / mesh.cell_sizes(end);

    current_data.left_contact  = -D * n_x_left  + mu * n_left  * E_left;
    current_data.right_contact = -D * n_x_right + mu * n_right * E_right;
    current_data.domain_average = mean(current_data.cell_mean);
end