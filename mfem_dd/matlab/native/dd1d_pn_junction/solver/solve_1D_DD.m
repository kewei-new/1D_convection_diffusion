function result = solve_1D_DD(mesh, basis_data, params)
%SOLVE_1D_DD Solve the 1D drift-diffusion problem.

    num_cells = mesh.num_cells;
    nloc = basis_data.nloc;
    total_dofs = num_cells * nloc;
    T_end = params.time.T_end;
    problem = params.problem;

    % -------------------------------------------------
    % 1. Initial projection / restart
    % -------------------------------------------------
    if isfield(params, 'restart') && isfield(params.restart, 'initial_n_coeff') ...
            && ~isempty(params.restart.initial_n_coeff)
        n_old = params.restart.initial_n_coeff(:);
    else
        n_old = L2_projection_1D(problem.initial.n, basis_data, mesh);
    end

    % -------------------------------------------------
    % 2. Build time-step size and matrices
    % -------------------------------------------------
    dt_base = evaluate_dt(mesh, params);

    diffusion_matrix = problem.physics.diffusion * ...
        assemble_ipdg_diffusion_matrix(mesh, basis_data, params);

    inv_mass = basis_data.inv_mass(:);
    mass_diag = reshape((1 ./ inv_mass) * mesh.cell_sizes, total_dofs, 1);

    [poisson_matrix_main, poisson_matrix_aux, poisson_rhs_matrix] = ...
        assemble_LDG_matrix(basis_data.gauss_data, basis_data.inv_mass, ...
            basis_data.p_order, num_cells, mesh.cell_sizes);

    poisson_data.main = poisson_matrix_main;
    poisson_data.aux = poisson_matrix_aux;
    poisson_data.rhs = poisson_rhs_matrix;
    poisson_data.solver = build_linear_solver(poisson_rhs_matrix);

    % -------------------------------------------------
    % 3. Time evolution
    % -------------------------------------------------
    store_history = isfield(params.output, 'store_history') && params.output.store_history;
    if store_history
        history_time = 0;
        history_current = NaN;
    else
        history_time = [];
        history_current = [];
    end

    snapshot_times = [];
    if isfield(params.output, 'snapshot_times') && ~isempty(params.output.snapshot_times)
        snapshot_times = params.output.snapshot_times(:).';
    end
    snapshot_coeff = cell(size(snapshot_times));
    next_snapshot_id = 1;

    time_now = 0;
    current_dt = NaN;
    current_mass_dt_diag = [];
    implicit_solver = [];

    while time_now < T_end
        dt = min(dt_base, T_end - time_now);

        if isempty(current_mass_dt_diag) || abs(dt - current_dt) > max(1e-14, eps(max(dt,1)))
            current_dt = dt;
            current_mass_dt_diag = mass_diag / dt;
            implicit_matrix = spdiags(current_mass_dt_diag, 0, total_dofs, total_dofs) ...
                - 0.5 * diffusion_matrix;
            implicit_solver = build_linear_solver(implicit_matrix);
        end

        mass_term = current_mass_dt_diag .* n_old;

        % stage 1
        [~, E_coeff] = solve_1D_LDG(mesh, basis_data, n_old, time_now, params, poisson_data);
        transport_rhs_0 = assemble_transport_rhs_single(mesh, basis_data, n_old, E_coeff, time_now, params);
        n_stage1 = solve_with_linear_solver(implicit_solver, 0.5 * transport_rhs_0 + mass_term);

        % stage 2
        [~, E_coeff] = solve_1D_LDG(mesh, basis_data, n_stage1, time_now + dt/2, params, poisson_data);
        transport_rhs_1 = assemble_transport_rhs_single(mesh, basis_data, n_stage1, E_coeff, time_now + dt/2, params);
        n_stage2 = solve_with_linear_solver(implicit_solver, ...
            11/18 * transport_rhs_0 + 1/18 * transport_rhs_1 + ...
            1/6 * diffusion_matrix * n_stage1 + mass_term);

        % stage 3
        [~, E_coeff] = solve_1D_LDG(mesh, basis_data, n_stage2, time_now + 2*dt/3, params, poisson_data);
        transport_rhs_2 = assemble_transport_rhs_single(mesh, basis_data, n_stage2, E_coeff, time_now + 2*dt/3, params);
        n_stage3 = solve_with_linear_solver(implicit_solver, ...
            5/6 * transport_rhs_0 - 5/6 * transport_rhs_1 + 1/2 * transport_rhs_2 ...
            - 1/2 * diffusion_matrix * n_stage1 + 1/2 * diffusion_matrix * n_stage2 + mass_term);

        % stage 4
        [~, E_coeff] = solve_1D_LDG(mesh, basis_data, n_stage3, time_now + dt/2, params, poisson_data);
        transport_rhs_3 = assemble_transport_rhs_single(mesh, basis_data, n_stage3, E_coeff, time_now + dt/2, params);
        n_new = solve_with_linear_solver(implicit_solver, ...
            1/4 * transport_rhs_0 + 7/4 * transport_rhs_1 + 3/4 * transport_rhs_2 - 7/4 * transport_rhs_3 ...
            + 3/2 * diffusion_matrix * n_stage1 - 3/2 * diffusion_matrix * n_stage2 + 1/2 * diffusion_matrix * n_stage3 ...
            + mass_term);

        time_new = time_now + dt;
        n_old = n_new;
        time_now = time_new;

        [phi_coeff_tmp, E_coeff_tmp] = solve_1D_LDG(mesh, basis_data, n_old, time_now, params, poisson_data);
        if store_history
            current_data = compute_current_1D(n_old, E_coeff_tmp, mesh, basis_data, params);
            history_time(end+1,1) = time_now; %#ok<AGROW>
            history_current(end+1,1) = current_data.right_contact; %#ok<AGROW>
        end
        % if store_history
        %     history_time(end+1,1) = time_now;
        %     history_current(end+1,1) = NaN;
        % end

        while next_snapshot_id <= numel(snapshot_times) && time_now >= snapshot_times(next_snapshot_id) - 1e-14
            snapshot_coeff{next_snapshot_id} = struct(...
                'time', time_now, ...
                'n_coeff', reshape(n_old, nloc, num_cells), ...
                'phi_coeff', phi_coeff_tmp, ...
                'E_coeff', E_coeff_tmp);
            next_snapshot_id = next_snapshot_id + 1;
        end
    end

    [phi_coeff, E_coeff] = solve_1D_LDG(mesh, basis_data, n_old, time_now, params, poisson_data);
    current_data = compute_current_1D(n_old, E_coeff, mesh, basis_data, params);
    % current_data = struct();

    result.n_coeff   = reshape(n_old, nloc, num_cells);
    result.phi_coeff = phi_coeff;
    result.E_coeff   = E_coeff;
    result.time      = time_now;
    result.dt        = current_dt;
    result.current   = current_data;
    result.history.time = history_time;
    result.history.current = history_current;
    result.snapshots = snapshot_coeff;
end
