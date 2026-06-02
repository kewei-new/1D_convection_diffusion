function snapshot = dd1d_native_pn_operator_snapshot(varargin)
%DD1D_NATIVE_PN_OPERATOR_SNAPSHOT Export compact PN device PDE evidence.
%
% This snapshot is a reference target for the native C++ PN-device port. It
% records operator, vector, first-step, and current summaries without storing
% large dense arrays.
opts = local_parse_options(varargin{:});
native_root = local_native_root();

old_path = path;
cleanup = onCleanup(@() path(old_path));
addpath(genpath(native_root));

params = build_pn_junction_params();
params.mesh.ng_base = opts.elements;
params.ipdg.method = char(opts.method);
params.ipdg = build_method_flags(params.ipdg);
params.output.verbose = opts.verbose;
params.output.compute_Linf = false;
params.output.store_history = false;
params.output.snapshot_times = [];
params.problem = build_problem_definition(params);

basis_data = build_basis_data_1D(params);
mesh = build_mesh_1D_uniform(params, params.mesh.ng_base);

n0 = L2_projection_1D(params.problem.initial.n, basis_data, mesh);
dt = evaluate_dt(mesh, params);
diffusion_matrix = params.problem.physics.diffusion * ...
    assemble_ipdg_diffusion_matrix(mesh, basis_data, params);
mass_diag = reshape((1 ./ basis_data.inv_mass(:)) * mesh.cell_sizes, ...
    mesh.num_cells * basis_data.nloc, 1);

[poisson_main, poisson_aux, poisson_rhs] = assemble_LDG_matrix( ...
    basis_data.gauss_data, basis_data.inv_mass, basis_data.p_order, ...
    mesh.num_cells, mesh.cell_sizes);
poisson_data = struct("main", poisson_main, "aux", poisson_aux, ...
    "rhs", poisson_rhs, "solver", build_linear_solver(poisson_rhs));

[phi0, E0] = solve_1D_LDG(mesh, basis_data, n0, 0.0, params, poisson_data);
rhs0 = assemble_transport_rhs_single(mesh, basis_data, n0, E0, 0.0, params);
current0 = compute_current_1D(n0, E0, mesh, basis_data, params);

implicit_matrix = spdiags(mass_diag / dt, 0, numel(n0), numel(n0)) ...
    - 0.5 * diffusion_matrix;
implicit_solver = build_linear_solver(implicit_matrix);
mass_term = (mass_diag / dt) .* n0;

n1 = solve_with_linear_solver(implicit_solver, 0.5 * rhs0 + mass_term);
[phi1, E1] = solve_1D_LDG(mesh, basis_data, n1, dt/2, params, poisson_data);
rhs1 = assemble_transport_rhs_single(mesh, basis_data, n1, E1, dt/2, params);

n2 = solve_with_linear_solver(implicit_solver, ...
    11/18 * rhs0 + 1/18 * rhs1 + 1/6 * diffusion_matrix * n1 + mass_term);
[phi2, E2] = solve_1D_LDG(mesh, basis_data, n2, 2*dt/3, params, poisson_data);
rhs2 = assemble_transport_rhs_single(mesh, basis_data, n2, E2, 2*dt/3, params);

n3 = solve_with_linear_solver(implicit_solver, ...
    5/6 * rhs0 - 5/6 * rhs1 + 1/2 * rhs2 ...
    - 1/2 * diffusion_matrix * n1 + 1/2 * diffusion_matrix * n2 + mass_term);
[phi3, E3] = solve_1D_LDG(mesh, basis_data, n3, dt/2, params, poisson_data);
rhs3 = assemble_transport_rhs_single(mesh, basis_data, n3, E3, dt/2, params);

n_step = solve_with_linear_solver(implicit_solver, ...
    1/4 * rhs0 + 7/4 * rhs1 + 3/4 * rhs2 - 7/4 * rhs3 ...
    + 3/2 * diffusion_matrix * n1 - 3/2 * diffusion_matrix * n2 ...
    + 1/2 * diffusion_matrix * n3 + mass_term);
[phi_step, E_step] = solve_1D_LDG(mesh, basis_data, n_step, dt, params, poisson_data);
current_step = compute_current_1D(n_step, E_step, mesh, basis_data, params);

snapshot = struct();
snapshot.case_name = "dd_pn_device";
snapshot.backend = "native_matlab_pn_operator_snapshot";
snapshot.generated_by = "dd1d_native_pn_operator_snapshot";
snapshot.native_root = string(native_root);
snapshot.method = string(params.ipdg.method);
snapshot.alpha = params.ipdg.alpha;
snapshot.beta = params.ipdg.beta;
snapshot.mesh = struct("elements", mesh.num_cells, "left", mesh.left, ...
    "right", mesh.right, "h", mesh.cell_sizes(1));
snapshot.space = struct("degree", basis_data.p_order, ...
    "local_dofs", basis_data.nloc, "total_dofs", numel(n0));
snapshot.time = struct("dt", dt, "stage1_time", dt/2, ...
    "stage2_time", 2*dt/3, "stage3_time", dt/2, "step_time", dt);
snapshot.physics = struct("mobility", params.problem.physics.mobility, ...
    "diffusion", params.problem.physics.diffusion, ...
    "poisson_scale", params.problem.physics.poisson_scale, ...
    "bias", params.device.bias);

snapshot.operators = struct();
snapshot.operators.diffusion = local_matrix_stats(diffusion_matrix);
snapshot.operators.poisson_main = local_matrix_stats(poisson_main);
snapshot.operators.poisson_aux = local_matrix_stats(poisson_aux);
snapshot.operators.poisson_rhs = local_matrix_stats(poisson_rhs);
snapshot.operators.implicit = local_matrix_stats(implicit_matrix);

snapshot.vectors = struct();
snapshot.vectors.mass_diag = local_vector_stats(mass_diag);
snapshot.vectors.n0 = local_vector_stats(n0);
snapshot.vectors.phi0 = local_vector_stats(phi0);
snapshot.vectors.E0 = local_vector_stats(E0);
snapshot.vectors.rhs0 = local_vector_stats(rhs0);
snapshot.vectors.n1 = local_vector_stats(n1);
snapshot.vectors.phi1 = local_vector_stats(phi1);
snapshot.vectors.E1 = local_vector_stats(E1);
snapshot.vectors.rhs1 = local_vector_stats(rhs1);
snapshot.vectors.n2 = local_vector_stats(n2);
snapshot.vectors.phi2 = local_vector_stats(phi2);
snapshot.vectors.E2 = local_vector_stats(E2);
snapshot.vectors.rhs2 = local_vector_stats(rhs2);
snapshot.vectors.n3 = local_vector_stats(n3);
snapshot.vectors.phi3 = local_vector_stats(phi3);
snapshot.vectors.E3 = local_vector_stats(E3);
snapshot.vectors.rhs3 = local_vector_stats(rhs3);
snapshot.vectors.n_step = local_vector_stats(n_step);
snapshot.vectors.phi_step = local_vector_stats(phi_step);
snapshot.vectors.E_step = local_vector_stats(E_step);
snapshot.current = struct("initial_right_contact", current0.right_contact, ...
    "initial_left_contact", current0.left_contact, ...
    "initial_domain_average", current0.domain_average, ...
    "step_right_contact", current_step.right_contact, ...
    "step_left_contact", current_step.left_contact, ...
    "step_domain_average", current_step.domain_average);

if strlength(opts.output_json) > 0
    parent = fileparts(opts.output_json);
    if strlength(parent) > 0 && ~exist(parent, "dir")
        mkdir(parent);
    end
    fid = fopen(opts.output_json, "w");
    if fid < 0
        error("dd1d_native_pn_operator_snapshot:cannotOpen", ...
            "Cannot open output JSON: %s", opts.output_json);
    end
    cleanup_file = onCleanup(@() fclose(fid));
    fprintf(fid, "%s", jsonencode(snapshot, PrettyPrint=true));
end
end

function stats = local_matrix_stats(A)
[i, j, v] = find(A);
stats = struct();
stats.rows = size(A, 1);
stats.cols = size(A, 2);
stats.nnz = nnz(A);
stats.fro_norm = norm(A, "fro");
stats.one_norm = norm(A, 1);
stats.inf_norm = norm(A, inf);
stats.value_sum = full(sum(v));
stats.abs_sum = full(sum(abs(v)));
stats.weighted_sum = full(sum((double(i) + 0.125 * double(j)) .* v));
end

function stats = local_vector_stats(x)
v = x(:);
weights = (1:numel(v)).';
stats = struct();
stats.size = numel(v);
stats.min = min(v);
stats.max = max(v);
stats.mean = mean(v);
stats.norm2 = norm(v);
stats.value_sum = sum(v);
stats.abs_sum = sum(abs(v));
stats.weighted_sum = sum(weights .* v);
end

function native_root = local_native_root()
root = fileparts(fileparts(mfilename("fullpath")));
native_root = fullfile(root, "native", "dd1d_pn_junction");
if ~exist(native_root, "dir")
    error("DD1D PN native source mirror is missing: %s", native_root);
end
end

function opts = local_parse_options(varargin)
opts = struct("elements", 160, "method", "SIPG", ...
    "verbose", false, "output_json", "");
if mod(numel(varargin), 2) ~= 0
    error("Options must be name/value pairs.");
end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k + 1};
    switch name
        case {"elements", "ng_base"}
            opts.elements = value;
        case "method"
            opts.method = string(value);
        case "verbose"
            opts.verbose = logical(value);
        case "output_json"
            opts.output_json = string(value);
        otherwise
            error("Unknown option: %s", name);
    end
end
end
