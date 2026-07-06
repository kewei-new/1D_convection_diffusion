function baseline = legacy_dd1d_baseline(varargin)
%LEGACY_DD1D_BASELINE Run or read the legacy DD1D smooth MMS baseline.
opts = local_parse_options(varargin{:});
root = fileparts(fileparts(fileparts(mfilename("fullpath"))));
repo_root = fileparts(root);
baseline_root = fullfile(root, "cases", "dd1d_smooth_mms", "legacy_baseline");
result_file = fullfile(baseline_root, ...
    sprintf("error_table_%s_p%d.txt", upper(opts.method), opts.p_order));

if opts.run_legacy
    legacy_root = fullfile(repo_root, "DD验阶", "DD1D_smooth", "V1");
    if ~isfolder(legacy_root)
        error("Legacy DD1D runtime folder has been removed after migration: %s", legacy_root);
    end
    baseline = local_run_legacy(legacy_root, opts);
else
    baseline = mfemdd.read_legacy_error_table(result_file);
    baseline.source = string(result_file);
    baseline.generated_by = "controlled_legacy_table";
    baseline.original_legacy_source = "DD验阶/DD1D_smooth/V1/result/error_table_SIPG_p3.txt";
end
end

function baseline = local_run_legacy(legacy_root, opts)
old_path = path;
cleanup = onCleanup(@() path(old_path));
addpath(genpath(legacy_root));

params = build_default_params();
params.ipdg.method = char(opts.method);
params.dg.p_order = opts.p_order;
params.dg.nloc = opts.p_order + 1;
params.mesh.refine_steps = opts.refine_steps;
params.mesh.ng_base = opts.ng_base;
params.output.verbose = opts.verbose;
params.ipdg = build_method_flags(params.ipdg);
params.ipdg.alpha = (params.dg.p_order + 1)^2;

basis_data = build_basis_data_1D(params);
errors = struct();
errors.n_L2 = zeros(params.mesh.refine_steps, 1);
errors.n_Linf = zeros(params.mesh.refine_steps, 1);
errors.phi_L2 = zeros(params.mesh.refine_steps, 1);
errors.phi_Linf = zeros(params.mesh.refine_steps, 1);
errors.E_L2 = zeros(params.mesh.refine_steps, 1);
errors.E_Linf = zeros(params.mesh.refine_steps, 1);
mesh_sizes = zeros(params.mesh.refine_steps, 1);

for refine_id = 1:params.mesh.refine_steps
    num_cells = params.mesh.ng_base * 2^(refine_id - 1);
    mesh = build_mesh_1D_uniform(params, num_cells);
    mesh_sizes(refine_id) = mesh.cell_sizes(1);
    result = solve_1D_DD(mesh, basis_data, params);
    errors.n_L2(refine_id) = evaluate_error_Lnorm(result.n_coeff, "exact_fun_n", 2, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.n_Linf(refine_id) = evaluate_error_absolute(result.n_coeff, "exact_fun_n", ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.phi_L2(refine_id) = evaluate_error_Lnorm(result.phi_coeff, "exact_fun_phi", 2, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.phi_Linf(refine_id) = evaluate_error_absolute(result.phi_coeff, "exact_fun_phi", ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.E_L2(refine_id) = evaluate_error_Lnorm(result.E_coeff, "exact_fun_phi_x", 2, ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
    errors.E_Linf(refine_id) = evaluate_error_absolute(result.E_coeff, "exact_fun_phi_x", ...
        mesh.num_cells, basis_data.p_order, mesh.cell_sizes, mesh.cell_centers, result.time);
end

baseline = struct();
baseline.method = string(params.ipdg.method);
baseline.p_order = params.dg.p_order;
baseline.alpha = params.ipdg.alpha;
baseline.beta = params.ipdg.beta;
baseline.columns = ["h", "n_L2", "n_L2_order", "n_Linf", "n_Linf_order", ...
    "phi_L2", "phi_L2_order", "phi_Linf", "phi_Linf_order", ...
    "E_L2", "E_L2_order", "E_Linf", "E_Linf_order"];
baseline.data = [mesh_sizes, errors.n_L2, local_orders(errors.n_L2), ...
    errors.n_Linf, local_orders(errors.n_Linf), ...
    errors.phi_L2, local_orders(errors.phi_L2), ...
    errors.phi_Linf, local_orders(errors.phi_Linf), ...
    errors.E_L2, local_orders(errors.E_L2), ...
    errors.E_Linf, local_orders(errors.E_Linf)];
baseline.source = "computed_legacy_adapter";
baseline.generated_by = "legacy_dd1d_baseline";
end

function orders = local_orders(errors)
orders = nan(size(errors));
for k = 2:numel(errors)
    if errors(k) > 0 && errors(k-1) > 0
        orders(k) = log2(errors(k-1) / errors(k));
    end
end
end

function opts = local_parse_options(varargin)
opts = struct("method", "SIPG", "p_order", 3, "run_legacy", false, ...
    "refine_steps", 5, "ng_base", 20, "verbose", false);
if mod(numel(varargin), 2) ~= 0
    error("Options must be name/value pairs.");
end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k+1};
    switch name
        case "method"
            opts.method = string(value);
        case "p_order"
            opts.p_order = value;
        case "run_legacy"
            opts.run_legacy = logical(value);
        case "refine_steps"
            opts.refine_steps = value;
        case "ng_base"
            opts.ng_base = value;
        case "verbose"
            opts.verbose = logical(value);
        otherwise
            error("Unknown option: %s", name);
    end
end
end
