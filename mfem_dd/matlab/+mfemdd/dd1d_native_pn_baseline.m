function baseline = dd1d_native_pn_baseline(varargin)
%DD1D_NATIVE_PN_BASELINE Recompute the legacy 1D PN junction experiment.
opts = local_parse_options(varargin{:});
native_root = local_native_root();

old_path = path;
cleanup = onCleanup(@() path(old_path));
addpath(genpath(native_root));

[params, basis_data, mesh] = local_build_problem(opts);

params_eq = params;
params_eq.device.bias = 0.0;
params_eq.problem = build_problem_definition(params_eq);
result_eq = solve_1D_DD(mesh, basis_data, params_eq);

bias_list = linspace(-1.0, 1.0, 9);
iv_data = zeros(numel(bias_list), 4);
prev_coeff = [];
iv_results = cell(numel(bias_list), 1);

for k = 1:numel(bias_list)
    params_iv = params;
    params_iv.device.bias = bias_list(k);
    params_iv.problem = build_problem_definition(params_iv);
    if ~isempty(prev_coeff)
        params_iv.restart.initial_n_coeff = prev_coeff(:);
    end

    result_iv = solve_1D_DD(mesh, basis_data, params_iv);
    charge_data = compute_space_charge_proxy_1D(result_iv.n_coeff, ...
        mesh, basis_data, params_iv);

    iv_data(k, :) = [bias_list(k), result_iv.current.right_contact, ...
        charge_data.Qmag, charge_data.Wdep];
    prev_coeff = result_iv.n_coeff(:);
    iv_results{k} = result_iv;
end

cv_data = local_cv_from_iv(iv_data);

params_fwd = params;
params_fwd.device.bias = 0.6;
params_fwd.time.T_end = 0.3;
params_fwd.problem = build_problem_definition(params_fwd);
res_fwd = solve_1D_DD(mesh, basis_data, params_fwd);

params_trans = params;
params_trans.device.bias = -0.8;
params_trans.time.T_end = 0.2;
params_trans.output.store_history = true;
params_trans.output.snapshot_times = [0.02, 0.05, 0.10, 0.20];
params_trans.restart.initial_n_coeff = res_fwd.n_coeff(:);
params_trans.problem = build_problem_definition(params_trans);
res_trans = solve_1D_DD(mesh, basis_data, params_trans);

transient_data = [res_trans.history.time, res_trans.history.current];

iv_raw = array2table(iv_data, "VariableNames", ...
    {'bias', 'right_current', 'Qmag', 'Wdep'});
cv_raw = array2table(cv_data, "VariableNames", {'bias', 'Qmag', 'Cqs'});
transient_raw = array2table(transient_data, "VariableNames", ...
    {'time', 'right_current'});

iv = local_legacy_csv_table(iv_raw);
cv = local_legacy_csv_table(cv_raw);
transient = local_legacy_csv_table(transient_raw);

summary = local_summary(iv, cv, transient);

baseline = struct();
baseline.case_name = "dd_pn_device";
baseline.device = "pn_junction";
baseline.generated_by = "dd1d_native_pn_baseline";
baseline.source = "mfem_dd_native_matlab";
baseline.native_root = string(native_root);
baseline.summary = summary;
baseline.iv = iv;
baseline.cv = cv;
baseline.transient = transient;
baseline.iv_raw = iv_raw;
baseline.cv_raw = cv_raw;
baseline.transient_raw = transient_raw;
baseline.output_format = "legacy_dlmwrite_precision";
baseline.params = params;
baseline.mesh = mesh;
baseline.equilibrium_result = result_eq;
baseline.iv_results = iv_results;
baseline.forward_result = res_fwd;
baseline.transient_result = res_trans;
end

function [params, basis_data, mesh] = local_build_problem(opts)
params = build_pn_junction_params();
params.ipdg.method = char(opts.method);
params.ipdg = build_method_flags(params.ipdg);
params.problem = build_problem_definition(params);
params.output.verbose = opts.verbose;
params.output.compute_Linf = false;
params.output.store_history = false;
params.output.snapshot_times = [];

basis_data = build_basis_data_1D(params);
mesh = build_mesh_1D_uniform(params, params.mesh.ng_base);
end

function cv_data = local_cv_from_iv(iv_data)
if size(iv_data, 1) >= 3
    cv_data = zeros(size(iv_data, 1) - 2, 3);
    for k = 2:size(iv_data, 1) - 1
        dQdV = (iv_data(k + 1, 3) - iv_data(k - 1, 3)) / ...
            (iv_data(k + 1, 1) - iv_data(k - 1, 1));
        cv_data(k - 1, :) = [iv_data(k, 1), iv_data(k, 3), dQdV];
    end
else
    cv_data = [iv_data(1, 1), iv_data(1, 3), NaN];
end
end

function rounded = local_legacy_csv_table(raw)
rounded = raw;
names = string(raw.Properties.VariableNames);
for k = 1:numel(names)
    name = names(k);
    rounded.(name) = arrayfun(@local_legacy_csv_number, raw.(name));
end
end

function value = local_legacy_csv_number(value)
if isnan(value) || isinf(value)
    return;
end
value = str2double(sprintf("%.5g", value));
end

function summary = local_summary(iv, cv, transient)
summary = struct();
summary.iv_rows = height(iv);
summary.cv_rows = height(cv);
summary.transient_rows = height(transient);
summary.iv_reverse_current_minus1v = local_lookup(iv, "bias", -1.0, "right_current");
summary.iv_zero_bias_current = local_lookup(iv, "bias", 0.0, "right_current");
summary.iv_forward_current_1v = local_lookup(iv, "bias", 1.0, "right_current");
summary.iv_zero_bias_qmag = local_lookup(iv, "bias", 0.0, "Qmag");
summary.cv_zero_bias_cqs = local_lookup(cv, "bias", 0.0, "Cqs");

finite_id = find(isfinite(transient.right_current), 1, "first");
summary.transient_first_finite_time = transient.time(finite_id);
summary.transient_first_finite_current = transient.right_current(finite_id);
summary.transient_terminal_time = transient.time(end);
summary.transient_terminal_current = transient.right_current(end);
end

function value = local_lookup(tbl, key_name, key_value, value_name)
keys = tbl.(key_name);
idx = find(abs(keys - key_value) <= 1.0e-12, 1, "first");
if isempty(idx)
    error("Could not find %s=%g in native PN table.", key_name, key_value);
end
value = tbl.(value_name)(idx);
end

function native_root = local_native_root()
root = fileparts(fileparts(mfilename("fullpath")));
native_root = fullfile(root, "native", "dd1d_pn_junction");
if ~exist(native_root, "dir")
    error("DD1D PN native source mirror is missing: %s", native_root);
end
end

function opts = local_parse_options(varargin)
opts = struct("method", "SIPG", "verbose", false);
if mod(numel(varargin), 2) ~= 0
    error("Options must be name/value pairs.");
end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k + 1};
    switch name
        case "method"
            opts.method = string(value);
        case "verbose"
            opts.verbose = logical(value);
        otherwise
            error("Unknown option: %s", name);
    end
end
end
