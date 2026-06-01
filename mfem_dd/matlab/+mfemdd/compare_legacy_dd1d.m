function report = compare_legacy_dd1d(varargin)
%COMPARE_LEGACY_DD1D Compare current MFEM-DD output with the legacy table.
opts = local_parse_options(varargin{:});
baseline = mfemdd.legacy_dd1d_baseline("method", opts.method, "p_order", opts.p_order);
current_table = mfemdd.dd1d_convergence_table("method", opts.method, ...
    "p_order", opts.p_order, ...
    "backend", opts.backend, ...
    "refine_steps", opts.refine_steps, ...
    "ng_base", opts.ng_base, ...
    "verbose", opts.verbose);
comparison_rows = min(size(current_table.data, 1), size(baseline.data, 1));
baseline_data = baseline.data(1:comparison_rows, :);
current_data = current_table.data(1:comparison_rows, :);
diff_data = abs(current_data - baseline_data);
diff_data(isnan(current_data) & isnan(baseline_data)) = 0.0;
diff_data(xor(isnan(current_data), isnan(baseline_data))) = Inf;
max_abs_diff = max(diff_data, [], 1, "omitnan");
max_rel_diff = local_max_rel_diff(baseline_data, current_data);
exact_compared_rows = all(max_abs_diff <= opts.abs_tol | isnan(max_abs_diff));
tolerance_pass = local_tolerance_pass(baseline.columns, baseline_data, current_data, opts);

legacy_first = local_row_struct(baseline.columns, baseline.data(1,:));
current_first = local_row_struct(current_table.columns, current_table.data(1,:));

report = struct();
report.case_name = "dd1d_smooth_mms";
report.method = baseline.method;
report.p_order = baseline.p_order;
report.legacy_source = baseline.source;
report.backend = opts.backend;
report.comparison_rows = comparison_rows;
report.legacy_rows = size(baseline.data, 1);
report.current_rows = size(current_table.data, 1);
report.legacy_first = legacy_first;
report.current_first = current_first;
report.columns = baseline.columns;
report.max_abs_diff = max_abs_diff;
report.max_abs_diff_by_column = local_diff_struct(baseline.columns, max_abs_diff);
report.max_rel_diff = max_rel_diff;
report.max_rel_diff_by_column = local_diff_struct(baseline.columns, max_rel_diff);
report.matches_compared_rows_exact = exact_compared_rows;
report.matches_full_table_exact = exact_compared_rows ...
    && report.comparison_rows == report.legacy_rows ...
    && report.current_rows == report.legacy_rows;
report.matches_compared_rows = tolerance_pass;
report.matches_full_table = report.matches_compared_rows ...
    && report.comparison_rows == report.legacy_rows ...
    && report.current_rows == report.legacy_rows;
report.tolerances = struct("h_abs_tol", opts.h_abs_tol, ...
    "error_rel_tol", opts.error_rel_tol, ...
    "order_abs_tol", opts.order_abs_tol, ...
    "abs_tol", opts.abs_tol);
report.legacy_table = baseline;
report.current_table = current_table;

if strlength(opts.output_json) > 0
    parent = fileparts(opts.output_json);
    if strlength(parent) > 0 && ~exist(parent, "dir")
        mkdir(parent);
    end
    fid = fopen(opts.output_json, "w");
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, "%s", jsonencode(report, PrettyPrint=true));
end
end

function row = local_row_struct(columns, values)
row = struct();
for k = 1:numel(columns)
    row.(char(columns(k))) = values(k);
end
end

function row = local_diff_struct(columns, values)
row = struct();
for k = 1:numel(columns)
    row.(char(columns(k))) = values(k);
end
end

function max_rel_diff = local_max_rel_diff(baseline_data, current_data)
rel_data = abs(current_data - baseline_data) ./ max(abs(baseline_data), realmin);
rel_data(isnan(current_data) & isnan(baseline_data)) = 0.0;
rel_data(xor(isnan(current_data), isnan(baseline_data))) = Inf;
max_rel_diff = max(rel_data, [], 1, "omitnan");
end

function pass = local_tolerance_pass(columns, baseline_data, current_data, opts)
pass = true;
for column_id = 1:numel(columns)
    name = columns(column_id);
    b = baseline_data(:, column_id);
    c = current_data(:, column_id);
    both_nan = isnan(b) & isnan(c);
    one_nan = xor(isnan(b), isnan(c));
    abs_diff = abs(c - b);

    if name == "h"
        column_pass = abs_diff <= opts.h_abs_tol;
    elseif contains(name, "_order")
        column_pass = abs_diff <= opts.order_abs_tol;
    else
        column_pass = abs_diff <= max(opts.abs_tol, opts.error_rel_tol .* max(abs(b), realmin));
    end

    column_pass(both_nan) = true;
    column_pass(one_nan) = false;
    pass = pass && all(column_pass);
end
end

function opts = local_parse_options(varargin)
opts = struct("method", "SIPG", "p_order", 3, "abs_tol", 1.0e-12, ...
    "output_json", "", "backend", "matlab_mfem", ...
    "refine_steps", 5, "ng_base", 20, "verbose", false, ...
    "h_abs_tol", 1.0e-6, "error_rel_tol", 1.0e-6, ...
    "order_abs_tol", 0.15);
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
        case "abs_tol"
            opts.abs_tol = value;
        case "output_json"
            opts.output_json = string(value);
        case "backend"
            opts.backend = string(value);
        case "refine_steps"
            opts.refine_steps = value;
        case "ng_base"
            opts.ng_base = value;
        case "verbose"
            opts.verbose = logical(value);
        case "h_abs_tol"
            opts.h_abs_tol = value;
        case "error_rel_tol"
            opts.error_rel_tol = value;
        case "order_abs_tol"
            opts.order_abs_tol = value;
        otherwise
            error("Unknown option: %s", name);
    end
end
end
