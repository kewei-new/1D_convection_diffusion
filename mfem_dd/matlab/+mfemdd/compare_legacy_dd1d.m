function report = compare_legacy_dd1d(varargin)
%COMPARE_LEGACY_DD1D Compare current MFEM-DD output with the legacy table.
opts = local_parse_options(varargin{:});
baseline = mfemdd.legacy_dd1d_baseline("method", opts.method, "p_order", opts.p_order);
current_table = mfemdd.dd1d_convergence_table("method", opts.method, ...
    "p_order", opts.p_order, "backend", opts.backend);
diff_data = abs(current_table.data - baseline.data);
diff_data(isnan(current_table.data) & isnan(baseline.data)) = 0.0;
diff_data(xor(isnan(current_table.data), isnan(baseline.data))) = Inf;
max_abs_diff = max(diff_data, [], 1, "omitnan");

legacy_first = local_row_struct(baseline.columns, baseline.data(1,:));
current_first = local_row_struct(current_table.columns, current_table.data(1,:));

report = struct();
report.case_name = "dd1d_smooth_mms";
report.method = baseline.method;
report.p_order = baseline.p_order;
report.legacy_source = baseline.source;
report.backend = opts.backend;
report.legacy_first = legacy_first;
report.current_first = current_first;
report.columns = baseline.columns;
report.max_abs_diff = max_abs_diff;
report.max_abs_diff_by_column = local_diff_struct(baseline.columns, max_abs_diff);
report.matches_full_table = all(max_abs_diff <= opts.abs_tol | isnan(max_abs_diff));
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

function opts = local_parse_options(varargin)
opts = struct("method", "SIPG", "p_order", 3, "abs_tol", 1.0e-12, ...
    "output_json", "", "backend", "matlab_mfem");
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
        otherwise
            error("Unknown option: %s", name);
    end
end
end
