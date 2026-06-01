function report = compare_legacy_pn1d(varargin)
%COMPARE_LEGACY_PN1D Compare unified device metrics with legacy PN CSVs.
opts = local_parse_options(varargin{:});
baseline = mfemdd.legacy_pn1d_baseline();
current = mfemdd.CaseRunner.run("dd_pn_device", ...
    "elements", opts.elements, ...
    "order", opts.order, ...
    "backend", opts.backend);

fields = ["iv_rows", "cv_rows", "transient_rows", ...
    "iv_reverse_current_minus1v", "iv_zero_bias_current", ...
    "iv_forward_current_1v", "iv_zero_bias_qmag", ...
    "cv_zero_bias_cqs", "transient_first_finite_time", ...
    "transient_first_finite_current", "transient_terminal_time", ...
    "transient_terminal_current"];

max_abs_diff = 0.0;
max_rel_diff = 0.0;
field_report = struct();
matches = true;
for k = 1:numel(fields)
    field = fields(k);
    legacy_value = baseline.summary.(field);
    if isfield(current, field)
        current_value = current.(field);
    else
        current_value = NaN;
    end
    abs_diff = abs(current_value - legacy_value);
    rel_diff = abs_diff / max(abs(legacy_value), realmin);
    pass = abs_diff <= max(opts.abs_tol, opts.quantity_relative .* max(abs(legacy_value), realmin));
    field_report.(field) = struct("legacy", legacy_value, ...
        "current", current_value, ...
        "abs_diff", abs_diff, ...
        "rel_diff", rel_diff, ...
        "pass", pass);
    matches = matches && pass;
    max_abs_diff = max(max_abs_diff, abs_diff);
    max_rel_diff = max(max_rel_diff, rel_diff);
end

report = struct();
report.case_name = "dd_pn_device";
report.device = "pn_junction";
report.backend = opts.backend;
report.legacy_source_root = baseline.source_root;
report.current_status = current.status;
report.fields = fields;
report.field_report = field_report;
report.max_abs_diff = max_abs_diff;
report.max_rel_diff = max_rel_diff;
report.matches_baseline = matches;
report.tolerances = struct("abs_tol", opts.abs_tol, ...
    "quantity_relative", opts.quantity_relative);

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

function opts = local_parse_options(varargin)
opts = struct("backend", "matlab_mfem", "elements", 16, "order", 1, ...
    "abs_tol", 1.0e-12, "quantity_relative", 1.0e-5, ...
    "output_json", "");
if mod(numel(varargin), 2) ~= 0
    error("Options must be name/value pairs.");
end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k+1};
    switch name
        case "backend"
            opts.backend = string(value);
        case "elements"
            opts.elements = value;
        case "order"
            opts.order = value;
        case "abs_tol"
            opts.abs_tol = value;
        case "quantity_relative"
            opts.quantity_relative = value;
        case "output_json"
            opts.output_json = string(value);
        otherwise
            error("Unknown option: %s", name);
    end
end
end
