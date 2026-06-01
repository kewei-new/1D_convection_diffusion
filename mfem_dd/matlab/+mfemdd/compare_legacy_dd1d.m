function report = compare_legacy_dd1d(varargin)
%COMPARE_LEGACY_DD1D Compare current MFEM-DD output with the legacy table.
opts = local_parse_options(varargin{:});
baseline = mfemdd.legacy_dd1d_baseline("method", opts.method, "p_order", opts.p_order);
current = mfemdd.CaseRunner.run("dd1d_smooth_mms", ...
    "elements", round(2*pi / baseline.data(1,1)), ...
    "order", opts.p_order, ...
    "backend", opts.backend);

legacy_first = struct();
legacy_first.h = baseline.data(1, 1);
legacy_first.n_L2 = baseline.data(1, 2);
legacy_first.phi_L2 = baseline.data(1, 6);
legacy_first.E_L2 = baseline.data(1, 10);

report = struct();
report.case_name = "dd1d_smooth_mms";
report.method = baseline.method;
report.p_order = baseline.p_order;
report.legacy_source = baseline.source;
report.legacy_first = legacy_first;
report.current = current;
report.n_L2_abs_diff = abs(current.n_l2_error - legacy_first.n_L2);
report.phi_L2_abs_diff = abs(current.phi_l2_error - legacy_first.phi_L2);
report.E_L2_abs_diff = abs(current.E_l2_error - legacy_first.E_L2);
report.matches_first_row = report.n_L2_abs_diff <= opts.abs_tol && ...
    report.phi_L2_abs_diff <= opts.abs_tol && ...
    report.E_L2_abs_diff <= opts.abs_tol;

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
