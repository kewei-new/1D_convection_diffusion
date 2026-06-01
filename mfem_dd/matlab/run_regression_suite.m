function suite = run_regression_suite(varargin)
%RUN_REGRESSION_SUITE Run the repository's first-stage regression cases.
startup_mfem_dd();
opts = local_parse_options(varargin{:});

repo_root = fileparts(fileparts(mfilename("fullpath")));
is_1d_repo = contains(lower(repo_root), "1d_convection_diffusion");
elements = local_select_elements(is_1d_repo, opts);
include_device = opts.include_device;
if isempty(include_device)
    include_device = lower(opts.backend) == "matlab_mfem";
end

if is_1d_repo
    cases = "dd1d_smooth_mms";
else
    cases = "dd2d_smooth_mms";
end
if include_device
    cases = [cases, "dd_pn_device"];
end

rows = struct([]);
idx = 0;
for case_name = cases
    for n = elements
        idx = idx + 1;
        rows(idx).case_name = char(case_name);
        rows(idx).elements = n;
        rows(idx).metrics = mfemdd.CaseRunner.run(case_name, ...
            "elements", n, ...
            "order", opts.order, ...
            "backend", opts.backend);
    end
end

suite = struct();
suite.generated_at = char(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss"));
suite.rows = rows;
suite.tolerances = mfemdd.CaseRunner.defaultTolerances();

out_dir = opts.output_dir;
if strlength(out_dir) == 0
    out_dir = fullfile(repo_root, "artifacts", "matlab_regression");
end
if ~exist(out_dir, "dir")
    mkdir(out_dir);
end
mfemdd.write_metrics_csv(rows, fullfile(out_dir, "metrics.csv"));
end

function elements = local_select_elements(is_1d_repo, opts)
if ~isempty(opts.elements)
    elements = opts.elements;
    return;
end
if lower(opts.backend) == "legacy_runtime" || lower(opts.backend) == "legacy_matlab"
    if is_1d_repo
        if opts.quick
            elements = [20, 40];
        else
            elements = [20, 40, 80, 160, 320];
        end
    else
        if opts.quick
            elements = [6, 12];
        else
            elements = [6, 12, 24, 48];
        end
    end
else
    if opts.quick
        elements = [8, 16];
    else
        elements = [8, 16, 32];
    end
end
end

function opts = local_parse_options(varargin)
opts = struct("quick", true, "order", 2, "output_dir", "", ...
    "backend", "matlab_mfem", "elements", [], "include_device", []);
if mod(numel(varargin), 2) ~= 0
    error("Options must be name/value pairs.");
end
for k = 1:2:numel(varargin)
    name = lower(string(varargin{k}));
    value = varargin{k+1};
    switch name
        case "quick"
            opts.quick = logical(value);
        case "order"
            opts.order = value;
        case "output_dir"
            opts.output_dir = string(value);
        case "backend"
            opts.backend = string(value);
        case "elements"
            opts.elements = value;
        case "include_device"
            opts.include_device = logical(value);
        otherwise
            error("Unknown option: %s", name);
    end
end
end
