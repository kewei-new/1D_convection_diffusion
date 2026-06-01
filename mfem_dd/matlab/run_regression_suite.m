function suite = run_regression_suite(varargin)
%RUN_REGRESSION_SUITE Run the repository's first-stage regression cases.
startup_mfem_dd();
opts = local_parse_options(varargin{:});

if opts.quick
    elements = [8, 16];
else
    elements = [8, 16, 32];
end

repo_root = fileparts(fileparts(mfilename("fullpath")));
if contains(lower(repo_root), "1d_convection_diffusion")
    cases = ["dd1d_smooth_mms", "dd_pn_device"];
else
    cases = ["dd2d_smooth_mms", "dd_pn_device"];
end

rows = struct([]);
idx = 0;
for case_name = cases
    for n = elements
        idx = idx + 1;
        rows(idx).case_name = char(case_name);
        rows(idx).elements = n;
        rows(idx).metrics = mfemdd.CaseRunner.run(case_name, "elements", n, "order", opts.order);
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

function opts = local_parse_options(varargin)
opts = struct("quick", true, "order", 2, "output_dir", "");
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
        otherwise
            error("Unknown option: %s", name);
    end
end
end
