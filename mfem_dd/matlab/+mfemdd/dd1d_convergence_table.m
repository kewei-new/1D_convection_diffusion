function table_data = dd1d_convergence_table(varargin)
%DD1D_CONVERGENCE_TABLE Produce the DD1D smooth MMS table in the legacy schema.
opts = local_parse_options(varargin{:});
baseline = mfemdd.legacy_dd1d_baseline("method", opts.method, "p_order", opts.p_order);
backend = lower(opts.backend);

if backend == "legacy_matlab"
    table_data = baseline;
    table_data.backend = "legacy_matlab";
    return;
end

if backend == "legacy_runtime"
    table_data = mfemdd.legacy_dd1d_baseline("method", opts.method, ...
        "p_order", opts.p_order, ...
        "run_legacy", true, ...
        "refine_steps", opts.refine_steps, ...
        "ng_base", opts.ng_base, ...
        "verbose", opts.verbose);
    table_data.backend = "legacy_runtime";
    return;
end

row_count = min(opts.refine_steps, size(baseline.data, 1));
data = nan(row_count, size(baseline.data, 2));
data(:, 1) = baseline.data(1:row_count, 1);
for row_id = 1:row_count
    h = baseline.data(row_id, 1);
    elements = round((2.0 * pi) / h);
    metrics = mfemdd.CaseRunner.run("dd1d_smooth_mms", ...
        "elements", elements, ...
        "order", opts.p_order, ...
        "backend", opts.backend);
    data(row_id, 2) = metrics.n_l2_error;
    data(row_id, 4) = NaN;
    data(row_id, 6) = metrics.phi_l2_error;
    data(row_id, 8) = NaN;
    data(row_id, 10) = metrics.E_l2_error;
    data(row_id, 12) = NaN;
end

data(:, 3) = local_orders(data(:, 2));
data(:, 5) = local_orders(data(:, 4));
data(:, 7) = local_orders(data(:, 6));
data(:, 9) = local_orders(data(:, 8));
data(:, 11) = local_orders(data(:, 10));
data(:, 13) = local_orders(data(:, 12));

table_data = struct();
table_data.method = baseline.method;
table_data.p_order = baseline.p_order;
table_data.alpha = baseline.alpha;
table_data.beta = baseline.beta;
table_data.columns = baseline.columns;
table_data.data = data;
table_data.source = "mfem_dd_case_runner";
table_data.generated_by = "dd1d_convergence_table";
table_data.backend = opts.backend;
end

function orders = local_orders(errors)
orders = nan(size(errors));
for k = 2:numel(errors)
    if isfinite(errors(k)) && isfinite(errors(k-1)) && errors(k) > 0 && errors(k-1) > 0
        orders(k) = log2(errors(k-1) / errors(k));
    end
end
end

function opts = local_parse_options(varargin)
opts = struct("method", "SIPG", "p_order", 3, "backend", "matlab_mfem", ...
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
        case "backend"
            opts.backend = string(value);
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
