function baseline = legacy_pn1d_baseline()
%LEGACY_PN1D_BASELINE Read the legacy 1D PN junction physical outputs.
repo_root = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
pn_root = local_single_dir(fullfile(repo_root, "DD*", "V3", ...
    "DD1D_pn_junction", "result", "pn_junction"));

iv_file = fullfile(pn_root, "iv_curve.csv");
cv_file = fullfile(pn_root, "cv_curve.csv");
transient_file = fullfile(pn_root, "transient_current.csv");

iv = readtable(iv_file);
cv = readtable(cv_file);
transient = readtable(transient_file);

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

baseline = struct();
baseline.case_name = "dd_pn_device";
baseline.device = "pn_junction";
baseline.generated_by = "stored_legacy_device_csv";
baseline.source_root = string(pn_root);
baseline.iv_source = string(iv_file);
baseline.cv_source = string(cv_file);
baseline.transient_source = string(transient_file);
baseline.summary = summary;
baseline.iv = iv;
baseline.cv = cv;
baseline.transient = transient;
end

function folder = local_single_dir(pattern)
matches = dir(pattern);
matches = matches([matches.isdir]);
if isempty(matches)
    error("Legacy PN1D baseline folder not found: %s", pattern);
end
folder = fullfile(matches(1).folder, matches(1).name);
end

function value = local_lookup(tbl, key_name, key_value, value_name)
keys = tbl.(key_name);
idx = find(abs(keys - key_value) <= 1.0e-12, 1, "first");
if isempty(idx)
    error("Could not find %s=%g in legacy table.", key_name, key_value);
end
value = tbl.(value_name)(idx);
end
