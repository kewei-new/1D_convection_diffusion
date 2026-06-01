function test_native_dd1d()
%TEST_NATIVE_DD1D Verify the self-contained MATLAB backend matches DD1D legacy output.
test_root = fileparts(mfilename("fullpath"));
addpath(fullfile(fileparts(test_root), "matlab"));
startup_mfem_dd();

metrics = run_case("dd1d_smooth_mms", ...
    "elements", 20, ...
    "order", 3, ...
    "backend", "matlab_mfem");
legacy = mfemdd.legacy_dd1d_baseline("method", "SIPG", "p_order", 3);
row = legacy.data(1, :);

assert(metrics.status == "native_matlab_mfem");
assert(abs(metrics.n_l2_error - row(2)) <= 1.0e-12);
assert(abs(metrics.n_linf_error - row(4)) <= 1.0e-12);
assert(abs(metrics.phi_l2_error - row(6)) <= 1.0e-12);
assert(abs(metrics.phi_linf_error - row(8)) <= 1.0e-12);
assert(abs(metrics.E_l2_error - row(10)) <= 1.0e-12);
assert(abs(metrics.E_linf_error - row(12)) <= 1.0e-12);

report = mfemdd.compare_legacy_dd1d("backend", "matlab_mfem", "refine_steps", 2);
assert(report.matches_compared_rows);
assert(report.current_table.generated_by == "dd1d_native_baseline");
end
