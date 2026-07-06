function test_legacy_runtime_dd1d()
%TEST_LEGACY_RUNTIME_DD1D Verify run_case reads a controlled legacy DD1D row.
test_root = fileparts(mfilename("fullpath"));
addpath(fullfile(fileparts(test_root), "matlab"));
startup_mfem_dd();

metrics = run_case("dd1d_smooth_mms", ...
    "elements", 40, ...
    "order", 3, ...
    "backend", "legacy_matlab");
legacy = mfemdd.legacy_dd1d_baseline("method", "SIPG", "p_order", 3);
row = legacy.data(2, :);

assert(metrics.status == "legacy_matlab_baseline");
assert(abs(metrics.n_l2_error - row(2)) <= 1.0e-12);
assert(abs(metrics.n_linf_error - row(4)) <= 1.0e-12);
assert(abs(metrics.phi_l2_error - row(6)) <= 1.0e-12);
assert(abs(metrics.phi_linf_error - row(8)) <= 1.0e-12);
assert(abs(metrics.E_l2_error - row(10)) <= 1.0e-12);
assert(abs(metrics.E_linf_error - row(12)) <= 1.0e-12);
end
