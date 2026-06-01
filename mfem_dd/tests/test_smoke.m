function test_smoke()
%TEST_SMOKE Fast MATLAB smoke checks for the MFEM-style implementation.
test_root = fileparts(mfilename("fullpath"));
addpath(fullfile(fileparts(test_root), "matlab"));
startup_mfem_dd();
m1 = run_case("dd1d_smooth_mms", "elements", 8, "order", 2, "backend", "mfem_projection");
assert(isfinite(m1.n_l2_error) && m1.n_l2_error >= 0.0);
assert(m1.dimension == 1);

m1_legacy = run_case("dd1d_smooth_mms", "order", 3, "backend", "legacy_matlab");
assert(m1_legacy.n_l2_error > 0.0);
assert(m1_legacy.status == "legacy_matlab_baseline");

legacy_report = mfemdd.compare_legacy_dd1d("backend", "legacy_matlab");
assert(legacy_report.matches_full_table);

m2 = run_case("dd2d_smooth_mms", "elements", 4, "order", 1, "backend", "mfem_projection");
assert(isfinite(m2.phi_l2_error) && m2.phi_l2_error >= 0.0);
assert(m2.dimension == 2);

suite = run_regression_suite("quick", true, "order", 1, "backend", "mfem_projection");
assert(numel(suite.rows) >= 2);
end
