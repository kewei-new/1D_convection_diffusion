function test_legacy_pn1d()
%TEST_LEGACY_PN1D Verify unified PN device metrics match legacy CSV output.
test_root = fileparts(mfilename("fullpath"));
addpath(fullfile(fileparts(test_root), "matlab"));
startup_mfem_dd();

baseline = mfemdd.legacy_pn1d_baseline();
assert(baseline.summary.iv_rows == 9);
assert(baseline.summary.cv_rows == 7);
assert(baseline.summary.transient_rows == 108);
assert(abs(baseline.summary.iv_reverse_current_minus1v + 18006.0) <= 1.0e-12);
assert(abs(baseline.summary.iv_zero_bias_current - 0.019603) <= 1.0e-12);
assert(abs(baseline.summary.iv_forward_current_1v - 17981.0) <= 1.0e-12);
assert(abs(baseline.summary.iv_zero_bias_qmag - 2612.5) <= 1.0e-12);
assert(abs(baseline.summary.cv_zero_bias_cqs + 24.835) <= 1.0e-12);
assert(abs(baseline.summary.transient_first_finite_current + 406000.0) <= 1.0e-8);

metrics = run_case("dd_pn_device", "backend", "matlab_mfem");
assert(metrics.status == "legacy_device_baseline");
assert(abs(metrics.iv_forward_current_1v - baseline.summary.iv_forward_current_1v) <= 1.0e-12);

report = mfemdd.compare_legacy_pn1d("backend", "matlab_mfem");
assert(report.matches_baseline);
end
