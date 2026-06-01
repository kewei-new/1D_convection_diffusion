function test_native_pn1d()
%TEST_NATIVE_PN1D Verify the native PN solve reproduces legacy CSV outputs.
test_root = fileparts(mfilename("fullpath"));
addpath(fullfile(fileparts(test_root), "matlab"));
startup_mfem_dd();

baseline = mfemdd.legacy_pn1d_baseline();
native = mfemdd.dd1d_native_pn_baseline("verbose", false);

assert(native.generated_by == "dd1d_native_pn_baseline");
assert(native.output_format == "legacy_dlmwrite_precision");
assert(native.summary.iv_rows == baseline.summary.iv_rows);
assert(native.summary.cv_rows == baseline.summary.cv_rows);
assert(native.summary.transient_rows == baseline.summary.transient_rows);
assert(max(abs(native.iv.right_current - baseline.iv.right_current)) == 0.0);
assert(max(abs(native.cv.Cqs - baseline.cv.Cqs)) == 0.0);
assert(max(abs(native.transient.right_current - baseline.transient.right_current)) == 0.0);

report = mfemdd.compare_legacy_pn1d("backend", "matlab_mfem");
assert(report.current_status == "native_device_mfem");
assert(report.matches_baseline);
assert(report.max_abs_diff == 0.0);
assert(report.max_rel_diff == 0.0);
end
