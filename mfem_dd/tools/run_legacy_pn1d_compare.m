function report = run_legacy_pn1d_compare(varargin)
%RUN_LEGACY_PN1D_COMPARE Entry point for PN device baseline comparison.
root = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(root, "matlab"));
out_file = fullfile(root, "artifacts", "legacy_pn1d_compare", "compare.json");
report = mfemdd.compare_legacy_pn1d("output_json", out_file, varargin{:});
disp(jsonencode(report, PrettyPrint=true));
end
