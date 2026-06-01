function report = run_legacy_dd1d_compare(varargin)
%RUN_LEGACY_DD1D_COMPARE Entry point for baseline-vs-current DD1D comparison.
root = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(root, "matlab"));
out_file = fullfile(root, "artifacts", "legacy_dd1d_compare", "compare.json");
report = mfemdd.compare_legacy_dd1d("output_json", out_file, varargin{:});
disp(jsonencode(report, PrettyPrint=true));
end
