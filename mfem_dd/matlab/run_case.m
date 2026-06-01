function metrics = run_case(case_name, varargin)
%RUN_CASE Run one canonical MFEM-DD case and return metrics.
startup_mfem_dd();
metrics = mfemdd.CaseRunner.run(case_name, varargin{:});
end
