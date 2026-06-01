# Migration Report

## Current status

This repository now has an MFEM-first scaffold under `mfem_dd/`.
No legacy folder is superseded or deleted in this commit.

The first baseline adapter for `dd1d_smooth_mms` is now available through:

```matlab
run_case("dd1d_smooth_mms", "order", 3, "backend", "legacy_matlab")
run_case("dd1d_smooth_mms", "elements", 20, "order", 3, "backend", "legacy_runtime")
```

The `legacy_matlab` backend reads the stored legacy table
`DD验阶/DD1D_smooth/V1/result/error_table_SIPG_p3.txt` through the unified
`mfem_dd` entry point. It is a baseline adapter, not a replacement solver.
The `legacy_runtime` backend runs the original MATLAB solver for the requested
mesh size from the same unified case API.
The default `matlab_mfem` backend now resolves to a self-contained MATLAB-native
port inside `mfem_dd` for `dd1d_smooth_mms`; the old runtime remains available
as `backend="legacy_runtime"` and the original projection scaffold remains
available as `backend="mfem_projection"`.

The comparison layer also supports:

```matlab
r = mfemdd.compare_legacy_dd1d("backend", "legacy_runtime", "refine_steps", 5);
```

This route runs the original MATLAB solver from the unified `mfem_dd` API and
compares the recomputed table against the stored legacy table. The default
`matlab_mfem` route instead runs the self-contained `mfemdd.dd1d_native_baseline`
implementation copied into the unified framework.

The regression suite can now emit legacy-runtime metrics through the same CSV
schema:

```matlab
cd mfem_dd/matlab
startup_mfem_dd
run_regression_suite("quick", true, ...
    "order", 3, ...
    "backend", "legacy_runtime", ...
    "include_device", false)
```

The generated `metrics.csv` includes optional legacy columns such as `n_Linf`,
`phi_Linf`, and `E_Linf`; missing non-applicable fields are written as `NaN`.

The C++ `dd1d_mms` app now has an explicit baseline mode:

```bash
dd1d_mms -n 20 -o 3 -b legacy_baseline
```

This mode writes the stored legacy MATLAB DD1D metrics into the same wide
`metrics.csv` schema. It is a C++ baseline adapter, not a C++ port of the
IPDG/LDG/IMEX solver.

The same C++ app now also has a native 1D smooth MMS backend:

```bash
dd1d_mms -n 20 -o 3 -b native_mfem
```

This path assembles and advances the 1D modal DG/IPDG diffusion, LDG Poisson,
and four-stage IMEX update in C++ while still linking against the external MFEM
runtime. It writes the same `metrics.csv` schema with `status=native_cpp_mfem`.

The 1D PN junction physical-output baseline is now exposed through:

```matlab
run_case("dd_pn_device", "backend", "matlab_mfem")
r = mfemdd.compare_legacy_pn1d("backend", "matlab_mfem")
```

It reads the legacy CSV outputs under
`DD*/V3/DD1D_pn_junction/result/pn_junction/`:

- `iv_curve.csv`
- `cv_curve.csv`
- `transient_current.csv`

The current `matlab_mfem` device path recomputes the PN junction solve from a
self-contained source mirror under `mfem_dd/matlab/native/dd1d_pn_junction`.
It then reports the IV/CV/transient tables using the same legacy `dlmwrite`
effective precision so the comparison is against the old files as actually
written, while the raw recomputed tables remain attached to the native baseline
object as `iv_raw`, `cv_raw`, and `transient_raw`.

The C++ `dd_device` app also has a 1D PN baseline mode:

```bash
dd_device -n 16 -o 1 -b legacy_baseline
```

This writes the same selected PN metrics into `metrics.csv`.

A full local validation gate is available:

```powershell
powershell -ExecutionPolicy Bypass -File mfem_dd/tools/run_legacy_validation.ps1
```

It runs the full MATLAB legacy-runtime table comparison, MATLAB PN device CSV
comparison, quick MATLAB legacy-runtime regression, C++ build, CTest, the C++
MMS legacy-baseline app for all 5 embedded legacy mesh rows, the C++ native
MFEM MMS backend against the MATLAB-native full-precision table, and the C++ PN
device legacy-baseline app. Generated reports are written under
`mfem_dd/artifacts/legacy_validation/`.

The validation gate also writes a machine-readable alignment audit:

```powershell
mfem_dd/tools/write_alignment_audit.ps1
```

The audit output is
`mfem_dd/artifacts/alignment_audit/alignment_audit.json`. It classifies each
piece of evidence as native, bridge, or baseline adapter and keeps the remaining
native gaps explicit, so passing validation is not mistaken for permission to
delete legacy folders.

## Tolerances

- Relative L2 difference against legacy baseline: `<= 1e-6`
- Convergence-order absolute difference: `<= 0.15`
- Key physical quantity relative difference: `<= 1e-5`
- Mesh-size absolute difference for printed legacy tables: `<= 1e-6`

## Case table

| Case | Legacy baseline | MATLAB status | C++ status | Delete legacy? | Notes |
| --- | --- | --- | --- | --- | --- |
| `dd1d_smooth_mms` | `DD验阶/DD1D_smooth/V1/result/error_table_SIPG_p3.txt` | default `matlab_mfem` uses self-contained `mfem_dd` MATLAB-native IPDG/LDG/IMEX code and matches all 5 mesh rows within numerical tolerances; explicit `legacy_runtime` remains as old-runtime cross-check; `legacy_matlab` matches the stored table exactly; `mfem_projection` remains scaffold-only | `dd1d_mms -b legacy_baseline` writes stored legacy metrics and validation checks all 5 embedded rows; `dd1d_mms -b native_mfem` runs the native C++ IPDG/LDG/IMEX backend and matches the MATLAB-native table within `2e-12` absolute tolerance | No | First migration target. |
| `dd2d_smooth_mms` | Pending: compare in 2D repository | Scaffold implemented | Scaffold implemented | No | Included for API symmetry. |
| `dd_pn_device` | `DD*/V3/DD1D_pn_junction/result/pn_junction/{iv_curve.csv,cv_curve.csv,transient_current.csv}` | default `matlab_mfem` recomputes the PN solve from `mfem_dd/matlab/native/dd1d_pn_junction` and matches the full legacy CSV outputs; `legacy_matlab` and `legacy_runtime` still read the stored CSV baseline; `mfem_projection` remains scaffold-only | `dd_device -b legacy_baseline` writes the same selected legacy PN metrics; native C++ device solve is not ported | No | MATLAB-native physical path closed; C++ physical path remains baseline-only. |

## Latest 1D Full-Table Comparison

Command:

```matlab
cd mfem_dd/tools
run_legacy_dd1d_compare("backend", "legacy_matlab")
run_legacy_dd1d_compare("backend", "legacy_runtime", "refine_steps", 5)
run_legacy_dd1d_compare("backend", "matlab_mfem", "refine_steps", 5)
cd ../tests
test_legacy_runtime_dd1d
test_native_dd1d
```

For a compact native-backend check:

```matlab
cd mfem_dd/matlab
startup_mfem_dd
r = mfemdd.compare_legacy_dd1d("backend", "matlab_mfem");
disp(r.matches_full_table)
disp(r.max_abs_diff_by_column)
```

Result:

| Backend | Table scope | Max notable differences | Pass |
| --- | --- | --- | --- |
| `legacy_matlab` | 5 mesh rows x 13 legacy columns | all columns `0` | Yes, exact stored-table match |
| `legacy_runtime` | 5 mesh rows x 13 legacy columns | `h=7.3464e-07`, error columns `<=4.1583e-13`, order columns `<=4.5600e-05` | Yes, within numerical tolerances |
| `matlab_mfem` | 5 mesh rows x 13 legacy columns | `h=7.3464e-07`, error columns `<=4.1583e-13`, order columns `<=4.5600e-05` | Yes, self-contained MATLAB-native backend |
| `mfem_projection` | 5 mesh rows x 13 legacy columns | `n_L2=7.8584e-03`, `phi_L2=1.2239e-02`, `E_L2=1.2239e-02`, Linf columns missing (`Inf`) | No |

Latest quick regression CSV checks:

| Command | Rows | First-row status | First-row `n_L2` |
| --- | --- | --- | --- |
| `run_regression_suite("quick", true, "order", 3, "backend", "legacy_runtime", "include_device", false)` | 2 | `legacy_runtime` | `4.507504048373722e-06` |
| `run_regression_suite("quick", true, "order", 3, "backend", "matlab_mfem", "include_device", false)` | 2 | `native_matlab_mfem` | `4.507504048373722e-06` |

Latest C++ baseline check:

| Command scope | Rows | Checked fields | First-row `n_L2` | Last-row `n_L2` | Pass |
| --- | --- | --- | --- | --- | --- |
| `dd1d_mms -n {20,40,80,160,320} -o 3 -b legacy_baseline` | 5 | `n_L2`, `n_Linf`, `phi_L2`, `phi_Linf`, `E_L2`, `E_Linf` | `4.507504e-06` | `6.906262e-10` | Yes, stored baseline adapter |

Latest C++ native check:

| Command scope | Rows | Compared against | Max abs diff | Max rel diff | Pass |
| --- | --- | --- | --- | --- | --- |
| `dd1d_mms -n {20,40,80,160,320} -o 3 -b native_mfem` | 5 | MATLAB-native full-precision table | `1.0264e-12` | `3.0016e-03` | Yes by `2e-12` absolute tolerance; relative difference is dominated by `1e-10`-scale finest-grid errors |

Latest PN/device physical-output check:

| Source | Metric | Legacy output | Native recomputed output | Pass |
| --- | --- | --- | --- | --- |
| `iv_curve.csv` | reverse current at `-1 V` | `-18006` | `-18006` | Yes |
| `iv_curve.csv` | zero-bias current | `0.019603` | `0.019603` | Yes |
| `iv_curve.csv` | forward current at `1 V` | `17981` | `17981` | Yes |
| `iv_curve.csv` | zero-bias `Qmag` | `2612.5` | `2612.5` | Yes |
| `cv_curve.csv` | zero-bias `Cqs` | `-24.835` | `-24.835` | Yes |
| `transient_current.csv` | first finite transient current | `-406000` at `t=0.001875` | `-406000` at `t=0.001875` | Yes |
| `transient_current.csv` | terminal transient current | `-13113` at `t=0.2` | `-13113` at `t=0.2` | Yes |

`mfemdd.compare_legacy_pn1d("backend", "matlab_mfem")` reports zero
difference for these selected metrics and for the full numeric contents of
`iv_curve.csv`, `cv_curve.csv`, and `transient_current.csv` after applying the
same legacy CSV output precision to the recomputed tables.
`dd_device -b legacy_baseline` emits the same selected values in C++ CSV form.

Latest validation-suite check:

| Command | MATLAB runtime status | MATLAB native status | MATLAB first-row `n_L2` | C++ baseline rows | C++ native rows |
| --- | --- | --- | --- | --- | --- |
| `mfem_dd/tools/run_legacy_validation.ps1` | passed | passed | `4.5075040483737219e-06` | 5 passed | 5 passed |

The same validation run also records PN/device baseline status:

| MATLAB PN status | C++ PN status | `iv_forward_current_1v` | `cv_zero_bias_cqs` |
| --- | --- | --- | --- |
| passed | passed | `17981` | `-24.835` |

The old H1/nodal projection scaffold is now only exposed as `mfem_projection`.
It uses the same 1D domain `[0, 2*pi]` and exact functions as the legacy smooth
test, but it is intentionally kept separate from the modal DG solver.

The MATLAB-native smooth MMS gap is now closed for 1D: `matlab_mfem` runs
self-contained `mfem_dd` code for L2 projection, modal DG/IPDG diffusion, LDG
Poisson, and the legacy four-stage IMEX update to `T=1`. The MATLAB-native PN
device gap is also closed for the single-carrier PN surrogate: `matlab_mfem`
recomputes the device experiments and matches the old IV/CV/transient CSV
outputs. The C++ native smooth MMS path is available through
`dd1d_mms -b native_mfem`; it follows the same modal DG/IPDG/LDG/IMEX algorithm
and is validated against the MATLAB-native table with a small absolute
tolerance. The remaining native gap in this repository is the C++ device-oriented
physical solve.

Only the selected 1D PN junction CSV outputs have been validated so far. 2D
device-oriented physical simulations and any independent native MATLAB/C++
device solve remain unvalidated.

## Deletion rule

A legacy folder can only be removed by a later dedicated commit after this table records matching legacy,
MATLAB, and C++ metrics.
