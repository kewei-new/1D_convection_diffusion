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

The comparison layer also supports:

```matlab
r = mfemdd.compare_legacy_dd1d("backend", "legacy_runtime", "refine_steps", 5);
```

This route runs the original MATLAB solver from the unified `mfem_dd` API and
compares the recomputed table against the stored legacy table. It is used as a
runtime adapter while the native MFEM-style solver is still being aligned.

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

## Tolerances

- Relative L2 difference against legacy baseline: `<= 1e-6`
- Convergence-order absolute difference: `<= 0.15`
- Key physical quantity relative difference: `<= 1e-5`
- Mesh-size absolute difference for printed legacy tables: `<= 1e-6`

## Case table

| Case | Legacy baseline | MATLAB status | C++ status | Delete legacy? | Notes |
| --- | --- | --- | --- | --- | --- |
| `dd1d_smooth_mms` | `DD验阶/DD1D_smooth/V1/result/error_table_SIPG_p3.txt` | `legacy_matlab` backend matches the full stored baseline table exactly; `legacy_runtime` recomputes all 5 mesh rows within numerical tolerances and is exposed through `run_case`; native `matlab_mfem` still differs | `dd1d_mms -b legacy_baseline` writes stored legacy metrics; native C++ projection still differs | No | First migration target. |
| `dd2d_smooth_mms` | Pending: compare in 2D repository | Scaffold implemented | Scaffold implemented | No | Included for API symmetry. |
| `dd_pn_device` | Pending: PN/device folders | Scaffold only | Scaffold only | No | Records doping/charge proxy only. |

## Latest 1D Full-Table Comparison

Command:

```matlab
cd mfem_dd/tools
run_legacy_dd1d_compare("backend", "legacy_matlab")
run_legacy_dd1d_compare("backend", "legacy_runtime", "refine_steps", 5)
cd ../tests
test_legacy_runtime_dd1d
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
| `matlab_mfem` | 5 mesh rows x 13 legacy columns | `n_L2=7.8584e-03`, `phi_L2=1.2239e-02`, `E_L2=1.2239e-02`, Linf columns missing (`Inf`) | No |

Latest quick legacy-runtime regression CSV check:

| Command | Rows | First-row status | First-row `n_L2` |
| --- | --- | --- | --- |
| `run_regression_suite("quick", true, "order", 3, "backend", "legacy_runtime", "include_device", false)` | 2 | `legacy_runtime` | `4.507504048373722e-06` |

Latest C++ baseline check:

| Command | First-row `n_L2` | First-row `phi_L2` | First-row `E_L2` | Pass |
| --- | --- | --- | --- | --- |
| `dd1d_mms -n 20 -o 3 -b legacy_baseline` | `4.507504e-06` | `5.547920e-06` | `5.658847e-06` | Yes, stored baseline adapter |

The native MATLAB backend now uses the same 1D domain `[0, 2*pi]` and exact
functions as the legacy smooth test, but it is still an H1/nodal scaffold rather
than the old modal DG + IPDG/LDG/IMEX solver. It must not be used to delete
legacy folders.

The first confirmed native gap is in
`mfemdd.CaseRunner.runSmoothMMS1D`: it projects the exact manufactured
solution into an H1 space and measures projection error. The legacy solver
instead performs L2 projection of the initial carrier density, assembles modal
DG/IPDG diffusion and LDG Poisson matrices, then advances the carrier equation
with the legacy four-stage IMEX method to `T=1`. Native replacement work should
start by porting those operators behind MFEM-style MATLAB classes rather than
adjusting tolerances around the current scaffold.

No PN/device-oriented physical simulation has been validated against legacy
outputs yet. Those cases remain scaffold-only until their old entry points,
baseline metrics, and comparable output files are identified.

## Deletion rule

A legacy folder can only be removed by a later dedicated commit after this table records matching legacy,
MATLAB, and C++ metrics.
