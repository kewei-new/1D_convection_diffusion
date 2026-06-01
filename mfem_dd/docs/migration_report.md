# Migration Report

## Current status

This repository now has an MFEM-first scaffold under `mfem_dd/`.
No legacy folder is superseded or deleted in this commit.

The first baseline adapter for `dd1d_smooth_mms` is now available through:

```matlab
run_case("dd1d_smooth_mms", "order", 3, "backend", "legacy_matlab")
```

This backend reads the stored legacy table
`DD验阶/DD1D_smooth/V1/result/error_table_SIPG_p3.txt` through the unified
`mfem_dd` entry point. It is a baseline adapter, not a replacement solver.

The comparison layer also supports:

```matlab
r = mfemdd.compare_legacy_dd1d("backend", "legacy_runtime", "refine_steps", 2);
```

This route runs the original MATLAB solver from the unified `mfem_dd` API and
compares the recomputed table against the stored legacy table. It is used as a
runtime adapter while the native MFEM-style solver is still being aligned.

## Tolerances

- Relative L2 difference against legacy baseline: `<= 1e-6`
- Convergence-order absolute difference: `<= 0.15`
- Key physical quantity relative difference: `<= 1e-5`
- Mesh-size absolute difference for printed legacy tables: `<= 1e-6`

## Case table

| Case | Legacy baseline | MATLAB status | C++ status | Delete legacy? | Notes |
| --- | --- | --- | --- | --- | --- |
| `dd1d_smooth_mms` | `DD验阶/DD1D_smooth/V1/result/error_table_SIPG_p3.txt` | `legacy_matlab` backend matches the full stored baseline table exactly; `legacy_runtime` recomputes the first 2 mesh rows within numerical tolerances; native `matlab_mfem` still differs | Scaffold implemented, not compared to legacy yet | No | First migration target. |
| `dd2d_smooth_mms` | Pending: compare in 2D repository | Scaffold implemented | Scaffold implemented | No | Included for API symmetry. |
| `dd_pn_device` | Pending: PN/device folders | Scaffold only | Scaffold only | No | Records doping/charge proxy only. |

## Latest 1D Full-Table Comparison

Command:

```matlab
cd mfem_dd/tools
run_legacy_dd1d_compare("backend", "legacy_matlab")
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
| `legacy_runtime` | first 2 mesh rows x 13 legacy columns | `h=7.3464e-07`, error columns `<=4.1583e-13`, order columns `<=4.2415e-05` | Yes, within numerical tolerances; not a full-table run yet |
| `matlab_mfem` | 5 mesh rows x 13 legacy columns | `n_L2=7.8584e-03`, `phi_L2=1.2239e-02`, `E_L2=1.2239e-02`, Linf columns missing (`Inf`) | No |

The native MATLAB backend now uses the same 1D domain `[0, 2*pi]` and exact
functions as the legacy smooth test, but it is still an H1/nodal scaffold rather
than the old modal DG + IPDG/LDG/IMEX solver. It must not be used to delete
legacy folders.

No PN/device-oriented physical simulation has been validated against legacy
outputs yet. Those cases remain scaffold-only until their old entry points,
baseline metrics, and comparable output files are identified.

## Deletion rule

A legacy folder can only be removed by a later dedicated commit after this table records matching legacy,
MATLAB, and C++ metrics.
