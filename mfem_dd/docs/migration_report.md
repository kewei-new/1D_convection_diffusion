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

## Tolerances

- Relative L2 difference against legacy baseline: `<= 1e-6`
- Convergence-order absolute difference: `<= 0.15`
- Key physical quantity relative difference: `<= 1e-5`

## Case table

| Case | Legacy baseline | MATLAB status | C++ status | Delete legacy? | Notes |
| --- | --- | --- | --- | --- | --- |
| `dd1d_smooth_mms` | `DD验阶/DD1D_smooth/V1/result/error_table_SIPG_p3.txt` | `legacy_matlab` backend matches the first stored baseline row exactly; native `matlab_mfem` still differs | Scaffold implemented, not compared to legacy yet | No | First migration target. |
| `dd2d_smooth_mms` | Pending: compare in 2D repository | Scaffold implemented | Scaffold implemented | No | Included for API symmetry. |
| `dd_pn_device` | Pending: PN/device folders | Scaffold only | Scaffold only | No | Records doping/charge proxy only. |

## Latest 1D comparison

Command:

```matlab
cd mfem_dd/tools
run_legacy_dd1d_compare("backend", "legacy_matlab")
run_legacy_dd1d_compare("backend", "matlab_mfem")
```

Result:

| Backend | n_L2 diff | phi_L2 diff | E_L2 diff | Pass |
| --- | ---: | ---: | ---: | --- |
| `legacy_matlab` | `0` | `0` | `0` | Yes, as a baseline adapter |
| `matlab_mfem` | `7.85395938243e-03` | `1.22332891129e-02` | `1.22331781859e-02` | No |

The native MATLAB backend now uses the same 1D domain `[0, 2*pi]` and exact
functions as the legacy smooth test, but it is still an H1/nodal scaffold rather
than the old modal DG + IPDG/LDG/IMEX solver. It must not be used to delete
legacy folders.

## Deletion rule

A legacy folder can only be removed by a later dedicated commit after this table records matching legacy,
MATLAB, and C++ metrics.
