# MFEM-First Drift-Diffusion Workspace

This directory is the new MFEM-first implementation surface for the drift-diffusion work in this repository.
The legacy folders remain in place until a case has a passing migration report.

## Layout

- `matlab/` contains the MATLAB reference implementation with MFEM-style objects.
- `cpp/` contains C++ MFEM applications linked against an external MFEM build.
- `cases/` contains canonical case definitions and notes.
- `tests/` contains smoke and regression entry points.
- `docs/` contains migration and validation reports.
- `tools/` contains helper scripts.
- `artifacts/` is the default generated-output location and is ignored except for `.gitkeep`.

## MATLAB quick start

```matlab
cd mfem_dd/matlab
startup_mfem_dd
metrics = run_case("dd1d_smooth_mms", "elements", 16, "order", 2)
suite = run_regression_suite("quick", true)
```

The MATLAB layer intentionally mirrors MFEM concepts: `Mesh`, `FiniteElementCollection`,
`FiniteElementSpace`, `GridFunction`, `Coefficient`, `BilinearForm`, `LinearForm`,
`Solver`, and `TimeStepper`.

## C++ quick start

Configure with a previously built MFEM tree. For example:

```powershell
cmake -S mfem_dd/cpp -B build/mfem_dd -DMFEM_DIR=E:/Projects/codes/2D_convection_diffusion/DD_model/模拟/V6/build
cmake --build build/mfem_dd
ctest --test-dir build/mfem_dd --output-on-failure
```

If `MFEM_DIR` is not set, the CMake project also checks `MFEM_ROOT`, `../mfem-4.9/build`,
and `../../mfem-4.9/build`.

## Migration rule

Do not delete a legacy case until:

1. The legacy baseline is recorded.
2. MATLAB metrics pass the tolerance in `docs/migration_report.md`.
3. C++ metrics pass the same tolerance.
4. The report marks the legacy folder as superseded.
