# MFEM Implementation Matrix

This matrix is the implementation rule for the classified legacy code. Prefer
the listed MFEM primitive first. Add a local routine only when the historical
formula cannot be represented by these primitives without changing the method.

| Legacy responsibility | Use directly from MFEM 4.9 | Local work only when required |
| --- | --- | --- |
| Uniform mesh and polynomial space | `Mesh::MakeCartesian1D`, `L2_FECollection`, `FiniteElementSpace`, `GridFunction` | nonuniform legacy mesh mapping |
| Coefficients and initial/boundary data | `FunctionCoefficient`, `VectorFunctionCoefficient`, `GridFunction::ProjectCoefficient` | only legacy modal-data adaptation |
| Mass and diffusion volume terms | `BilinearForm`, `MassIntegrator`, `DiffusionIntegrator` | none for the standard terms |
| IPDG/SIPG face terms | `DGDiffusionIntegrator` plus MFEM boundary integrators | a matched face term only when the legacy penalty/flux differs |
| LDG derivative coupling | `MixedBilinearForm`, `MixedScalarDerivativeIntegrator` | LDG alternating face flux and its boundary convention |
| Carrier convection and DG flux | `ConvectionIntegrator`, `DGTraceIntegrator` | Lax--Friedrichs or limiter-dependent flux when not equivalent |
| Explicit RK | `ForwardEulerSolver`, `RK2Solver`, `RK3SSPSolver`, `RK4Solver` | none for matching explicit tableaux |
| IMEX-RK | existing assembled mass/diffusion/transport operators | table-driven stage update; MFEM 4.9 has no matching generic IMEX-RK solver |
| Algebraic solve | `CGSolver`, `GMRESSolver`, `GSSmoother`, `BlockILU`, `UMFPackSolver` | only preconditioner selection and tolerance policy |
| Quadrature | `IntegrationRules`, `IntegrationRule` | legacy rule reproduction only when an exact rule differs |
| Limiter | MFEM vectors/spaces and element/face access | limiter formula; no generic DG limiter should be assumed available |
| Output and postprocessing | `GridFunction`, MFEM mesh/field evaluation | legacy-error norms or plots not expressible by the standard output path |

The matrix does not assert numerical equivalence. Each retained equation and
boundary combination must first match its selected legacy output, then may be
used as the canonical implementation for exact source duplicates.
