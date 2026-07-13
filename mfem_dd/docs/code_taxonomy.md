# Legacy Code Taxonomy

The inventory is organized by numerical responsibility rather than historical
folders. The working categories are: space discretization (IPDG/SIPG, LDG and
DG spaces), time discretization (IMEX-RK, TVD/SSP-RK and other RK methods),
equation type (elliptic, parabolic/transport and drift-diffusion), boundary
condition (periodic, Dirichlet, Neumann, Robin and mixed), limiter, and shared
infrastructure (quadrature, mesh, algebraic solves, postprocessing and APIs).

Run:

```powershell
powershell -ExecutionPolicy Bypass -File mfem_dd\tools\inventory_legacy_code.ps1
```

The generated `source_inventory.csv` is a provisional multi-tag file list.
`exact_duplicate_groups.csv` identifies byte-identical copies; one canonical
implementation can be validated before its copies. `overlap_candidates.csv`
lists same-name, different-content candidates requiring a manual numerical
containment decision. No overlap candidate is treated as deletable by this
tool.
