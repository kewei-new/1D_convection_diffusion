# First-Pass Inventory Review

Generated inventory: 2026-07-13.

| Scope | Result |
| --- | ---: |
| Source files | 1,737 |
| Exact duplicate groups | 329 |
| Files in exact duplicate groups | 1,443 |
| Same-name, different-content candidates | 1,120 |
| Initially unclassified files | 196 |

## Verified First-Pass Coverage

Byte-identical groups can share one numerical verification target. This is a
verification reduction rule only; no source copy is deleted from this result.

- Reference basis and its derivative each have 27 identical copies.
- One-dimensional Gauss quadrature has 19 identical copies.
- Inverse mass, uniform-grid, and Gauss-Lobatto helpers occur in 14--18 exact
  copies.
- The LDG matrix assembly helper has 10 exact copies.

For these groups the future MFEM result must be compared with the selected
canonical legacy implementation once, then with every exposed case using that
implementation. Identical copies require no separate numerical rerun.

## Not Yet Mergeable

Same-name files with different hashes are not containment evidence. In
particular, `main.m`, matrix assembly files, boundary treatments, source and
initial conditions may encode different equation, boundary, flux, limiter, or
time-step settings. They remain pending until the review queue records those
differences.

## Next Classification Order

1. Classify the 196 unclassified files by numerical role.
2. Review non-identical boundary and assembly variants before solver drivers.
3. Map each retained spatial operator to existing MFEM forms/integrators, each
   time method to the common time-stepper, and each limiter to a distinct
   explicitly tested component.
4. Only then select one canonical source per behavior and create a numerical
   regression case.
