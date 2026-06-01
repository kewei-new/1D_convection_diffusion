PN-junction extension on top of the DD1D_smooth code
====================================================

Main additions
--------------
1. config/build_pn_junction_params.m
   Device-oriented parameter builder.
2. model/nd_profile.m
   Smooth doping profile Nd(x).
3. model/build_problem_definition.m
   Now supports two modes:
      - smooth_test            : original manufactured-solution mode
      - pn_junction_device     : device mode
4. solver/assemble_transport_rhs_single.m
   Added Dirichlet contact treatment for carrier transport.
5. solver/solve_1D_DD.m
   Added restart, history recording, and terminal current output.
6. postprocess/compute_current_1D.m
   Computes J = -D n_x + mu n E.
7. postprocess/compute_space_charge_proxy_1D.m
   Provides quasi-static charge / depletion-width proxies.
8. main/main_pn_junction_experiments.m
   Runs four experiments:
      (a) equilibrium
      (b) I-V sweep
      (c) quasi-static C-V via dQ/dV
      (d) forward-to-reverse switching transient

How to run
----------
In MATLAB, from the project root:

    cd main
    main_pn_junction_experiments

Notes
-----
- This extension stays as close as possible to the original single-carrier
  IMEX + IPDG/LDG structure and interface.
- Therefore it is a device-oriented single-carrier p-n junction surrogate,
  not a full bipolar DD implementation with both electrons and holes.
- If you want a strict bipolar DD model, the next step is to duplicate the
  carrier equation into (n,p), add recombination R(n,p), and change the
  Poisson right-hand side from (n-Nd) to (p-n+Nd).
