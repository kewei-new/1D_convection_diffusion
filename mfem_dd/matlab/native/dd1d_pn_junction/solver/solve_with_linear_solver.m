function x = solve_with_linear_solver(solver, rhs)
%SOLVE_WITH_LINEAR_SOLVER Solve a linear system using cached factorization.

    if solver.use_decomposition
        x = solver.factor \ rhs;
    else
        x = solver.A \ rhs;
    end
end
