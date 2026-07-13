function solver = build_linear_solver(A)
%BUILD_LINEAR_SOLVER Cache a sparse linear solve when possible.

    solver.A = A;
    solver.use_decomposition = false;

    if exist('decomposition', 'builtin') == 5 || exist('decomposition', 'file') == 2
        solver.factor = decomposition(A, 'lu');
        solver.use_decomposition = true;
    end
end
