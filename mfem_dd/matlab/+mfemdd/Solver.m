classdef Solver
    %SOLVER Linear solver facade matching the MFEM workflow.
    methods (Static)
        function x = Solve(A, b)
            x = A \ b;
        end
    end
end
