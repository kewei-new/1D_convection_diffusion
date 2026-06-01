classdef LinearForm
    %LINEARFORM Minimal right-hand side object.
    properties
        FESpace mfemdd.FiniteElementSpace
        Vector double
    end

    methods
        function obj = LinearForm(fes)
            obj.FESpace = fes;
            obj.Vector = zeros(fes.GetTrueVSize(), 1);
        end

        function obj = ProjectCoefficient(obj, coeff)
            obj.Vector = coeff.Eval(obj.FESpace.Mesh.Nodes);
        end
    end
end
