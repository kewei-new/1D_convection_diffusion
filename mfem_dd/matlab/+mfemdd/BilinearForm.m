classdef BilinearForm
    %BILINEARFORM Minimal assembled operator placeholder.
    properties
        FESpace mfemdd.FiniteElementSpace
        Integrators cell = {}
        Matrix double = []
    end

    methods
        function obj = BilinearForm(fes)
            obj.FESpace = fes;
        end

        function obj = AddDomainIntegrator(obj, name, coefficient)
            if nargin < 3
                coefficient = [];
            end
            obj.Integrators{end+1} = struct("name", string(name), "coefficient", coefficient);
        end

        function obj = Assemble(obj)
            n = obj.FESpace.GetTrueVSize();
            obj.Matrix = speye(n);
        end
    end
end
