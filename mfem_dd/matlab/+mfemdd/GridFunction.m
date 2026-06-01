classdef GridFunction
    %GRIDFUNCTION Nodal field defined on a FiniteElementSpace.
    properties
        FESpace mfemdd.FiniteElementSpace
        Values double
        Name string = "grid_function"
    end

    methods
        function obj = GridFunction(fes, name)
            obj.FESpace = fes;
            obj.Values = zeros(fes.GetTrueVSize(), 1);
            if nargin > 1
                obj.Name = string(name);
            end
        end

        function obj = ProjectCoefficient(obj, coeff)
            obj.Values = coeff.Eval(obj.FESpace.Mesh.Nodes);
        end

        function err = ComputeL2Error(obj, coeff)
            err = sqrt(obj.integrateSquaredDifference(coeff));
        end

        function norm_value = ComputeL2NormCoefficient(obj, coeff)
            norm_value = sqrt(obj.integrateSquaredCoefficient(coeff));
        end
    end

    methods (Access = private)
        function total = integrateSquaredDifference(obj, coeff)
            mesh = obj.FESpace.Mesh;
            total = 0.0;
            if mesh.Dim == 1
                gp = [-1/sqrt(3); 1/sqrt(3)];
                gw = [1; 1];
                for e = 1:mesh.GetNE()
                    ids = mesh.Elements(e,:);
                    x0 = mesh.Nodes(ids(1),1);
                    x1 = mesh.Nodes(ids(2),1);
                    vals = obj.Values(ids);
                    for q = 1:numel(gp)
                        xi = gp(q);
                        x = 0.5 * ((1 - xi) * x0 + (1 + xi) * x1);
                        uh = 0.5 * ((1 - xi) * vals(1) + (1 + xi) * vals(2));
                        exact = coeff.Eval(x);
                        total = total + gw(q) * 0.5 * (x1 - x0) * (uh - exact)^2;
                    end
                end
            else
                gp = [-1/sqrt(3), 1/sqrt(3)];
                gw = [1, 1];
                for e = 1:mesh.GetNE()
                    ids = mesh.Elements(e,:);
                    xy = mesh.Nodes(ids,:);
                    vals = obj.Values(ids);
                    x0 = xy(1,1); x1 = xy(2,1);
                    y0 = xy(1,2); y1 = xy(4,2);
                    jac = 0.25 * (x1 - x0) * (y1 - y0);
                    for i = 1:2
                        for j = 1:2
                            xi = gp(i); eta = gp(j);
                            n = 0.25 * [(1-xi)*(1-eta); (1+xi)*(1-eta); (1+xi)*(1+eta); (1-xi)*(1+eta)];
                            p = n.' * xy;
                            uh = n.' * vals;
                            exact = coeff.Eval(p);
                            total = total + gw(i) * gw(j) * jac * (uh - exact)^2;
                        end
                    end
                end
            end
        end

        function total = integrateSquaredCoefficient(obj, coeff)
            mesh = obj.FESpace.Mesh;
            total = 0.0;
            if mesh.Dim == 1
                gp = [-1/sqrt(3); 1/sqrt(3)];
                gw = [1; 1];
                for e = 1:mesh.GetNE()
                    ids = mesh.Elements(e,:);
                    x0 = mesh.Nodes(ids(1),1);
                    x1 = mesh.Nodes(ids(2),1);
                    for q = 1:numel(gp)
                        xi = gp(q);
                        x = 0.5 * ((1 - xi) * x0 + (1 + xi) * x1);
                        exact = coeff.Eval(x);
                        total = total + gw(q) * 0.5 * (x1 - x0) * exact^2;
                    end
                end
            else
                gp = [-1/sqrt(3), 1/sqrt(3)];
                gw = [1, 1];
                for e = 1:mesh.GetNE()
                    ids = mesh.Elements(e,:);
                    xy = mesh.Nodes(ids,:);
                    x0 = xy(1,1); x1 = xy(2,1);
                    y0 = xy(1,2); y1 = xy(4,2);
                    jac = 0.25 * (x1 - x0) * (y1 - y0);
                    for i = 1:2
                        for j = 1:2
                            xi = gp(i); eta = gp(j);
                            n = 0.25 * [(1-xi)*(1-eta); (1+xi)*(1-eta); (1+xi)*(1+eta); (1-xi)*(1+eta)];
                            p = n.' * xy;
                            exact = coeff.Eval(p);
                            total = total + gw(i) * gw(j) * jac * exact^2;
                        end
                    end
                end
            end
        end
    end
end
