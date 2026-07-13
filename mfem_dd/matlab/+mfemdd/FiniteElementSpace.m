classdef FiniteElementSpace
    %FINITEELEMENTSPACE Bind a mesh, element collection, and vector dimension.
    properties
        Mesh mfemdd.Mesh
        FEC mfemdd.FiniteElementCollection
        VDim (1,1) double = 1
    end

    methods
        function obj = FiniteElementSpace(mesh, fec, vdim)
            if nargin < 3
                vdim = 1;
            end
            obj.Mesh = mesh;
            obj.FEC = fec;
            obj.VDim = vdim;
        end

        function n = GetTrueVSize(obj)
            n = size(obj.Mesh.Nodes, 1) * obj.VDim;
        end
    end
end
