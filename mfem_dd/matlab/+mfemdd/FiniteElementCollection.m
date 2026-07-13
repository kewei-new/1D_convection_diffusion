classdef FiniteElementCollection
    %FINITEELEMENTCOLLECTION MFEM-style finite element metadata.
    properties
        Family string
        Order (1,1) double
        Dim (1,1) double
    end

    methods
        function obj = FiniteElementCollection(family, order, dim)
            if nargin < 1
                family = "H1";
            end
            if nargin < 2
                order = 1;
            end
            if nargin < 3
                dim = 1;
            end
            obj.Family = string(family);
            obj.Order = order;
            obj.Dim = dim;
        end
    end
end
