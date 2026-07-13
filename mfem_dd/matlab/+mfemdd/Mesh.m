classdef Mesh
    %MESH Minimal MFEM-style structured mesh for MATLAB reference runs.
    properties
        Dim (1,1) double
        Nodes double
        Elements double
        ElementType string
        Extents double
    end

    methods
        function obj = Mesh(dim, nodes, elements, element_type, extents)
            obj.Dim = dim;
            obj.Nodes = nodes;
            obj.Elements = elements;
            obj.ElementType = string(element_type);
            obj.Extents = extents;
        end

        function h = HMax(obj)
            if obj.Dim == 1
                h = max(diff(obj.Nodes(:,1)));
            else
                xs = unique(obj.Nodes(:,1));
                ys = unique(obj.Nodes(:,2));
                h = max([max(diff(xs)), max(diff(ys))]);
            end
        end

        function n = GetNE(obj)
            n = size(obj.Elements, 1);
        end
    end

    methods (Static)
        function mesh = MakeCartesian1D(elements, length_x)
            if nargin < 2
                length_x = 1.0;
            end
            nodes = linspace(0.0, length_x, elements + 1).';
            conn = [(1:elements).', (2:elements+1).'];
            mesh = mfemdd.Mesh(1, nodes, conn, "segment", [length_x, 0.0]);
        end

        function mesh = MakeCartesian2D(nx, ny, length_x, length_y)
            if nargin < 3
                length_x = 1.0;
            end
            if nargin < 4
                length_y = 1.0;
            end
            [X, Y] = meshgrid(linspace(0.0, length_x, nx + 1), linspace(0.0, length_y, ny + 1));
            nodes = [X(:), Y(:)];
            id = @(i,j) j * (nx + 1) + i + 1;
            conn = zeros(nx * ny, 4);
            row = 0;
            for j = 0:ny-1
                for i = 0:nx-1
                    row = row + 1;
                    conn(row,:) = [id(i,j), id(i+1,j), id(i+1,j+1), id(i,j+1)];
                end
            end
            mesh = mfemdd.Mesh(2, nodes, conn, "quadrilateral", [length_x, length_y]);
        end
    end
end
