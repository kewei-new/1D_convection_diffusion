classdef Coefficient
    %COEFFICIENT Wrap a scalar function handle with MFEM-like Eval semantics.
    properties
        Handle function_handle
        Name string = "coefficient"
    end

    methods
        function obj = Coefficient(handle, name)
            obj.Handle = handle;
            if nargin > 1
                obj.Name = string(name);
            end
        end

        function values = Eval(obj, points)
            values = obj.Handle(points);
            values = values(:);
        end
    end
end
