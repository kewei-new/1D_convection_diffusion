function e = evaluate_error_absolute(uh, exact_fun, arg3, arg4, arg5, arg6, arg7)
%EVALUATE_ERROR_ABSOLUTE Compute DG L-infinity error.
%
% New call:
%   e = evaluate_error_absolute(uh, exact_fun, mesh, basis_data, t)
%
% Backward-compatible old call:
%   e = evaluate_error_absolute(uh, exact_fun, ng, mp, h, mid_points, t)

    if nargin == 5 && isstruct(arg3) && isstruct(arg4)
        mesh = arg3;
        basis_data = arg4;
        t = arg5;

        ng = mesh.num_cells;
        mp = basis_data.p_order;
        h = mesh.cell_sizes;
        mid_points = mesh.cell_centers;
        Gp = basis_data.gauss_points;
        basis_values = basis_data.basis_values;
    else
        ng = double(arg3);
        mp = double(arg4);
        h = arg5;
        mid_points = arg6;
        t = arg7;

        if ~isscalar(ng) || ~isscalar(mp)
            error('evaluate_error_absolute: ng and mp must be scalar.');
        end

        ng = round(ng);
        mp = round(mp);

        Gauss_point = generate_Gauss_lobatto(mp);
        Gp = Gauss_point(:,1);
        basis_values = zeros(length(Gp), mp+1);
        for m = 0:mp
            basis_values(:,m+1) = reference_basis(Gp, m);
        end
    end

    if isvector(uh)
        uh = reshape(uh, mp+1, ng);
    else
        [r, c] = size(uh);
        if r == mp+1 && c == ng
            % already correct shape
        elseif r == ng && c == mp+1
            uh = uh.';
        else
            error('evaluate_error_absolute: size of uh is incompatible with mp and ng.');
        end
    end

    e = 0;
    for n = 1:ng
        x = Gp * h(n) + mid_points(n);
        uh_value = basis_values * uh(:,n);
        uex = feval(exact_fun, x, t);
        e = max(e, max(abs(uh_value - uex)));
    end
end
