function e = evaluate_error_Lnorm(uh, exact_fun, L_order, arg4, arg5, arg6, arg7, arg8)
%EVALUATE_ERROR_LNORM Compute DG error in L^p norm.
%
% New call:
%   e = evaluate_error_Lnorm(uh, exact_fun, L_order, mesh, basis_data, t)
%
% Backward-compatible old call:
%   e = evaluate_error_Lnorm(uh, exact_fun, L_order, ng, mp, h, mid_points, t)

    if nargin == 6 && isstruct(arg4) && isstruct(arg5)
        mesh = arg4;
        basis_data = arg5;
        t = arg6;

        ng = mesh.num_cells;
        mp = basis_data.p_order;
        h = mesh.cell_sizes;
        mid_points = mesh.cell_centers;
        Gp = basis_data.gauss_points;
        Gw = basis_data.gauss_weights;
        basis_values = basis_data.basis_values;
    else
        ng = double(arg4);
        mp = double(arg5);
        h = arg6;
        mid_points = arg7;
        t = arg8;

        if ~isscalar(ng) || ~isscalar(mp)
            error('evaluate_error_Lnorm: ng and mp must be scalar.');
        end

        ng = round(ng);
        mp = round(mp);

        Gauss_point = generate_Gauss_lobatto(mp);
        Gp = Gauss_point(:,1);
        Gw = Gauss_point(:,2);

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
            error('evaluate_error_Lnorm: size of uh is incompatible with mp and ng.');
        end
    end

    e = 0;

    for n = 1:ng
        x = Gp * h(n) + mid_points(n);
        uh_value = basis_values * uh(:,n);
        uex = feval(exact_fun, x, t);

        if L_order == 2
            e = e + sum(Gw .* (uh_value - uex).^2) * h(n);
        else
            e = e + sum(Gw .* abs(uh_value - uex).^L_order) * h(n);
        end
    end

    e = e^(1/L_order);
end
