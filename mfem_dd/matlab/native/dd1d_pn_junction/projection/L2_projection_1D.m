function uh0 = L2_projection_1D(initial_fun, arg2, arg3, arg4, arg5, arg6)
%L2_PROJECTION_1D Evaluate the discrete initial condition.
%
% New call:
%   uh0 = L2_projection_1D(initial_fun, basis_data, mesh)
%
% Backward-compatible old call:
%   uh0 = L2_projection_1D(initial_fun, inv_a, ng, mp, h, mid_points)

    if nargin == 3 && isstruct(arg2) && isstruct(arg3)
        basis_data = arg2;
        mesh = arg3;

        inv_a = basis_data.inv_mass(:);
        Gp = basis_data.proj_points;
        Gw = basis_data.proj_weights;
        basis_values = basis_data.proj_basis_values;

        x = Gp * mesh.cell_sizes + mesh.cell_centers;
        ini_value = feval(initial_fun, x);
        projected = basis_values' * bsxfun(@times, Gw, ini_value);
        projected = bsxfun(@times, inv_a, projected);
        uh0 = projected(:);
        return;
    end

    inv_a = arg2;
    ng = arg3;
    mp = arg4;
    h = arg5;
    mid_points = arg6;

    Gauss_coefficient = generate_1D_gaussian_quadrature_r(4);
    Gp = Gauss_coefficient(:,1);
    Gw = Gauss_coefficient(:,2);

    uh0 = zeros(mp+1,ng);
    for n = 1:ng
        Gp_local = Gp*h(n) + mid_points(n);
        ini_value = feval(initial_fun,Gp_local);

        for m = 0:mp
            uh0(m+1,n) = inv_a(m+1) * Gw' * (ini_value .* reference_basis(Gp,m));
        end
    end

    uh0 = reshape(uh0,(mp+1)*ng,1);
end
