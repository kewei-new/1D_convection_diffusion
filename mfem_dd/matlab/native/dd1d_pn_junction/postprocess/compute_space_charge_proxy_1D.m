function charge_data = compute_space_charge_proxy_1D(n_coeff, mesh, basis_data, params)
%COMPUTE_SPACE_CHARGE_PROXY_1D Quasi-static space-charge and depletion proxies.
%
% Qmag = 0.5 * integral |n - Nd|
% Wdep = measure of cells where |n-Nd| exceeds a relative threshold.

    n_coeff = reshape(n_coeff, basis_data.nloc, mesh.num_cells);
    x = basis_data.gauss_points * mesh.cell_sizes + mesh.cell_centers;
    n_val = basis_data.basis_values * n_coeff;
    Nd = params.problem.source.poisson(x, 0);
    rho = Nd - n_val;

    Qmag = 0.5 * sum(sum(abs(rho) .* (basis_data.gauss_weights * mesh.cell_sizes)));

    rho_cell = mean(abs(rho), 1);
    threshold = 0.05 * max(rho_cell);
    active = rho_cell > max(threshold, 1e-12);
    if any(active)
        Wdep = sum(mesh.cell_sizes(active));
    else
        Wdep = 0;
    end

    charge_data.Qmag = Qmag;
    charge_data.Wdep = Wdep;
    charge_data.rho_cell = rho_cell(:);
end
