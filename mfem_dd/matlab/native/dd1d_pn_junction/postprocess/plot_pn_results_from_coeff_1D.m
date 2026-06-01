function plot_pn_results_from_coeff_1D(result, mesh, basis_data, params, tag, outdir)

    x = mesh.cell_centers(:);

    phi = evaluate_DG_on_centers(result.phi_coeff, mesh.num_cells, basis_data.p_order);
    E   = -evaluate_DG_on_centers(result.E_coeff,   mesh.num_cells, basis_data.p_order);
    n   = evaluate_DG_on_centers(result.n_coeff,    mesh.num_cells, basis_data.p_order);
    doping = nd_profile(x, params);

    plot_pn_results_1D(x, doping, n, phi, E, outdir, tag);
end