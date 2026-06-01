function basis_data = build_basis_data_1D(params)
%BUILD_BASIS_DATA_1D Precompute basis-related data.

    p_order = params.dg.p_order;
    nloc = params.dg.nloc;

    gauss_data = generate_Gauss_lobatto(p_order);
    Gp = gauss_data(:,1);
    Gw = gauss_data(:,2);

    basis_values = zeros(length(Gp), nloc);
    basis_grad_values = zeros(length(Gp), nloc);
    face_values = zeros(2, nloc);
    face_grad_values = zeros(2, nloc);

    for m = 0:p_order
        basis_values(:,m+1) = reference_basis(Gp, m);
        basis_grad_values(:,m+1) = reference_basis_x(Gp, m);

        face_values(1,m+1) = reference_basis(-0.5, m);
        face_values(2,m+1) = reference_basis( 0.5, m);

        face_grad_values(1,m+1) = reference_basis_x(-0.5, m);
        face_grad_values(2,m+1) = reference_basis_x( 0.5, m);
    end

    proj_gauss_data = generate_1D_gaussian_quadrature_r(4);
    proj_Gp = proj_gauss_data(:,1);
    proj_basis_values = zeros(length(proj_Gp), nloc);
    for m = 0:p_order
        proj_basis_values(:,m+1) = reference_basis(proj_Gp, m);
    end

    basis_data.p_order = p_order;
    basis_data.nloc = nloc;

    basis_data.gauss_data = gauss_data;
    basis_data.gauss_points = Gp;
    basis_data.gauss_weights = Gw;
    basis_data.basis_values = basis_values;
    basis_data.basis_grad_values = basis_grad_values;
    basis_data.face_values = face_values;
    basis_data.face_grad_values = face_grad_values;

    basis_data.proj_gauss_data = proj_gauss_data;
    basis_data.proj_points = proj_Gp;
    basis_data.proj_weights = proj_gauss_data(:,2);
    basis_data.proj_basis_values = proj_basis_values;

    basis_data.inv_mass = generate_1D_invmass(p_order);
end