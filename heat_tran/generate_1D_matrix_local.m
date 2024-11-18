function M = generate_1D_matrix_local(Gauss_coefficient,inv_a,mp,h)

M1 = generate_1D_matrix_r(Gauss_coefficient,mp,1,1);

M = inv_a*M1*h^(-2);
