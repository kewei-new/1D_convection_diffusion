function uh0 = L2_projection_1D(initial_fun,mass_M_r,ng,mp,P,T)
% evaluate discrete initial condition

Gauss_coefficient = generate_1D_gaussian_quadrature_r(4);

right_term = evaluate_F_multiply_bais(initial_fun,Gauss_coefficient,ng,mp,P,T);

% 左右两端都会出一个h_lcoal,相互抵消了
uh0_temp = mass_M_r\right_term;

% 需要拉成长条
uh0 = reshape(uh0_temp,ng*(mp+1),1);
