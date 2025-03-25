function result = generate_1D_vector_FEM(f_fun,uh,Gauss_coefficient,P,T,Tb,basis_type)
% Solving the right term using FEM
% The parameter uh needs to be judged on the basis of the question as to whether it needs to be entered or not, and here it needs to be
result = evaluate_F_multiply_bais_FEM(f_fun,uh,Gauss_coefficient,P,T,Tb,0,basis_type);
