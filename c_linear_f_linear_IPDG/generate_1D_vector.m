function result = generate_1D_vector(f_fun,Gauss_coefficient,ng,mp2,P,T)
% evaluate the right hands term

temp = evaluate_F_multiply_bais(f_fun,Gauss_coefficient,ng,mp2,P,T);

result = reshape(temp,ng*(mp2+1),1);

