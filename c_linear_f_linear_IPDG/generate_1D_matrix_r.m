function result = generate_1D_matrix_r(Gauss_coefficient,mp1,mp2,trial_basis_type,trial_drivative,test_basis_type,test_derivative)
% evaluate the integral of basis_fun*basis_fun in [-0.5,0.5]
% trial_drivative, test_derivatuve: use to control the derivatuve of the
% basis function

result = sparse(mp2+1,mp1+1);
GP_r = Gauss_coefficient(:,1);
GW_r = Gauss_coefficient(:,2);
Gpn = size(Gauss_coefficient,1);


for i = 0:mp2
    for j = 0:mp1
        for k =1:Gpn

            result(i+1,j+1) = result(i+1,j+1) + GW_r(k)...
            *reference_basis(GP_r(k),j,trial_basis_type,trial_drivative)...
            *reference_basis(GP_r(k),i,test_basis_type,test_derivative);

        end
    end
end