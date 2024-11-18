function result = generate_1D_matrix_r(Gauss_coefficient,mp,trial_drivative,test_derivative)
% evaluate the integral of basis_fun*basis_fun in [-0.5,0.5]
% trial_drivative, test_derivatuve: use to control the derivatuve of the
% basis function

result = zeros(mp+1,mp+1);
GP_r = Gauss_coefficient(:,1);
GW_r = Gauss_coefficient(:,2);
Gpn = size(Gauss_coefficient,1);


for i = 0:mp
    for j = 0:mp
        for k =1:Gpn

            result(i+1,j+1) = result(i+1,j+1) + GW_r(k)...
            *reference_basis(GP_r(k),j,trial_drivative)...
            *reference_basis(GP_r(k),i,test_derivative);

        end
    end
end