function generate_1D_matrix_with_coefficient(coefficient_fun,phih_local,Gauss_coefficient,mp1,mp2,trial_basis_type,trial_drivative,test_basis_type,test_derivative)

result = sparse(mp2+1,mp1+1);
GP_r = Gauss_coefficient(:,1);
GW_r = Gauss_coefficient(:,2);
Gpn = size(Gauss_coefficient,1);


for i = 0:mp2
    for j = 0:mp1
        for k =1:Gpn

            result(i+1,j+1) = result(i+1,j+1) + GW_r(k)...
            *feval(coefficient_fun,GP_r(k),phih_local,mp1)...
            *reference_basis(GP_r(k),j,trial_basis_type,trial_drivative)...
            *reference_basis(GP_r(k),i,test_basis_type,test_derivative);

        end
    end
end