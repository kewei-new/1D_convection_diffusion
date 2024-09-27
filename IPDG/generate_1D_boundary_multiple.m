function result = generate_1D_boundary_multiple(mp,trial_basis_type,trial_derivative,test_basis_type)
% evaluate basis fun mutiply on the boundary, such as
% u(0.5)*u(0.5),u(-0.5)*u(-0.5),u(0.5)*u(-0.5)
% because in the DG method, in the interior boundary,usually involves
% u(right_point)^{-}*u(right_point)^{-},u(left_point)^{+}*u(left_point)^+},u(right_point)^{-}*u(left_point)^{+}
% Abbreviate them to brr,bll,brl
% result(:,:,1)==brr,result(:,:,2)==bll,result(:,:,3)==brl,result(:,:,4)==blr

result = zeros(mp+1,mp+1,4);

for i = 0:mp

    for j = 0:mp
    
        result(i+1,j+1,1) = reference_basis(0.5,j,trial_basis_type,trial_derivative)*reference_basis(0.5,i,test_basis_type,0);
        result(i+1,j+1,2) = reference_basis(-0.5,j,trial_basis_type,trial_derivative)*reference_basis(-0.5,i,test_basis_type,0);
        result(i+1,j+1,3) = reference_basis(0.5,j,trial_basis_type,trial_derivative)*reference_basis(-0.5,i,test_basis_type,0);
        result(i+1,j+1,4) = reference_basis(-0.5,j,trial_basis_type,trial_derivative)*reference_basis(0.5,i,test_basis_type,0);

    end
end


end