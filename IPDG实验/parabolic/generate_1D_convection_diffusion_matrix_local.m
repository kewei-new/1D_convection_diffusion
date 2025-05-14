function result = generate_1D_convection_diffusion_matrix_local(inv_matrix_A,matrix_B,matrix_C,ng,mp1,mp2,P,T,term_type)
% generate related convection/diffusion integral in local space
% matrix_A:mass matrix generate in[-0.5,0.5]
% matrix_B:stiffiness matrix generate in[-0.5,0.5]
% matrix_C:(uh,uhx) matrix generate in[-0.5,0.5]
% term_type:convection or diffusion terms; convection=11, diffusion=12


% in any element(intoa ng elements),it's a mp2*mp1(test*trial) Matrix
result = zeros(mp2+1,mp1+1,ng);

for n = 1:ng
    
    local_points = P(:,T(:,n));
    h_local = max(local_points) - min(local_points);

    if term_type == 11

        result(:,:,n) = inv_matrix_A*matrix_C/h_local;

    elseif term_type == 12
    
        result(:,:,n) = inv_matrix_A*matrix_B/(h_local)^2;
    end

end

