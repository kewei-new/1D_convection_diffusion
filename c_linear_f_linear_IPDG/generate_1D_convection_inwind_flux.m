function [H,R,L]=generate_1D_convection_inwind_flux(inv_matrix_A,matrix_D,ng,mp,P,T)

H = zeros(mp+1,mp+1,ng);
R = zeros(mp+1,mp+1,ng);
L = zeros(mp+1,mp+1,ng);

% evaluate h
h_local = P(:,T(2,:)) - P(:,T(1,:));

% evaluate 2 to N-1 element
for n = 2:ng
    
    H(:,:,n) = inv_matrix_A*matrix_D(:,:,1)/h_local(n);
    L(:,:,n-1) = -inv_matrix_A*matrix_D(:,:,3)/h_local(n);

end

n = 1;
H(:,:,n) = inv_matrix_A*matrix_D(:,:,1)/h_local(n);
L(:,:,ng) = -inv_matrix_A*matrix_D(:,:,3)/h_local(n);
