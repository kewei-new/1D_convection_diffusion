function [H,R,L]=generate_1D_penalty_local(inv_matrix_A,matrix_D,theta,ng,mp,P,T)

H = zeros(mp+1,mp+1,ng);
R = zeros(mp+1,mp+1,ng);
L = zeros(mp+1,mp+1,ng);

% evaluate h
h_local = P(:,T(2,:)) - P(:,T(1,:));

% evaluate 2 to N-1 element
for n = 2:ng-1
    
    H(:,:,n) = inv_matrix_A*(theta*(matrix_D(:,:,1)+matrix_D(:,:,2)))/h_local(n);
    R(:,:,n+1) = inv_matrix_A*(-theta*matrix_D(:,:,4))/h_local(n);
    L(:,:,n-1) = inv_matrix_A*(-theta*matrix_D(:,:,3))/h_local(n);

end

% evaluate 1,N element
% 1
n=1;
H(:,:,n) = inv_matrix_A*(theta*(matrix_D(:,:,1)+matrix_D(:,:,2)))/h_local(n);
R(:,:,n+1) = inv_matrix_A*(-theta*matrix_D(:,:,4))/h_local(n);
L(:,:,ng) = inv_matrix_A*(-theta*matrix_D(:,:,3))/h_local(n);
% N
n=ng;
H(:,:,n) = inv_matrix_A*(theta*(matrix_D(:,:,1)+matrix_D(:,:,2)))/h_local(n);
R(:,:,1) = inv_matrix_A*(-theta*matrix_D(:,:,4))/h_local(n);
L(:,:,n-1) = inv_matrix_A*(-theta*matrix_D(:,:,3))/h_local(n);