function [L,H,R] = generate_1D_variable_coefficient_flux_local(inv_matrix_A,matrix_D,phi,mp,ng,P,T,Tb)

H = zeros(mp+1,mp+1,ng);
R = zeros(mp+1,mp+1,ng);
L = zeros(mp+1,mp+1,ng);

% evaluate h
h_local = P(:,T(2,:)) - P(:,T(1,:));
phi = phi(Tb);
% evaluate 2 to N-1 element
for n = 2:ng-1
    
    H(:,:,n) = 1/2*inv_matrix_A*(matrix_D(:,:,1)*evaluate_local_FEM(0.5,phi(:,n),101,1)-matrix_D(:,:,2)*evaluate_local_FEM(-0.5,phi(:,n),101,1))/h_local(n)^2;
    R(:,:,n+1) = 1/2*inv_matrix_A*matrix_D(:,:,4)*evaluate_local_FEM(-0.5,phi(:,n+1),101,1)/h_local(n)/h_local(n+1);
    L(:,:,n-1) = -1/2*inv_matrix_A*matrix_D(:,:,3)*evaluate_local_FEM(0.5,phi(:,n-1),101,1)/h_local(n)/h_local(n-1);


end

n = 1;
H(:,:,n) = 1/2*inv_matrix_A*(matrix_D(:,:,1)*evaluate_local_FEM(0.5,phi(:,n),101,1)-matrix_D(:,:,2)*evaluate_local_FEM(-0.5,phi(:,n),101,1))/h_local(n)^2;
R(:,:,n+1) = 1/2*inv_matrix_A*matrix_D(:,:,4)*evaluate_local_FEM(-0.5,phi(:,n+1),101,1)/h_local(n)/h_local(n+1);
L(:,:,ng) = -1/2*inv_matrix_A*matrix_D(:,:,3)*evaluate_local_FEM(0.5,phi(:,ng),101,1)/h_local(n)/h_local(ng);

n = ng;
H(:,:,n) = 1/2*inv_matrix_A*(matrix_D(:,:,1)*evaluate_local_FEM(0.5,phi(:,n),101,1)-matrix_D(:,:,2)*evaluate_local_FEM(-0.5,phi(:,n),101,1))/h_local(n)^2;
L(:,:,n-1) = -1/2*inv_matrix_A*matrix_D(:,:,3)*evaluate_local_FEM(0.5,phi(:,n-1),101,1)/h_local(n)/h_local(n-1);
R(:,:,1) = 1/2*inv_matrix_A*matrix_D(:,:,4)*evaluate_local_FEM(-0.5,phi(:,1),101,1)/h_local(n)/h_local(1);