function [H,R,L] = generate_1D_convection_diffusion_flux_local(inv_matrix_A,matrix_D,matrix_E,theta,ng,mp,P,T,term_type)
% generate the flux from the convection/diffusion term (f(u)_x) / (g(u)_xx)
% [flux(f(u))*v] /[flux(g(u_x))*v]
% consider mass matrix in the same time
H = zeros(mp+1,mp+1,ng);
R = zeros(mp+1,mp+1,ng);
L = zeros(mp+1,mp+1,ng);

% evaluate h
h_local = P(:,T(2,:)) - P(:,T(1,:));

% evaluate 2 to N-1 element
for n = 2:ng-1
    
    if term_type==11
        H(:,:,n) = 1/2*inv_matrix_A*(matrix_D(:,:,1)-matrix_D(:,:,2))/h_local(n);
        R(:,:,n+1) = 1/2*inv_matrix_A*matrix_D(:,:,4)/h_local(n);
        L(:,:,n-1) = -1/2*inv_matrix_A*matrix_D(:,:,3)/h_local(n);

    elseif term_type==12
        H(:,:,n) = inv_matrix_A*(-1/2*(matrix_E(:,:,1)-matrix_E(:,:,2))/h_local(n)+theta*(matrix_D(:,:,1)+matrix_D(:,:,2)))/h_local(n);
        R(:,:,n+1) = inv_matrix_A*(-1/2*matrix_E(:,:,4)/h_local(n+1)-theta*matrix_D(:,:,4))/h_local(n);
        L(:,:,n-1) = inv_matrix_A*(1/2*matrix_E(:,:,3)/h_local(n-1)-theta*matrix_D(:,:,3))/h_local(n);
    end


end

% evaluate 1,N element
% 1
n=1;
if term_type==11
    H(:,:,n) = 1/2*inv_matrix_A*(matrix_D(:,:,1)-matrix_D(:,:,2))/h_local(n);
    R(:,:,n+1) = 1/2*inv_matrix_A*matrix_D(:,:,4)/h_local(n);
elseif term_type==12
    H(:,:,n) = inv_matrix_A*(-1/2*(matrix_E(:,:,1)-matrix_E(:,:,2))/h_local(n)+theta*(matrix_D(:,:,1)+matrix_D(:,:,2)))/h_local(n);
    R(:,:,n+1) = inv_matrix_A*(-1/2*matrix_E(:,:,4)/h_local(n+1)-theta*matrix_D(:,:,4))/h_local(n);
end
% N
n=ng;
if term_type==11
    H(:,:,n) = 1/2*inv_matrix_A*(matrix_D(:,:,1)-matrix_D(:,:,2))/h_local(n);
    L(:,:,n-1) = -1/2*inv_matrix_A*matrix_D(:,:,3)/h_local(n);
elseif term_type==12
    H(:,:,n) = inv_matrix_A*(-1/2*(matrix_E(:,:,1)-matrix_E(:,:,2))/h_local(n)+theta*(matrix_D(:,:,1)+matrix_D(:,:,2)))/h_local(n);
    L(:,:,n-1) = inv_matrix_A*(1/2*matrix_E(:,:,3)/h_local(n-1)-theta*matrix_D(:,:,3))/h_local(n);
end
