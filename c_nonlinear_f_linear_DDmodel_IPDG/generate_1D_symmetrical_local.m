function [H,R,L]=generate_1D_symmetrical_local(inv_matrix_A,matrix_F,ng,mp,P,T)
% 生成对称项
% 但这个函数是有点冗余的，后续可以考虑如何更改
% 先写出来，然后对比一下看看是不是对称的

H = zeros(mp+1,mp+1,ng);
R = zeros(mp+1,mp+1,ng);
L = zeros(mp+1,mp+1,ng);

h_local = P(:,T(2,:)) - P(:,T(1,:));

for n = 2:ng-1

        H(:,:,n) = inv_matrix_A*(1/2*(matrix_F(:,:,1)-matrix_F(:,:,2))/h_local(n)/h_local(n));
        R(:,:,n) = inv_matrix_A*(-1/2*matrix_F(:,:,4)/h_local(n-1))/h_local(n);
        L(:,:,n) = inv_matrix_A*(1/2*matrix_F(:,:,3)/h_local(n+1))/h_local(n);

end

n = 1;
H(:,:,n) = inv_matrix_A*(1/2*(matrix_F(:,:,1)-matrix_F(:,:,2))/h_local(n)/h_local(n));
L(:,:,n) = inv_matrix_A*(1/2*matrix_F(:,:,3)/h_local(n+1))/h_local(n);
R(:,:,n) = inv_matrix_A*(-1/2*matrix_F(:,:,4)/h_local(ng))/h_local(n);

n = ng;
H(:,:,n) = inv_matrix_A*(1/2*(matrix_F(:,:,1)-matrix_F(:,:,2))/h_local(n)/h_local(n));
R(:,:,n) = inv_matrix_A*(-1/2*matrix_F(:,:,4)/h_local(n-1))/h_local(n);
L(:,:,n) = inv_matrix_A*(1/2*matrix_F(:,:,3)/h_local(1))/h_local(n);