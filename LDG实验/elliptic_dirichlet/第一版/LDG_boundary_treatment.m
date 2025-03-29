function [LDG_matrix,f_values] = LDG_boundary_treatment(LDG_matrix,inv_mass,mp,ng,H,R,L,h,mid_points)
% 这个部分放在外面是为了节省开销
% 边界处flux取值如下：
% phi(0-) = g_a
% Q(N+) = Q(N-) - c*(phi(N-) - g_b)

f_values = zeros(mp+1,ng);
M = diag(inv_mass)*h^(-1);
b1 = zeros(mp+1,1);
b2 = zeros(mp+1,1);

% 罚项系数
alpha = mp/h;

% 求(g_a,v(-1/2)) (g_b,v(1/2))
for i = 0:mp

    b1(i+1,1) = LDG_bound(mid_points(1)-h/2)*reference_basis(-1/2,i);
    b2(i+1,1) = LDG_bound(mid_points(ng)+h/2)*reference_basis(1/2,i);

end

% 求(u(1/2),v(1/2))
M1 = zeros(mp+1,mp+1);
for i = 0:mp
    for j = 0:mp

        M1(i+1,j+1) = reference_basis(1/2,j)*reference_basis(1/2,i);

    end
end


% n==1
n = 1;
LDG_matrix((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) = H(:,:,1)*H(:,:,2) + R(:,:,1)*L(:,:,2);
LDG_matrix((n-1)*(mp+1)+1:n*(mp+1),(n)*(mp+1)+1:(n+1)*(mp+1)) = R(:,:,1)*H(:,:,2);
f_values(:,n) = H(:,:,1)*M*b1;

% n==ng
n = ng;
H(:,:,2) = H(:,:,2)-M*M1; %第N单元的右边界用g_b替换
LDG_matrix((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) = H(:,:,1)*H(:,:,2) + M1*H(:,:,2) - alpha*M1;
LDG_matrix((n-1)*(mp+1)+1:n*(mp+1),(n-2)*(mp+1)+1:(n-1)*(mp+1)) = H(:,:,1)*L(:,:,2) + M1*L(:,:,2);
f_values(:,n) = -alpha*b2-H(:,:,1)*M*b2-M1*M*b2;
