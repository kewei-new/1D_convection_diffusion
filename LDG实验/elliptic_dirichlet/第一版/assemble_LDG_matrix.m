function [LDG_matrix,f_values,Q_matrix]=assemble_LDG_matrix(inv_mass,ng,mp,mid_points,h,H,R,L)
% 把组装矩阵和边界处理放在了一起
alpha = mp/h;

M1 = sparse(ng*(mp+1),ng*(mp+1));
M2 = sparse(ng*(mp+1),ng*(mp+1));

%% 组装(phix,v) = (q,v)
% 得到Q = M2*phi + b2 ，其中B2是含边界项的向量
% 先组装M2
for n = 2:ng
    
    M2((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) = H(:,:,2);
    M2((n-1)*(mp+1)+1:n*(mp+1),(n-2)*(mp+1)+1:(n-1)*(mp+1)) = L(:,:,2);

end
% 组装第1个单元，处理第ng个单元
n = 1;
M2((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) = H(:,:,2);

n = ng;
temp1 = zeros(mp+1,mp+1);
for i = 0:mp
    for j = 0:mp
        temp1(i+1,j+1) = reference_basis(0.5,j)*reference_basis(0.5,i);
    end
end

temp11 = 1/h*diag(inv_mass)*temp1;

M2((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) = M2((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) - temp11;

% 再组装b2
b2 = zeros(ng*(mp+1),1);
n = 1;
g_a = LDG_bound(mid_points(n)-h/2);
for m = 0:mp
    b2((n-1)*(mp+1)+1+m,1) =-1/h*inv_mass(m+1)*g_a*reference_basis(-1/2,m);
end

n = ng;
g_b = LDG_bound(mid_points(n)+h/2);
for m = 0:mp
    b2((n-1)*(mp+1)+1+m,1) = 1/h*inv_mass(m+1)*g_b*reference_basis(1/2,m);
end

%% 组装(qx,v) = (f,v)
% 得到M1*Q + Mn*phi + b1 = F, Mn除了第ng个单元其余单元都是0，是dirichlet为了保证矩阵非奇异加上的
% 先组装M1

for n = 1:ng-1

    M1((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) = H(:,:,1);
    M1((n-1)*(mp+1)+1:n*(mp+1),n*(mp+1)+1:(n+1)*(mp+1)) = R(:,:,1);

end

% 组装第ng个单元
n = ng;
M1((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) = H(:,:,1) + temp1;

% 组装b1
b1 = zeros(ng*(mp+1),1);
n = ng;
for m = 0:mp
    b1((n-1)*(mp+1)+1+m,1) = alpha*g_b*reference_basis(1/2,m);
end

% 组装Mn
Mn = sparse(ng*(mp+1),ng*(mp+1));
n = ng;
Mn((n-1)*(mp+1)+1:n*(mp+1),(n-1)*(mp+1)+1:n*(mp+1)) =-alpha*temp1;

%% 组装总纲
% 得到M1*(M2*phi + b2) + Mn*phi + b1 = F
% 整理可得(M1*M2+Mn)*phi = F - b1 - M1*b2

LDG_matrix = M1*M2 + Mn;
f_values = -b1 - M1*b2;
f_values = reshape(f_values,mp+1,ng);
% 为了计算Q, 把这两部分打包一下返回
Q_matrix = [M2,b2];

end