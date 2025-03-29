function [phi,E]=test(inv_mass,ng,mid_points,h,mp)
u0=1;u1=-1;
c0=mp/h;
Gauss_coefficient=generate_Gauss_lobatto(mp); % 六点高斯积分
gauss_p = Gauss_coefficient(:,1);
gauss_c = Gauss_coefficient(:,2);
gn=length(gauss_p);
%% A逆
% 质量矩阵的逆
M=sparse((mp+1)*ng,(mp+1)*ng); 
A = diag(inv_mass)*1/h;
for i=1:ng
    M((i-1)*(mp+1)+1:i*(mp+1),(i-1)*(mp+1)+1:i*(mp+1))=A;
end
%% B
B=zeros((mp+1)*ng,(mp+1)*ng);
B1=zeros(mp+1,mp+1);% (u,vx)
for k=1:mp+1
    for j=1:mp+1
        for z=1:gn
            B1(k,j)=B1(k,j)+gauss_c(z)*reference_basis(gauss_p(z),j-1)*reference_basis_x(gauss_p(z),k-1);
        end
    end
end
B2=zeros(mp+1,mp+1);% u(1)v(1)
for k=1:mp+1
    for j=1:mp+1
        B2(k,j)=reference_basis(1/2,j-1)*reference_basis(1/2,k-1);
    end
end    
B3=zeros(mp+1,mp+1);% u(1)v(-1)
for k=1:mp+1
    for j=1:mp+1
        B3(k,j)=reference_basis(1/2,j-1)*reference_basis(-1/2,k-1);
    end
end 
% 转载矩阵 解决(phix,p) = (Q,p)
for i=1:ng  %B(i,j) = -(u,vx)
    B((i-1)*(mp+1)+1:i*(mp+1),(i-1)*(mp+1)+1:i*(mp+1))=-B1;
end
for i=1:ng-1 %除了最右单元  B(i,j) = -(u,vx) + u(1)v(1)
    B((i-1)*(mp+1)+1:i*(mp+1),(i-1)*(mp+1)+1:i*(mp+1))=B((i-1)*(mp+1)+1:i*(mp+1),(i-1)*(mp+1)+1:i*(mp+1))+B2;
end
for i=2:ng %装载左侧单元矩阵  B(i,j-1) = -u(1)v(-1)
    B((i-1)*(mp+1)+1:i*(mp+1),(i-2)*(mp+1)+1:(i-1)*(mp+1))=B((i-1)*(mp+1)+1:i*(mp+1),(i-2)*(mp+1)+1:(i-1)*(mp+1))-B3;
end
% B*phi = Q
B=M*B;
%% U0
% 边界项处理 相当于g_a*v(0+),g_b*v(N-)
U0=zeros((mp+1)*ng,1);
for k=1:mp+1
    U0(k,1)=-u0*reference_basis(-1/2,k-1);
    U0(k+(mp+1)*(ng-1),1)=u1*reference_basis(1/2,k-1);
end
U0=M*U0;
%% C
% Q=M*BU+M*U0;
% CQ+DU=F; C(M*BU+M*U0)+DU=F; U Q
C=zeros((mp+1)*ng,(mp+1)*ng);
C1=zeros(mp+1,mp+1);
for k=1:mp+1
    for j=1:mp+1
        C1(k,j)=reference_basis(-1/2,k-1)*reference_basis(-1/2,j-1); %u(-1)v(-1)
    end
end 
for i=1:ng
    C((i-1)*(mp+1)+1:i*(mp+1),(i-1)*(mp+1)+1:i*(mp+1))=B1+C1; % C(i,j) = u(-1)v(-1) + (u,vx)
end
for i=1:ng-1 %C(i,j+1) = -u(-1)v(1)
    C((i-1)*(mp+1)+1:i*(mp+1),i*(mp+1)+1:(i+1)*(mp+1))=C((i-1)*(mp+1)+1:i*(mp+1),i*(mp+1)+1:(i+1)*(mp+1))-B3';
end
% C(N,N) = u(-1)v(-1) + (u,vx) - u(1)v(1)
C(1+(mp+1)*(ng-1):(mp+1)*ng,1+(mp+1)*(ng-1):(mp+1)*ng)=C(1+(mp+1)*(ng-1):(mp+1)*ng,1+(mp+1)*(ng-1):(mp+1)*ng)-B2;
% D*phi = C*B*phi = C*Q
D=C*B;

% U1 = C*U0 = C*M*U0 = C*M*u1*ref_basis(0.5)
U1=C*U0;
%%
% D(N,N) = C*B(N,N) + m/h*u(1)v(1)
D(1+(mp+1)*(ng-1):(mp+1)*ng,1+(mp+1)*(ng-1):(mp+1)*ng)=D(1+(mp+1)*(ng-1):(mp+1)*ng,1+(mp+1)*(ng-1):(mp+1)*ng)+c0*B2;
for k=1:mp+1
    % U1(N) = C*M*u1*ref_basis(0.5) - m/h*u1*ref_basis(0.5)
    U1(k+(mp+1)*(ng-1),1)=U1(k+(mp+1)*(ng-1),1)-c0*u1*reference_basis(1/2,k-1);
end
%% F
F=zeros((mp+1)*ng,1);
 for j=1:ng    
    for k=1:mp+1
        for z=1:gn
            % F(i) = (fi,v)
            F(k+(mp+1)*(j-1),1)=F(k+(mp+1)*(j-1),1)+gauss_c(z)*(-cos(h*gauss_p(z)+mid_points(j)))*reference_basis(gauss_p(z),k-1);
        end
    end
 end


%DU = F - U1
U=D\(h*F-U1);
uh=zeros(ng+1,1);
for i=1:ng
    for j=1:mp+1
        uh(i,1)=uh(i,1)+U((i-1)*(mp+1)+j,1)*reference_basis(-1/2,j-1);
    end
end
for j=1:mp+1
    uh(ng+1,1)=uh(ng+1,1)+U((ng-1)*(mp+1)+j,1)*reference_basis(1/2,j-1);
end
phi=reshape(U,mp+1,ng);
E=reshape(B*U+U0,mp+1,ng);