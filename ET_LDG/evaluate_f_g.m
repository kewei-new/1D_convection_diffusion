function result = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1,uh2,mp,h,mid_points,xm,xe,xk,T_0,epsilon,emf)
% 本函数考虑的是第n个有限元上的情况
% 基础变量名称:u1,u2,phix,k,E,T,muE,muEE,DE,DEE
% 后缀带x，例如Ex,Tx，表示一阶求导
% 后缀为lbc,表示左侧单元的端点值，一般出现是由于flux取负，会出现u^{n}(-1-)v(-1+)=u^{n-1}( 1-)v(-1+)
% 后缀为rbc,表示右侧单元的端点值，一般出现是由于flux取正，会出现u^{n}( 1+)v( 1-)=u^{n+1}(-1+)v( 1-)
% 后缀为rl,lr分别表示右侧单元的左端点，以及左侧单元的右端点，主要是懒得两个一起输入了
% 注：(1)Ec变量由于只出现在h(u)之中，没写函数名;
% (2)lbc和rbc是个二元组，第一个位置是左端点-0.5处的值，第二个位置是右端点0.5处的值
ng = size(uh1,2) - 2;
f = zeros(mp+1,ng+2,2);% f(:,1)==f1,f(:,2)==f2
g = zeros(mp+1,ng+2,2);% g(:,1)==g1,g(:,2)==g2
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);


u1 = ones(length(Gp),ng+2);
u2 = ones(length(Gp),ng+2);
phix = ones(length(Gp),ng+2);
xmu0 = zeros(length(Gp),ng+2);
for k = 1:ng+2

   u1(:,k) = evaluate_uh_DG(Gp,uh1(:,k),mp);
   u2(:,k) = evaluate_uh_DG(Gp,uh2(:,k),mp);
   phix(:,k) = evaluate_uh_DG(Gp,phixh(:,k),mp);
   xmu0(:,k) = mu0(Gp*h+mid_points(k));
   
end

u1(2,1) = init_fun1(0);
u2(2,1) = init_fun2(0);
phix(2,1) = phix(1,2);

u1(1,ng+2) = init_fun1(0.6);
u2(1,ng+2) = init_fun2(0.6);
phix(1,ng+2) = phix(2,ng+1);



% 已检查
n = u1/xe;
E = u2./n * xm;
T = 2./5./xk*(sqrt(abs(1.+10.*E/3.))-1.);

% 已检查
muE = xmu0*T_0./T;
muEE = 3/2*xmu0*xk*T_0.*(1-5/4*xk*T/xe);
DE = xk*xmu0*T_0;
DEE = 3/2*xmu0*xk^2*T_0.*T/xe.*(1-5/4*xk/xe*T);

% 已检查
% evaluate -(f1,vx) + (f1,v)_e
% f1 = phix.*u1.*muE

f1 = phix.*u1.*muE;

for k = 2:ng+1
    for m = 0:mp
    
        f(m+1,k,1) = inv_mass(m+1)*(-Gw'*(f1(:,k).*reference_basis_x(Gp,m))...
                  +1/2*((f1(2,k)+f1(1,k+1))*reference_basis(0.5,m)-(f1(2,k-1)+f1(1,k))*reference_basis(-0.5,m))...
                  -1/2*emf*((u1(1,k+1)-u1(2,k))*reference_basis(0.5,m)-(u1(1,k)-u1(2,k-1))*reference_basis(-0.5,m)))*h^(-1);

        % f(m+1,k,1) = inv_mass(m+1)*(-Gw'*(f1(:,k).*reference_basis_x(Gp,m))...
        %           +1/2*((f1(2,k)+f1(1,k+1))*reference_basis(0.5,m)-(f1(2,k-1)+f1(1,k))*reference_basis(-0.5,m))...
        %           )*h^(-1);
              
    end
end

% 已检查
% evaluate -(f2,vx) + (f2,v)_e
% f2 = phix.*n.*(muEE+DE)

f2 = phix.*n.*(muEE+DE);

for k = 2:ng+1
    for m = 0:mp
    
        f(m+1,k,2) = inv_mass(m+1)*(-Gw'*(f2(:,k).*reference_basis_x(Gp,m))...
                  +1/2*((f2(2,k)+f2(1,k+1))*reference_basis(0.5,m)-(f2(2,k-1)+f2(1,k))*reference_basis(-0.5,m))...
                  -1/2*emf*((u2(1,k+1)-u2(2,k))*reference_basis(0.5,m)-(u2(1,k)-u2(2,k-1))*reference_basis(-0.5,m)))*h^(-1);
       
        % f(m+1,k,2) = inv_mass(m+1)*(-Gw'*(f2(:,k).*reference_basis_x(Gp,m))...
        %           +1/2*((f2(2,k)+f2(1,k+1))*reference_basis(0.5,m)-(f2(2,k-1)+f2(1,k))*reference_basis(-0.5,m))...
        %           )*h^(-1);
              
    end
end

% evaluate -(g1xx,v)
% g1 = DE.*n
qh = zeros(mp+1,ng,2);

q1 = zeros(length(Gp),ng+1);
g1 = DE.*n;

for k = 2:ng+1
    for m = 0:mp
    
        qh(m+1,k-1,1) =inv_mass(m+1)*(-Gw'*(g1(:,k).*reference_basis_x(Gp,m))...
              +g1(2,k)*reference_basis(0.5,m)-g1(2,k-1)*reference_basis(-0.5,m))*h^(-1);
    
    end
    q1(:,k-1) = evaluate_uh_DG(Gp, qh(:,k-1,1),mp);
end
q1(1,ng+1) = q1(2,ng);

for k =1:ng
    for m = 0:mp
    
        g(m+1,k+1,1) =inv_mass(m+1)*(-Gw'*(q1(:,k).*reference_basis_x(Gp,m))...
              -q1(1,k)*reference_basis(-0.5,m)+q1(1,k+1)*reference_basis(0.5,m))*h^(-1);
    
    end
end

% evaluate -(g2x,vx) + (g2x,v)e
% g2 = DEE.*n
g2 = DEE.*n;
q2 = zeros(length(Gp),ng+1);

for k = 2:ng+1
    for m = 0:mp
    
        qh(m+1,k-1,2) =inv_mass(m+1)*(-Gw'*(g2(:,k).*reference_basis_x(Gp,m))...
              +g2(2,k)*reference_basis(0.5,m)-g2(2,k-1)*reference_basis(-0.5,m))*h^(-1);
    
    end
    q2(:,k-1) = evaluate_uh_DG(Gp, qh(:,k-1,2),mp);
end
q2(1,ng+1) = q2(2,ng);

for k =1:ng
    for m = 0:mp
    
        g(m+1,k+1,2) =inv_mass(m+1)*(-Gw'*(q2(:,k).*reference_basis_x(Gp,m))...
              +q2(1,k+1)*reference_basis(0.5,m)-q2(1,k)*reference_basis(-0.5,m))*h^(-1);
    
    end
end


hu = zeros(mp+1,ng+2);

Ec = 15/4*xk.*((1+1/2*xk/xe*T).*T-T_0);
for k = 2:ng+1
    for m = 0:mp
    
        hu(m+1,k) = inv_mass(m+1)*(Gw'*((xe*n(:,k).*muE(:,k).*phix(:,k).^2 ...
                        +xe/epsilon*(n(:,k)-nd(Gp*h+mid_points(k))).*n(:,k).*DE(:,k)...
                        -n(:,k).*Ec(:,k)).*reference_basis(Gp,m)));
    
    end
end


result = zeros(mp+1,ng+2,2);
result(:,:,1) = g(:,:,1) - f(:,:,1);
result(:,:,2) = g(:,:,2) - f(:,:,2) + hu;
