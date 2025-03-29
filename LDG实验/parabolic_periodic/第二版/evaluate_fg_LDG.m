function result = evaluate_fg_LDG(Gauss_coefficient,inv_mass,uh,h,mid_points,tm)

[mo,ng] = size(uh);
mp = mo-1;

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

%% 初始化
% fu:对流项
% gu:扩散项
% hu:常数项
fu = zeros(mp+1,ng);
gu = zeros(mp+1,ng);
hu = zeros(mp+1,ng);

u1 = zeros(length(Gp),ng+2);
for n = 2:ng+1
    u1(:,n) = evaluate_uh_DG(uh(:,n-1),Gp,mp);
end
% dirichlet边界处理
n = 2;
u1(2,1) = Dir_boundary(mid_points(n-1)-h(n-1)/2,tm);
n = ng+1;
u1(1,ng+2) = Dir_boundary(mid_points(n-1)+h(n-1)/2,tm);

%% fu = u
for n = 2:ng+1 
    for m = 0:mp

        fu(m+1,n-1) = h(n-1)^(-1)*inv_mass(m+1)*(-Gw'*(f_fun(u1(:,n)).*reference_basis_x(Gp,m))...
                             +1/2*((f_fun(u1(2,n))+f_fun(u1(1,n+1)))*reference_basis(0.5,m) - (f_fun(u1(2,n-1))+f_fun(u1(1,n)))*reference_basis(-0.5,m))...
                             -1/2*((u1(1,n+1)-u1(2,n))*reference_basis(0.5,m) - (u1(1,n)-u1(2,n-1))*reference_basis(-0.5,m)));
        % fu(m+1,n-1) = h(n-1)^(-1)*inv_mass(m+1)*(-Gw'*(f_fun(u1(:,n)).*reference_basis_x(Gp,m))...
        %                                          +f_fun(u1(2,n))*reference_basis(0.5,m) - f_fun(u1(2,n-1))*reference_basis(-0.5,m));
    end
end

%% gu = u
% LDG求解，先求qh再求uh
% (q,v) = (gux,v) = -(gu,vx) + (gu,v)e
qh = zeros(mp+1,ng);
q1 = zeros(length(Gp),ng+2);
for n = 2:ng+1
    for m =  0:mp

        qh(m+1,n-1) = h(n-1)^(-1)*inv_mass(m+1)*(-Gw'*(g_fun(u1(:,n)).*reference_basis_x(Gp,m))...
                             +g_fun(u1(2,n))*reference_basis(0.5,m) - g_fun(u1(2,n-1))*reference_basis(-0.5,m));
    
    end
    q1(:,n) = evaluate_uh_DG(qh(:,n-1),Gp,mp);
end

q1(1,ng+2) = q1(1,2);
% (guxx,v) = (qx,v) = -(q,vx) + (q,v)e
for n = 2:ng+1
    for m = 0:mp
    
        gu(m+1,n-1) =h(n-1)^(-1)*inv_mass(m+1)*(-Gw'*(q1(:,n).*reference_basis_x(Gp,m))...
                   +q1(1,n+1)*reference_basis(0.5,m) - q1(1,n)*reference_basis(-0.5,m));

    end
end

%% hu = cos(t)cos(x) - sin(t)sin(x) + sin(t)cos(x)
for n = 2:ng+1
    Gp_local = Gp*h(n-1) + mid_points(n-1);
    for m = 0:mp
        
        hu(m+1,n-1) = inv_mass(m+1)*Gw'*(h_fun(Gp_local,u1(:,n),tm).*reference_basis(Gp,m));

    end
end

result = gu - fu + hu;
end









