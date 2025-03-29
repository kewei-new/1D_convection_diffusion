function [uh,qh]=solve_1D_LDG(Gauss_coefficient,inv_mass,T,h,mid_points,mp,ng)
% inv_mass是质量矩阵的逆
% ut + ux = uxx + f

% u = sin(t)cos(x), f = cos(t)cos(x) - sin(t)sin(x) + sin(t)cos(x)
% u = exp(-t)sin(x-t), f = 0

% u(0) = exp(-t)sin(t) u(2*pi) = exp(-t)sin(t)
% u(0) = sin(t) u(2*pi)=sin(t);

% ut + q = qx + f
% q = ux
% (ut,v) + (q,v) = -(q,vx) + (q,v)e + (f,v)
% (q,p) = -(u,px) + (u,p)e

%% generate matrix
emf = 1;
emg = 1;

%% L2 projection
uh = zeros(mp+1,ng,3);
uh(:,:,1) = L2_projection_1D(inv_mass,ng,mp,h,mid_points);

%% time_evolution
% result 存储第一个时间步的uh
% result1 存储第二个时间步的uh
qh = zeros(mp+1,ng);

% tic
t = 0;
dt = evaluate_dt(emf,emg,h,mp);
% plot_DG('exact_fun',0,uh(:,:,1),ng,mp,mid_points)
while t < T
    
    if t+dt>T
        dt = T-t;
    end
    % print_progress(t, T);
    for io = 1:3
        uh = RK3_LDG(Gauss_coefficient,inv_mass,uh,h,mid_points,t,dt,io);
    end

    t = t+dt;
    % plot_DG('exact_fun',t,uh(:,:,1),ng,mp,mid_points)
end

uh = uh(:,:,1);

% 求qh
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
u1 = zeros(length(Gp),ng+2);
for n = 2:ng+1
    u1(:,n) = evaluate_uh_DG(uh(:,n-1),Gp,mp);
end
% dirichlet边界处理
n = 2;
u1(2,1) = Dir_boundary(mid_points(n-1)-h(n-1)/2,T);
n = ng+1;
u1(1,ng+2) = Dir_boundary(mid_points(n-1)+h(n-1)/2,T);

qh = zeros(mp+1,ng);
for n = 2:ng+1
    for m =  0:mp

        qh(m+1,n-1) = h(n-1)^(-1)*inv_mass(m+1)*(-Gw'*(g_fun(u1(:,n)).*reference_basis_x(Gp,m))...
                             +g_fun(u1(2,n))*reference_basis(0.5,m) - g_fun(u1(2,n-1))*reference_basis(-0.5,m));
    
    end
end
