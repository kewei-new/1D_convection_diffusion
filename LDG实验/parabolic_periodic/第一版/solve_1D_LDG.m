function [uh,qh]=solve_1D_LDG(Gauss_coefficient,inv_mass,T,h,mid_points,mp,ng)
% inv_mass是质量矩阵的逆
% ut + ux = uxx + f

% u = sin(x-t), f = sin(x-t)
% u = exp(-t)sin(x-t), f = 0

% ut + q = qx + f
% q = ux
% (ut,v) + (q,v) = -(q,vx) + (q,v)e + (f,v)
% (q,p) = -(u,px) + (u,p)e

%% generate matrix
emf = 1;
emg = 1;
% 矩阵大小为 (*,*,2)
[H,R,L] = generate_1D_matrix(Gauss_coefficient,mp,inv_mass,h);

%% L2 projection
uh = L2_projection_1D(inv_mass,ng,mp,h,mid_points);
uh_m = [uh(:,ng),uh,uh(:,1)];

%% time_evolution
% result 存储第一个时间步的uh
% result1 存储第二个时间步的uh
result = zeros(mp+1,ng+2);
result1 = zeros(mp+1,ng+2);
qh = zeros(mp+1,ng);


% tic
t = 0;
f_term = generate_1D_vector(Gauss_coefficient,ng,mp,h,mid_points,inv_mass,t);
dt = evaluate_dt(emf,emg,h,mp);
while t < T

    if t+dt>T
        dt = T-t;
    end

    % 使用RK方法求解该时间步下的结果，先求q再求u
    % io == 1
    for n = 2:ng+1
        qh(:,n-1) = H(:,:,2)*uh_m(:,n)+L(:,:,2)*uh_m(:,n-1);
    end
    qh_m = [qh(:,ng),qh,qh(:,1)];
    
    for n = 2:ng+1
        result(:,n) = uh_m(:,n) + dt*(-H(:,:,2)*uh_m(:,n)-L(:,:,2)*uh_m(:,n-1)+R(:,:,1)*qh_m(:,n+1)+H(:,:,1)*qh_m(:,n)+f_term(:,n));
    end
    result(:,1) = result(:,ng+1);
    result(:,ng+2) = result(:,2);

    % io = 2
    for n = 2:ng+1
        qh(:,n-1) = H(:,:,2)*result(:,n)+L(:,:,2)*result(:,n-1);
    end
    qh_m = [qh(:,ng),qh,qh(:,1)];
    
    f_term = generate_1D_vector(Gauss_coefficient,ng,mp,h,mid_points,inv_mass,t+dt/2);
    for n = 2:ng+1
        result1(:,n) = 3/4*uh_m(:,n) + 1/4*(result(:,n)+dt*(-H(:,:,2)*result(:,n)-L(:,:,2)*result(:,n-1)+R(:,:,1)*qh_m(:,n+1)+H(:,:,1)*qh_m(:,n)+f_term(:,n)));
    end
    result1(:,1) = result1(:,ng+1);
    result1(:,ng+2) = result1(:,2);

    % io = 3
    for n = 2:ng+1
        qh(:,n-1) = H(:,:,2)*result1(:,n)+L(:,:,2)*result1(:,n-1);
    end
    qh_m = [qh(:,ng),qh,qh(:,1)];
    
    f_term = generate_1D_vector(Gauss_coefficient,ng,mp,h,mid_points,inv_mass,t+dt);
    for n = 2:ng+1
        uh_m(:,n) = 1/3*uh_m(:,n) + 2/3*(result1(:,n)+dt*(-H(:,:,2)*result1(:,n)-L(:,:,2)*result1(:,n-1)+R(:,:,1)*qh_m(:,n+1)+H(:,:,1)*qh_m(:,n)+f_term(:,n)));
    end
    uh_m(:,1) = uh_m(:,ng+1);
    uh_m(:,ng+2) = uh_m(:,2);

    t = t+dt;

end

% % 求 t = T时刻的qh
% for n = 2:ng+1
%     qh(:,n-1) = H(:,:,2)*uh_m(:,n)+L(:,:,2)*uh_m(:,n-1);
% end
% 求 t = T时刻的uh
uh = uh_m(:,2:ng+1);
