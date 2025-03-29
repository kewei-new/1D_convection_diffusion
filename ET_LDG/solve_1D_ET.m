function  [u1,u2,phih,phix] =  solve_1D_ET(Gauss_coefficient,inv_mass,ng,mp,h,mid_points,T_end)

%% 输入已知参数
% 为了进行区分，会在一些参数前加x，例如：参数e，写作xe
% 变量：xe,xm,xk,xepsilon,xni

xm = 0.26*0.9109;
xe = 0.1602;
xk = 0.138046e-4;
epsilon = 11.7*8.85418;
% xni = 0.014;
T_0 = 300;

% xnd = nd(X)  调用函数nd可以生成xnd
% xmu0 = mu0(X) 调用函数mu0可以生成xnd
% 都可以实现R^n到R^n的转化

%% L2_projection
% u=[u1,u2];u1=xe*n;u2=n*E/m;
u1 = L2_projection_1D('init_fun1',inv_mass,ng,mp,h,mid_points);
u2 = L2_projection_1D('init_fun2',inv_mass,ng,mp,h,mid_points);
% periodically treatment
uh1 = [u1(:,ng),u1,u1(:,1)];
uh2 = [u2(:,ng),u2,u2(:,1)];
mid_points_p = [mid_points(ng),mid_points,mid_points(1)];

%% LDG constant matrix
%  use LDG mtthods to obtain electric field -/phix
% 由于不是每个方程都是possion方程这么简单，考虑变系数，就把求系数这一块放在里面了
[H,R,L] = generate_1D_LDG_matrix(Gauss_coefficient,mp,inv_mass,h);
[LDG_matrix,f_values,Q_matrix] = assemble_LDG_matrix(inv_mass,ng,mp,mid_points,h,H,R,L);

%% evaluate dt 
% need to evaluate emf(f(u)jacobian matrix's eigenvalue) & emg(g(u)jacobian's eigenvalue)
num = 0;
[phih,phix] = solve_1D_LDG(Gauss_coefficient,LDG_matrix,Q_matrix,f_values,u1,h,mid_points,mp,ng,num);
phixh = [phix(:,ng),phix,phix(:,1)];
[emf,emg]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1,uh2,phixh);
dt = evaluate_dt(emf,emg,h,mp);

%% time evolution
t_m = 0;
uh1_temp = zeros(mp+1,ng+2,2);
uh2_temp = zeros(mp+1,ng+2,2);

% 误差
u1_last = zeros(1,ng);
u2_last = zeros(1,ng);
for k = 1:ng
    u1_last(k) = evaluate_uh_DG(0,uh1(:,k+1),mp);
    u2_last(k) = evaluate_uh_DG(0,uh2(:,k+1),mp);
end
error = 100;

while (t_m<T_end) && (error>1e-6)
    
    num = num+1;
    print_progress(t_m, T_end);

    if t_m + dt>=T_end
        dt = T_end - t_m;
    end
    
    % 使用RK方法求解该时间步下的结果
    % io = 1
    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1,uh2,mp,h,mid_points_p,xm,xe,xk,T_0,epsilon,emf);
    uh1_temp(:,:,1) = uh1 + dt*(f_plus_g_add_h(:,:,1));
    uh2_temp(:,:,1) = uh2 + dt*(f_plus_g_add_h(:,:,2));
    
    u1 = uh1_temp(:,2:ng+1,1);

    % io = 2
    [~,phix] = solve_1D_LDG(Gauss_coefficient,LDG_matrix,Q_matrix,f_values,u1,h,mid_points,mp,ng,num);
    % mid_values = plot_1D_ET(mid_points,uh1_temp(:,2:ng+1,1),uh2_temp(:,2:ng+1,1),phih,phix);
    phixh = [phix(:,ng),phix,phix(:,1)];
    % [emf,~]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1_temp(:,:,1),uh2_temp(:,:,1),phixh);

    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1_temp(:,:,1),uh2_temp(:,:,1),mp,h,mid_points_p,xm,xe,xk,T_0,epsilon,emf);
    uh1_temp(:,:,2) = 3/4*uh1 +1/4*uh1_temp(:,:,1) + 1/4*dt*(f_plus_g_add_h(:,:,1));
    uh2_temp(:,:,2) = 3/4*uh2 +1/4*uh2_temp(:,:,1) + 1/4*dt*(f_plus_g_add_h(:,:,2));

    u1 = uh1_temp(:,2:ng+1,2);

    % io = 3
    [~,phix] = solve_1D_LDG(Gauss_coefficient,LDG_matrix,Q_matrix,f_values,u1,h,mid_points,mp,ng,num);
    phixh = [phix(:,ng),phix,phix(:,1)];
    % [emf,~]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1_temp(:,:,2),uh2_temp(:,:,2),phixh);

    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1_temp(:,:,2),uh2_temp(:,:,2),mp,h,mid_points_p,xm,xe,xk,T_0,epsilon,emf);
    uh1 = 1/3*uh1 +2/3*uh1_temp(:,:,2) + 2/3*dt*(f_plus_g_add_h(:,:,1));
    uh2 = 1/3*uh2 +2/3*uh2_temp(:,:,2) + 2/3*dt*(f_plus_g_add_h(:,:,2));

    t_m = t_m+dt;
    % 
    u1 = uh1(:,2:ng+1);
    [phih,phix] = solve_1D_LDG(Gauss_coefficient,LDG_matrix,Q_matrix,f_values,u1,h,mid_points,mp,ng,num);
    phixh = [phix(:,ng),phix,phix(:,1)];

    % 计算下一个时间步步长
    [emf,emg]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1,uh2,phixh);
    dt = evaluate_dt(emf,emg,h,mp);
    
    % 误差
    error = 0;
    for k = 1:ng
            u1_local = evaluate_uh_DG(0,uh1(:,k+1),mp);
            u2_local = evaluate_uh_DG(0,uh2(:,k+1),mp);
            error = error + abs(u1_local - u1_last(k)) + abs(u2_local - u2_last(k));
            u1_last(k) = u1_local;
            u2_last(k) = u2_local;
    end
    error = error/(2*ng);
    
    % mid_values = plot_1D_ET(mid_points,uh1(:,2:ng+1),uh2(:,2:ng+1),phih,phix);
end

u1 = uh1(:,2:ng+1);
u2 = uh2(:,2:ng+1);

if t_m==T_end
    fprintf('没有收敛到稳态！')
end

%% 能正常运行的话，删掉这里
% mid_values = plot_1D_ET(mid_points,uh1(:,2:ng+1),uh2(:,2:ng+1),phih,phix);
% mid_values = plot_1D_ET(mid_points,uh1_temp(:,2:ng+1,1),uh2_temp(:,2:ng+1,1),phih,phix);


