function  [u1,phih,phix] =  solve_1D_ET(Gauss_coefficient,inv_mass,ng,mp,h,mid_points,T_end)
% nt - (En)_x = n_xx + f1
% \phi_xx = n - f2

% n = sin(x)*cos(t)
% \phi = sin(t)*cos(x)
% E = sin(t)sin(x)

% phixx = -sin(t)cos(x)
% f2 = sin(t)cos(x)+sin(x)cos(t)
% f1 = -sin(x)sin(t) - sin(t)cos(x)sin(x)cos(t) - cos(x)cos(t)sin(t)sin(x)
% + sin(x)cos(t) = 

if mp==1
    C_div_h = 5/h;
elseif mp==2
    C_div_h = 10/h;
elseif mp==3
    C_div_h = 22/h;
end


%% L2_projection
% u=[u1,u2];u1=xe*n;u2=n*E/m;
u1 = L2_projection_1D('init_fun1',inv_mass,ng,mp,h,mid_points);
% periodically treatment
uh1 = [u1(:,ng),u1,u1(:,1)];
mid_points_p = [mid_points(ng),mid_points,mid_points(1)];

%% LDG constant matrix
%  use LDG mtthods to obtain electric field -/phix
% 由于不是每个方程都是possion方程这么简单，考虑变系数，就把求系数这一块放在里面了
[H,R,L] = generate_1D_LDG_matrix(Gauss_coefficient,mp,inv_mass,h);
[M1,M2,Mn] = assemble_LDG_matrix(inv_mass,ng,mp,mid_points,h,H,R,L,0);

%% evaluate dt 
% need to evaluate emf(f(u)jacobian matrix's eigenvalue) & emg(g(u)jacobian's eigenvalue)
[phih,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,M1,M2,Mn,u1,h,mid_points,mp,ng,0);
phixh = [phix(:,ng),phix,phix(:,1)];
emf = 1;
emg = 1;
% [emf,emg]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1,uh2,phixh);
dt = evaluate_dt(emf,emg,h,mp);

%% time evolution
t_m = 0;
while  t_m<T_end
    
    print_progress(t_m, T_end);

    if t_m + dt>=T_end
        dt = T_end - t_m;
    end
    
    % 使用RK方法求解该时间步下的结果
    % io = 1
    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1,mp,h,mid_points_p,emf,C_div_h,t_m);
    uh1_temp1 = uh1 + dt*(f_plus_g_add_h);
    % uh1_temp1(:,1) = uh1_temp1(:,ng+1);
    % uh1_temp1(:,ng+2) = uh1_temp1(:,2);
    u1 = uh1_temp1(:,2:ng+1);

    % io = 2
    [~,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,M1,M2,Mn,u1,h,mid_points,mp,ng,t_m+dt/2);
    % mid_values = plot_1D_ET(mid_points,uh1_temp(:,2:ng+1,1),uh2_temp(:,2:ng+1,1),phih,phix);
    phixh = [phix(:,ng),phix,phix(:,1)];
    % [emf,~]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1_temp(:,:,1),uh2_temp(:,:,1),phixh);

    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1_temp1,mp,h,mid_points_p,emf,C_div_h,t_m+dt/2);
    uh1_temp2 = 3/4*uh1 +1/4*uh1_temp1 + 1/4*dt*(f_plus_g_add_h);
    % uh1_temp2(:,1) = uh1_temp2(:,ng+1);
    % uh1_temp2(:,ng+2) = uh1_temp2(:,2);
    u1 = uh1_temp2(:,2:ng+1);

    % io = 3
    [~,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,M1,M2,Mn,u1,h,mid_points,mp,ng,t_m+dt);
    phixh = [phix(:,ng),phix,phix(:,1)];
    % [emf,~]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1_temp(:,:,2),uh2_temp(:,:,2),phixh);

    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phixh,uh1_temp2,mp,h,mid_points_p,emf,C_div_h,t_m+dt);
    uh1 = 1/3*uh1 +2/3*uh1_temp2 + 2/3*dt*(f_plus_g_add_h);
    % uh1(:,1) = uh1(:,ng+1);
    % uh1(:,ng+2) = uh1(:,2);
    
    t_m = t_m+dt;
    % 
    u1 = uh1(:,2:ng+1);
    [phih,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,M1,M2,Mn,u1,h,mid_points,mp,ng,t_m);
    phixh = [phix(:,ng),phix,phix(:,1)];

    % % 计算下一个时间步步长
    % [emf,emg]=evaluate_fluid_velocity(ng,mp,mid_points_p,h,uh1,uh2,phixh);
    % dt = evaluate_dt(emf,emg,h,mp);
    

    % mid_values = plot_1D_ET(mid_points,uh1(:,2:ng+1),phih,phix,t_m);
end

u1 = uh1(:,2:ng+1);

% if t_m==T_end
%     fprintf('没有收敛到稳态！')
% end

%% 能正常运行的话，删掉这里
% mid_values = plot_1D_ET(mid_points,uh1(:,2:ng+1),uh2(:,2:ng+1),phih,phix);
% mid_values = plot_1D_ET(mid_points,uh1_temp(:,2:ng+1,1),uh2_temp(:,2:ng+1,1),phih,phix);


