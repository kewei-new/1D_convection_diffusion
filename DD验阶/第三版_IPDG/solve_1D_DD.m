function  [u1,phih,phix] =  solve_1D_DD(Gauss_coefficient,inv_mass,ng,mp,h,mid_points,T_end)
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
    C_div_h = 3;
elseif mp==2
    C_div_h = 10;
elseif mp==3
    C_div_h = 22;
end


%% L2_projection
u1 = L2_projection_1D('init_fun1',inv_mass,ng,mp,h,mid_points);

%% evaluate dt 
% need to evaluate emf(f(u)jacobian matrix's eigenvalue) & emg(g(u)jacobian's eigenvalue)
[phih,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,h,mid_points,mp,ng,u1,0);
emf = 1;
emg = 1;
dt = evaluate_dt(emf,emg,h,mp);

%% time evolution
t_m = 0;
while  t_m<T_end
    
    % print_progress(t_m, T_end);

    if t_m + dt>=T_end
        dt = T_end - t_m;
    end
    
    % 使用RK方法求解该时间步下的结果
    % io = 1
    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phix,u1,mp,h,mid_points,emf,C_div_h,t_m);
    uh1_temp1 = u1 + dt*(f_plus_g_add_h);

    % io = 2
    [~,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,h,mid_points,mp,ng,uh1_temp1,t_m+dt/2);
    % [emf,~]=evaluate_fluid_velocity(ng,mp,mid_points,h,uh1_temp(:,:,1),uh2_temp(:,:,1),phixh);

    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phix,uh1_temp1,mp,h,mid_points,emf,C_div_h,t_m+dt/2);
    uh1_temp2 = 3/4*u1 +1/4*uh1_temp1 + 1/4*dt*(f_plus_g_add_h);

    % io = 3
    [~,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,h,mid_points,mp,ng,uh1_temp2,t_m+dt);
    % [emf,~]=evaluate_fluid_velocity(ng,mp,mid_points,h,uh1_temp(:,:,2),uh2_temp(:,:,2),phixh);

    f_plus_g_add_h = evaluate_f_g(Gauss_coefficient,inv_mass,phix,uh1_temp2,mp,h,mid_points,emf,C_div_h,t_m+dt);
    u1 = 1/3*u1 +2/3*uh1_temp2 + 2/3*dt*(f_plus_g_add_h);
    % uh1(:,1) = uh1(:,ng+1);
    % uh1(:,ng+2) = uh1(:,2);
    
    t_m = t_m+dt;
    % 
    [phih,phix] = solve_1D_LDG(Gauss_coefficient,inv_mass,h,mid_points,mp,ng,u1,t_m);


    

    % mid_values = plot_1D_ET(mid_points,u1,phih,phix,t_m);
end


% if t_m==T_end
%     fprintf('没有收敛到稳态！')
% end

%% 能正常运行的话，删掉这里
% mid_values = plot_1D_ET(mid_points,uh1(:,2:ng+1),uh2(:,2:ng+1),phih,phix);
% mid_values = plot_1D_ET(mid_points,uh1_temp(:,2:ng+1,1),uh2_temp(:,2:ng+1,1),phih,phix);


