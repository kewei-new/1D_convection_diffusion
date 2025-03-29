% main.m
% 使用IPDG方法，求解energy transition模型
clear
clc
% 计算区域
% left:左端点；right: 右端点
% T_end:终止时刻

left = 0; right = 2*pi;
T_end = 1;

%% 输入离散相关参数
% ng:网格数
% mo:有限元空间阶数，mp:有限元多项式次数
ng_base = 30;
mo =3;
mp = mo-1;
nt = 3;
%% Gauss quadrature
% 生成数值积分相关系数，这里使用Gauss-lobatto
% 根据有限元理论精度为mp+1阶，Gauss积分精度至少为mp+1
Gauss_coefficient = generate_Gauss_lobatto(mp);

%% inverse mass matrix
inv_mass = generate_1D_invmass(mp);
error2_Last = 1;
errors_Last = 1;
error_order = zeros(nt,2);
for k = 1:nt
    
    ng = ng_base * 2^(k-1);
    %% 剖分
    [h,mid_points] = generate_1D_uniform_grid(left,right,ng);
    
    %% 画图
    % solve_1D_ET:主要过程
    [u1,phih,phix] = solve_1D_DD(Gauss_coefficient,inv_mass,ng,mp,h,mid_points,T_end);
    
    % 求各个单元中点的值并画图
    % mid_values = plot_1D_ET(mid_points,u1,phih,phix,T_end);
    %% evvaluate error
    error_s = evaluate_error_absolute('exact_fun_phi',T_end,phih,h,mid_points,ng,mp);
    error_2 = evaluate_error_Lnorm('exact_fun_phi',2,T_end,phih,ng,mp,h,mid_points);

    error_order(k,:) = [log(error2_Last/error_2)/log(2),log(errors_Last/error_s)/log(2)];
    error2_Last = error_2;
    errors_Last = error_s;
end





