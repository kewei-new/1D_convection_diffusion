% main.m
% 使用IPDG方法，求解energy transition模型
clear
clc
% 计算区域
% left:左端点；right: 右端点
% T_end:终止时刻

left = 0; right = 0.6;
T_end = 20;

%% 输入离散相关参数
% ng:网格数
% mo:有限元空间阶数，mp:有限元多项式次数
ng = 100;
mo =2;
mp = mo-1;

%% Gauss quadrature
% 生成数值积分相关系数，这里使用Gauss-lobatto
% 根据有限元理论精度为mp+1阶，Gauss积分精度至少为mp+1
Gauss_coefficient = generate_Gauss_lobatto(mp);

%% 剖分
[h,mid_points] = generate_1D_uniform_grid(left,right,ng);
% [h,mid_points] = generate_1D_uneven_grid(left,right,ng);

%% inverse mass matrix
inv_mass = generate_1D_invmass(mp);

%% 画图
% solve_1D_ET:主要过程
[u1,u2,phih,phix] = solve_1D_ET(Gauss_coefficient,inv_mass,ng,mp,h,mid_points,T_end);

% 求各个单元中点的值并画图
mid_values = plot_1D_ET(mid_points,u1,u2,phih,phix);







