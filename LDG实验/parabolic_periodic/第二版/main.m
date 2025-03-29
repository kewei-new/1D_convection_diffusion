% pure_LDG main.m
% 计算误差和误差阶
clc
clear

left = 0;
right = 2*pi;
T_end = 0.5;

ng_base = 40;
mp = 2;

Gauss_coefficient = generate_1D_gaussian_quadrature_r(mp);
inv_mass = generate_1D_inverse_mass_matrix(mp); 

uh_error_order = zeros(2,2);
qh_error_order = zeros(2,2);
uh_error = zeros(3,2);
qh_error = zeros(3,2);

for k = 1:3

    ng = ng_base * 2^(k-1);

    [h,mid_points] = generate_uniform_grid(left,right,ng);


    [uh,qh]=solve_1D_LDG(Gauss_coefficient,inv_mass,T_end,h,mid_points,mp,ng);
    uh_L2_error = evaluate_error_Lnorm('exact_fun',2,T_end,uh,ng,mp,h,mid_points);
    uh_Linf_error = evaluate_error_absolute('exact_fun',T_end,uh,h,mid_points,ng,mp);
    qh_L2_error = evaluate_error_Lnorm('exact_fun_x',2,T_end,qh,ng,mp,h,mid_points);
    qh_Linf_error = evaluate_error_absolute('exact_fun_x',T_end,qh,h,mid_points,ng,mp);
    
    uh_error(k,1) = uh_L2_error;
    uh_error(k,2) = uh_Linf_error;
    qh_error(k,1) = qh_L2_error;
    qh_error(k,2) = qh_Linf_error;

    % plot_DG('exact_fun_x',T_end,qh,ng,mp,mid_points)
end

plot_DG('exact_fun_x',T_end,qh,ng,mp,mid_points)

for i = 2:3
    for j = 1:2
        uh_error_order(i,j) = log(uh_error(i-1,j)/uh_error(i,j))/log(2);
        qh_error_order(i,j) = log(qh_error(i-1,j)/qh_error(i,j))/log(2);
    end
end

