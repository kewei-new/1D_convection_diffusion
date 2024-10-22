function [uh,t] = time_evoulation_with_RK_variable_coefficient(Gauss_coefficient,M,inv_a,d,T_end,uh0,dt,mp,mt,P,T,trial_basis_type,test_basis_type)
% time evolution with RK method with variable coefficient
% T_end:the end moment of time evolution
% dt:time step of every itertion

% maxiter:max iteration number
% nt:current iteration number
% t_e:current moment

maxiter = 100000000;
% 当前时刻
t_e = 0;

% 用来存储，无实际意义
t = t_e;
uh = [];

xmu = 0.75;
ng = length(T);

for nt=1:maxiter

    % 达到末时刻就停止
    if t_e == T_end
        break;
    end
    
    % 控制不能大于末时刻
    if t_e + dt>=T_end
        dt = T_end - t_e;
        t_e = T_end;
    else
        t_e = t_e + dt;
    end
    
    % 使用RK方法求解该时间步下的结果
    result = zeros(length(uh0),mt-1);

    uhtemp = uh0;

    if mt==3

    % io = 1
    [phi,Tb] = solve_poisson_1D_FEM(mp,P,T,Gauss_coefficient,uh0);
    M2 = generate_1D_variable_coefficient_matrix_local('evaluate_local_FEM',inv_a,Gauss_coefficient,P,T,trial_basis_type,0,test_basis_type,1,Tb,phi);
    [L,H,R] = generate_1D_variable_coefficient_flux_local(inv_a,d,phi,mp,ng,P,T,Tb);
    A = xmu*assemble_1D_matrix_with_variable(ng,mp,M2,L,H,R);
    b = zeros(ng*(mp+1),1);
    % b = generate_1D_vector('f1_fun',t,Gauss_coefficient,ng,mp,P,T);
    % boundary_nodes = [1;1;1];
    % [A,b] = boundary_treatment_1D(boundary_nodes,A+M1,b,ng,mp);
    result(:,1) = uh0 + dt*rhs_variable_coefficient(uh0,A,M,b);

    % io = 2
    [phi,Tb] = solve_poisson_1D_FEM(mp,P,T,Gauss_coefficient,result(:,1));
    M2 = generate_1D_variable_coefficient_matrix_local('evaluate_local_FEM',inv_a,Gauss_coefficient,P,T,trial_basis_type,0,test_basis_type,1,Tb,phi);
    [L,H,R] = generate_1D_variable_coefficient_flux_local(inv_a,d,phi,mp,ng,P,T,Tb);
    A = xmu*assemble_1D_matrix_with_variable(ng,mp,M2,L,H,R);
    b = zeros(ng*(mp+1),1);
    % b = generate_1D_vector('f1_fun',t+dt/2,Gauss_coefficient,ng,mp,P,T);
    % boundary_nodes = [1;1;1];
    % [A,b] = boundary_treatment_1D(boundary_nodes,A+M1,b,ng,mp);
    result(:,2) = 3/4*uh0 + 1/4*result(:,1) + 1/4*dt*rhs_variable_coefficient(result(:,1),A,M,b);

    % io = 3
    [phi,Tb] = solve_poisson_1D_FEM(mp,P,T,Gauss_coefficient,result(:,2));
    M2 = generate_1D_variable_coefficient_matrix_local('evaluate_local_FEM',inv_a,Gauss_coefficient,P,T,trial_basis_type,0,test_basis_type,1,Tb,phi);
    [L,H,R] = generate_1D_variable_coefficient_flux_local(inv_a,d,phi,mp,ng,P,T,Tb);
    A = xmu*assemble_1D_matrix_with_variable(ng,mp,M2,L,H,R);
    b = zeros(ng*(mp+1),1);
    % b = generate_1D_vector('f1_fun',t+dt,Gauss_coefficient,ng,mp,P,T);
    % boundary_nodes = [1;1;1];
    % [A,b] = boundary_treatment_1D(boundary_nodes,A+M1,b,ng,mp);
    uh0 = 1/3*uh0 + 2/3*result(:,2) + 2/3*dt*rhs_variable_coefficient(result(:,2),A,M,b);
    
    end

    % 每100次步长记录一下
    if mod(nt-1,100) == 0
        result_list = [uh,uhtemp];
        uh = result_list;
        T_list = [t,t_e];
        t = T_list;
    end

end

end
