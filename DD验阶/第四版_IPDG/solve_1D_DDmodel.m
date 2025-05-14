function error = solve_1D_DDmodel(mo,k)
% nt - (En)_x = n_xx + f1
% \phi_xx = n - f2

% n = sin(x)*cos(t)
% \phi = sin(t)*cos(x)
% E = sin(t)sin(x)


%% nargin
if nargin < 1
    k = 2;
    mo= 2;
elseif nargin < 2
    k = 1;
end

%% parameters
% global xk epsilon xe xm T_0 xmu ni vbias

left = 0;
right = 2*pi;
T_end = 0.5;
ng = 10*2^(k-1);
mp = mo-1;
mt = 3;

% FEM basis type
if mp == 1
    basis_type = 101;
else
    basis_type = 102;
end

tic
%% reference data
Gauss_coefficient = generate_1D_gaussian_quadrature_r(mp);
inv_mass = generate_1D_inverse_mass_matrix(mp);

%% segments
[P_partition,T_partition] = generate_1D_grid(left,right,ng);
h = max(diff(P_partition));
mid_points = (P_partition(T_partition(1,:))+P_partition(T_partition(2,:)))/2;

% FEM partition
if basis_type == 101
    Pb = P_partition;
    Tb = T_partition;
    % 周期边界条件
    Tb = [Tb(:,ng),Tb,Tb(:,1)];
else
    [Pb,Tb] = generate_1D_FEM_partition(P_partition,T_partition);
    % 周期边界条件
    Tb = [Tb(:,ng),Tb,Tb(:,1)];
end

% C_penalty:penalty coefficient of IPDG, The larger the better, but too large leads to unstable solutions
% and the maximum value is related to mp:For mp=1, C_penalty=2; for mp=2, C_penalty=10; for mp=3, C_penalty=22
if mp == 1
   C_penalty = 2;
elseif mp == 2
    C_penalty = 8;
elseif mp == 3
    C_penalty = 22;
end
theta = C_penalty/h;

%% initial condition with L2-projection
n_m = L2_projection_1D(inv_mass,ng,mp,h,mid_points);
% plot_n_DG(0,n_m,P_partition,T_partition);
% periodic boundary
nh = n_m;
n_m = [n_m(:,ng),n_m,n_m(:,1)];

%% time evolution
emf=1;
emg=1;
cflc = generate_cfl_condition(mp,11);
cfld = generate_cfl_condition(mp,12);
dt = evaluate_dt(cflc,cfld,emg,emf,P_partition,T_partition);
t_m = 0;

% t1 = toc
while t_m < T_end
    
    if t_m + dt>=T_end
        dt = T_end - t_m;
    end
    
    % 使用RK方法求解该时间步下的结果
    % io = 1
    [phih,~] = solve_poisson_1D_FEM(t_m,mp,ng,P_partition,T_partition,Pb,Tb,Gauss_coefficient,h,mid_points,nh,basis_type);
    res = evaluate_f_g(Gauss_coefficient,inv_mass,Tb,phih,n_m,mp,h,mid_points,emf,theta,t_m);
    result = n_m + dt*res;
    result(:,1) = result(:,ng+1);
    result(:,ng+2) = result(:,2);
    nh = result(:,2:ng+1);
    
    % t2 = toc-t1

    % io = 2
    [phih,~] = solve_poisson_1D_FEM(t_m+dt/2,mp,ng,P_partition,T_partition,Pb,Tb,Gauss_coefficient,h,mid_points,nh,basis_type);
    res = evaluate_f_g(Gauss_coefficient,inv_mass,Tb,phih,result,mp,h,mid_points,emf,theta,t_m+dt/2);
    result1 = 3/4*n_m + 1/4*(result+ dt*res);
    result1(:,1) = result1(:,ng+1);
    result1(:,ng+2) = result1(:,2);
    nh = result1(:,2:ng+1);

    % io = 3
    [phih,~] = solve_poisson_1D_FEM(t_m+dt,mp,ng,P_partition,T_partition,Pb,Tb,Gauss_coefficient,h,mid_points,nh,basis_type);
    res = evaluate_f_g(Gauss_coefficient,inv_mass,Tb,phih,result1,mp,h,mid_points,emf,theta,t_m+dt);
    n_m = 1/3*n_m + 2/3*(result1 + dt*res);
    n_m(:,1) = n_m(:,ng+1);
    n_m(:,ng+2) = n_m(:,2);
    nh = n_m(:,2:ng+1);

    t_m = t_m+dt;

end
% t3 = toc;


% evaluate errors

error1 = evaluate_error_Lnorm('exact_fun',1,T_end,nh,ng,mp,h,mid_points);
error2 = evaluate_error_Lnorm('exact_fun',2,T_end,nh,ng,mp,h,mid_points);
errorf = evaluate_error_absolute('exact_fun',T_end,nh,P_partition,T_partition,ng,mp);

error=[error1,error2,errorf];

subplot(2,2,1)
plot_phi_FEM(T_end,phih,P_partition,T_partition,Tb,basis_type,0);
subplot(2,2,2)
plot_phi_FEM(T_end,phih,P_partition,T_partition,Tb,basis_type,1);
subplot(2,2,[3,4])
plot_n_DG(T_end,nh,P_partition,T_partition);

