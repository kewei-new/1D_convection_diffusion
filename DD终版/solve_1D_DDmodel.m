function solve_1D_DDmodel(mo,k)

%% nargin
if nargin < 1
    k = 3;
    mo= 3;
elseif nargin < 2
    k = 1;
end

%% parameters
% global xk epsilon xe xm T_0 xmu ni vbias
xk=0.138046e-4; epsilon=11.7*8.85418;xe=0.1602;xm=0.26*0.9109;T_0=300;xmu=0.75;
ni=0.014;vbias=1.5;

left = 0;
right = 0.6;
T_end = 0.6;
ng = 80*2^(k-1);
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
inv_a = generate_1D_inverse_mass_matrix(mp);

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
    C_penalty = 10;
elseif mp == 3
    C_penalty = 22;
end
theta = C_penalty/h;

%% initial condition with L2-projection
n_m = L2_projection_1D(inv_a,ng,mp,h,mid_points);
% periodic boundary
nh = n_m;
n_m = [n_m(:,ng),n_m,n_m(:,1)];

%% evaluate constant matrix
cha = 1;

M = assemble_1D_matrix_from_integral(Gauss_coefficient,inv_a,h,mp);
[H1,R1,L1] = generate_1D_center_flux_local(inv_a,h,mp,cha);
[H2,R2,L2] = generate_1D_jump_flux_local(inv_a,h,mp);
% 外面加负号，是因为RK要把这些项移到右侧
H = -(M+(-H1+theta*H2))*xmu*xk*T_0/xe;
L = -(-L1+theta*L2)*xmu*xk*T_0/xe;
R = -(-R1+theta*R2)*xmu*xk*T_0/xe;

%% time evolution
emf=xmu;
emg=0.8*xmu*xk*T_0/xe;
cflc = generate_cfl_condition(mp,11);
cfld = generate_cfl_condition(mp,12);
dt = evaluate_dt(cflc,cfld,emg,emf,P_partition,T_partition);
t_m = 0;
result = zeros(size(n_m,1),size(n_m,2));
result1 = zeros(size(n_m,1),size(n_m,2));

t1 = toc


while t_m < T_end
    
    if t_m + dt>=T_end
        dt = T_end - t_m;
    end
    
    % 使用RK方法求解该时间步下的结果
    % io = 1
    [phi,E] = solve_poisson_1D_FEM(mp,ng,P_partition,T_partition,Pb,Tb,Gauss_coefficient,h,mid_points,nh,basis_type);
    for n = 2:ng+1
        % M2 = assemble_1D_coefficient_matrix_from_integral(Gauss_coefficient,inv_a,phi,Tb,h,mp,n-1,xmu);
        [H3,R3,L3] = generate_1D_coefficient_center_flux_local(Gauss_coefficient,inv_a,phi,E,Tb,mp,h,n,xmu,basis_type);
        result(:,n) = n_m(:,n) + dt*((L+L3)*n_m(:,n-1)+(H+H3)*n_m(:,n)+(R+R3)*n_m(:,n+1));
    end
    result(:,1) = result(:,ng+1);
    result(:,ng+2) = result(:,2);
    nh = result(:,2:ng+1);

    t2 = toc-t1;

    % io = 2
    [phi,E] = solve_poisson_1D_FEM(mp,ng,P_partition,T_partition,Pb,Tb,Gauss_coefficient,h,mid_points,nh,basis_type);
    for n = 2:ng+1
        % M2 = assemble_1D_coefficient_matrix_from_integral(Gauss_coefficient,inv_a,phi,Tb,h,mp,n-1,xmu);
        [H3,R3,L3] = generate_1D_coefficient_center_flux_local(Gauss_coefficient,inv_a,phi,E,Tb,mp,h,n,xmu,basis_type);
        result1(:,n) = 3/4*n_m(:,n) + 1/4*(result(:,n) + dt*((L+L3)*result(:,n-1)+(H+H3)*result(:,n)+(R+R3)*result(:,n+1)));
    end
    result1(:,1) = result1(:,ng+1);
    result1(:,ng+2) = result1(:,2);
    nh = result1(:,2:ng+1);

    % io = 3
    [phi,E] = solve_poisson_1D_FEM(mp,ng,P_partition,T_partition,Pb,Tb,Gauss_coefficient,h,mid_points,nh,basis_type);
    for n = 2:ng+1
        % M2 = assemble_1D_coefficient_matrix_from_integral(Gauss_coefficient,inv_a,phi,Tb,h,mp,n-1,xmu);
        [H3,R3,L3] = generate_1D_coefficient_center_flux_local(Gauss_coefficient,inv_a,phi,E,Tb,mp,h,n,xmu,basis_type);
        n_m(:,n) = 1/3*n_m(:,n) + 2/3*(result1(:,n) + dt*((L+L3)*result1(:,n-1)+(H+H3)*result1(:,n)+(R+R3)*result1(:,n+1)));
    end
    n_m(:,1) = n_m(:,ng+1);
    n_m(:,ng+2) = n_m(:,2);
    nh = n_m(:,2:ng+1);

    t_m = t_m+dt;

end

t3 = toc
subplot(2,2,1)
plot_phi_FEM(phi,P_partition,T_partition,Tb,basis_type,0);
subplot(2,2,2)
plot_phi_FEM(phi,P_partition,T_partition,Tb,basis_type,1);
subplot(2,2,[3,4])
plot_n_DG(n_m,P_partition,T_partition)
