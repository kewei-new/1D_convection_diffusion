function error=solver_1D_CD_IPDG(k)
% 用来求解一维convection-diffusion，ut+ux-uxx=0,u(x,0)=sinx,周期边界条件
% uh:多项式近似解的系数向量
% k:用来配合计算误差阶(check_error_1D.m)

if nargin < 1
    k = 1;
end

%% 参数
% left,right:boundary
% emg:Convection velocity, the coefficient of the convection term
% emf:diffusion velocity,the coefficient of the diffusion term
% T_end:the last moment in time
% ng:num grids
% mp:hightest order of polynomial
% mt:Choose between RK3 and RK4 (mt = 3 for RK3, mt = 4 for RK4)

format long;
left = 0;
right = 2*pi;
emg = 1;
emf = 1;
T_end = 2*pi;
ng = 80*2^(k-1);
mp = 1;
mt = 3;

% trial/test_basis_type: trial/test function's basis type
% ==101:orthogonal legendre basis with highest order's coefficient is 1
test_basis_type = 101;
trial_basis_type = 101;

%% Gauss-lobatto and mass-matrix(legendre & [-0.5,0.5])
% GP_r:Gauss ponits in [-0.5,0.5]
% GW_r:Gauss quadrature weights in [-0.5,0.5]
% Gauss_coefficient(:,1)==GP_r;Gauss_coefficient(:,2)==GW_r;
% a:inverse of mass matrix evaluate the basis in [-0.5,0.5]
% b:stiffness matrix evaluate the basis in [-0.5,0.5]
% c:evaluate inner product (u,ux) in [-0.5,0.5]
% d:evaluate basis fun mutiply on the boundary with derivative 0,and d is a
% 4-element set which shape is (mp+1,mp+1,4) 
% e:evaluate basis fun mutiply on the boundary with derivative 1,and e is a
% 4-element set which shape is (mp+1,mp+1,4) 

Gauss_coefficient = generate_1D_gaussian_quadrature_r(mp);

% after the first evaluation, can save the result in txt, this part is
% independent from other work,such like use Gauss quadrature formula
a = generate_1D_matrix_r(Gauss_coefficient,mp,mp,trial_basis_type,0,test_basis_type,0);
b = generate_1D_matrix_r(Gauss_coefficient,mp,mp,trial_basis_type,1,test_basis_type,1);
c = generate_1D_matrix_r(Gauss_coefficient,mp,mp,trial_basis_type,0,test_basis_type,1);
d = generate_1D_boundary_multiple(mp,trial_basis_type,0,test_basis_type,0);
e = generate_1D_boundary_multiple(mp,trial_basis_type,1,test_basis_type,0);

%% Segment & related coefficients
% P_partition:剖分节点坐标,size is (1,ng+1)
% T_partition:剖分索引,size is (2,ng)
% T_partition(1,*):left point of the *th segment
% T_partition(2,*):right point of the *th segment
% h:max_diam(S_h), where S_h is the space of segment
% theta: interior penalty coefficient in IPDG method

[P_partition,T_partition] = generate_grid_1D(left,right,ng);
h = max(diff(P_partition));
theta = 2/h;

%% initial condition with L2-projection
% uh0：coefficient of the discrete initial value conditions

uh0 = L2_projection_1D('initial_condition',a,ng,mp,P_partition,T_partition);
% plot_solution('initial_condition',0,uh0,ng,mp,P_partition,T_partition)
% error1 = evaluate_error_Lnorm('initial_condition',1,0,uh0,ng,mp,P_partition,T_partition)

%% evaluate the left terms & right terms

% left terms:generate coefficient of the unknown variables
% not consider interior boundaries,would like to get such scheme: uht + B/A * uh - C/A * uh
% M1=B/A; M2=-C/A;
% consider interior boundaries, it could influent right segment,
% left segment and here(in) segemnt; so we use R,L,H to express these situations.
% R1,L1,H1:generate by convection;R2,L2,H2:generate by diffusion
inv_a = inv(a);

M1 = generate_1D_convection_diffusion_matrix_local(inv_a,b,c,ng,mp,mp,P_partition,T_partition,12);
M2 = -generate_1D_convection_diffusion_matrix_local(inv_a,b,c,ng,mp,mp,P_partition,T_partition,11);
[H1,R1,L1] = generate_1D_convection_diffusion_flux_local(inv_a,d,e,theta,ng,mp,P_partition,T_partition,11);
[H2,R2,L2] = generate_1D_convection_diffusion_flux_local(inv_a,d,e,theta,ng,mp,P_partition,T_partition,12);

M = M1+M2;
H = H1+H2;
R = R1+R2;
L = L1+L2;

% general matrix
% uht + A*uh = b
A = assemble_general_1D_matrix(M,H,R,L,ng,mp);

% right terms:involving all known value
% evaluate b = (f,v)

% 这里有问题，需要考虑质量矩阵，但是这个问题右端为0不需要改
rs = generate_1D_vector("f_fun",Gauss_coefficient,ng,mp,P_partition,T_partition);

%% boundary condition treatment
% boundary_nodes: involves boundary types, nodes in boundary
% 这里是周期边界条件，为了一般化程序所以加了这样一步，周期边界最方便是在assemble那一步直接把1，n这两个单元的左右补全
boundary_nodes = [1;1;1];
[A,rs] = boundary_treatment_1D(boundary_nodes,A,rs,ng,mp);

%% time evoulation
% choose CFL condition
cflc = generate_cfl_condition(mp,11);
cfld = generate_cfl_condition(mp,12);
dt = evaluate_dt(cflc,cfld,emg,emf,P_partition,T_partition);

% Runge-Kutta scheme for advancing in time
% solute uht = f(uh) = -A*uh + b
% 
rhs = @(uh)(rs-A*uh);
[uh_list,t_list] = time_evoulation_with_RK(T_end,uh0,rhs,dt,mt);

%% evaluate error(involving L1,L2,L\infinite)

error1 = evaluate_error_Lnorm('exact_fun',1,t_list(end),uh_list(:,end),ng,mp,P_partition,T_partition);
error2 = evaluate_error_Lnorm('exact_fun',2,t_list(end),uh_list(:,end),ng,mp,P_partition,T_partition);
errorf = evaluate_error_absolute('exact_fun',t_list(end),uh_list(:,end),P_partition,T_partition,ng,mp);

error=[error1,error2,errorf];
plot_solution('exact_fun',t_list(end),uh_list(:,end),ng,mp,P_partition,T_partition)



