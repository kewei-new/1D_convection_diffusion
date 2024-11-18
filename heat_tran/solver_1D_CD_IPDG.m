function error=solver_1D_CD_IPDG(k,mp)
% 用来求解一维convection-diffusion，ut-uxx=f(x,t),u(x,0)=sinx,周期边界条件
% f(x,t) = sin(x-t)-cos(x-t);
% u(x,t) = sin(x-t);
% k:用来配合计算误差阶(check_error_1D.m)

% ut + 2*ux = f(x,t)
% f(x,t) = cos(x-t)
% ut = -cos(x-t),ux = cos(x-t)

% ut + ux = 0

if nargin < 1
    k = 2;
    mp = 3;
end

%% 参数
% left,right:boundary
% emg:Convection velocity, the coefficient of the convection term
% emf:diffusion velocity,the coefficient of the diffusion term
% T_end:the last moment in time
% ng:num grids
% mp:hightest order of polynomial
% mt:Choose between RK3 and RK4 (mt = 3 for RK3, mt = 4 for RK4)

left = 0;
right = 2*pi;
emg = 10;
emf = 0;
T_end = 1;
ng = 20*2^(k-1);
mt = 3;

%% Gauss-lobatto and mass-matrix(legendre & [-0.5,0.5])
% GP_r:Gauss ponits in [-0.5,0.5]
% GW_r:Gauss quadrature weights in [-0.5,0.5]
% Gauss_coefficient(:,1)==GP_r;Gauss_coefficient(:,2)==GW_r;
% a:inverse of mass matrix evaluate the basis in [-0.5,0.5]
% inv_a:inverse of mass matrix
Gauss_coefficient = generate_1D_gaussian_quadrature_r(mp);

a = generate_1D_matrix_r(Gauss_coefficient,mp,0,0);

% 这里取正交基，除对角线元素均为0
inv_a = inv(a);

%% Segment & related coefficients
% P_partition:剖分节点坐标,size is (1,ng+1)
% T_partition:剖分索引,size is (2,ng)
% T_partition(1,*):left point of the *th segment
% T_partition(2,*):right point of the *th segment
% h:max_diam(S_h), where S_h is the space of segment
% theta: interior penalty coefficient in IPDG method

[P,T] = generate_1D_grid(left,right,ng);
h = max(diff(P));
mid_points = (P(:,T(1,:))+P(:,T(2,:)))/2;
if mp == 1
   C_penalty = 2;
elseif mp == 2
    C_penalty = 8;
elseif mp == 3
    C_penalty = 10;
end
theta =  C_penalty/h;

%% initial condition with L2-projection
% uh0：coefficient of the discrete initial value conditions

uh0 = L2_projection_1D('initial_condition',Gauss_coefficient,inv_a,mid_points,mp,ng,h,0);
% plot_solution('initial_condition',0,uh0,ng,mp,P,T)
% error1 = evaluate_error_Lnorm('initial_condition',1,0,uh0,ng,mp,P,T)

%% evaluate the left terms & right terms
% M = generate_1D_matrix_local(Gauss_coefficient,inv_a,mp,h);
% [H1,R1,L1] = generate_1D_flux(inv_a,mp,h);
% [H2,R2,L2] = generate_1D_symmetrical_flux(inv_a,mp,h);
% [H3,R3,L3] = generate_1D_penalty_flux(inv_a,mp,h);
% 
% H = H1-H2-theta*H3-M;
% R = R1-R2-theta*R3;
% L = L1-L2-theta*L3;

temp = generate_1D_matrix_r(Gauss_coefficient,mp,0,1);
M = 2*inv_a*temp/h;
H = zeros(mp+1,mp+1);
R = zeros(mp+1,mp+1);
L = zeros(mp+1,mp+1);
for m1 = 0:mp
    % test

    for m2 = 0:mp
        % trial

        H(m1+1,m2+1) = 2*reference_basis(0.5,m2,0)*reference_basis(0.5,m1,0);
        L(m1+1,m2+1) = -2*reference_basis(0.5,m2,0)*reference_basis(-0.5,m1,0);
    end

end
H = -inv_a*H*h^(-1)+M;
R = -inv_a*R*h^(-1);
L = -inv_a*L*h^(-1);


%% time evoulation
% choose CFL condition
cflc = generate_cfl_condition(mp,11);
cfld = generate_cfl_condition(mp,12);
dt = evaluate_dt(cflc,cfld,emg,emf,h);

% Runge-Kutta scheme for advancing in time
% solute uht = f(uh) = -A*uh + b
t_e = 0;
tic
b = generate_1D_vector(t_e,Gauss_coefficient,inv_a,mid_points,mp,ng,h);
while t_e < T_end

    if t_e + dt > T_end

        dt = T_end - t_e;

    end
    
    uh = zeros(mp+1,ng);
    uh1 = zeros(mp+1,ng);

    % io = 1
    for n = 1:ng

        if n == 1      
            u = uh0(:,n) + dt*(H*uh0(:,n) + L*uh0(:,ng) + R*uh0(:,n+1) + b(:,n));       
        elseif n == ng
            u = uh0(:,n) + dt*(H*uh0(:,n) + L*uh0(:,n-1) + R*uh0(:,1) + b(:,n));  
        else
            u = uh0(:,n) + dt*(H*uh0(:,n) + L*uh0(:,n-1) + R*uh0(:,n+1) + b(:,n));
        end
        uh(:,n) = u;
    end

    % io = 2
    b = generate_1D_vector(t_e+dt/2,Gauss_coefficient,inv_a,mid_points,mp,ng,h);
    for n = 1:ng
        if n == 1            
            u = 3/4*uh0(:,n) + 1/4*uh(:,n) + 1/4*dt*(H*uh(:,n) + L*uh(:,ng) + R*uh(:,n+1) + b(:,n));      
        elseif n == ng
            u = 3/4*uh0(:,n) + 1/4*uh(:,n) + 1/4*dt*(H*uh(:,n) + L*uh(:,n-1) + R*uh(:,1) + b(:,n)); 
        else
            u = 3/4*uh0(:,n) + 1/4*uh(:,n) + 1/4*dt*(H*uh(:,n) + L*uh(:,n-1) + R*uh(:,n+1) + b(:,n));
        end
        uh1(:,n) = u;
    end

    % io = 3
    b = generate_1D_vector(t_e+dt,Gauss_coefficient,inv_a,mid_points,mp,ng,h);
    for n = 1:ng
        if n == 1
            uh0(:,n) = 1/3*uh0(:,n) + 2/3*uh1(:,n) + 2/3*dt*(H*uh1(:,n) + L*uh1(:,ng) + R*uh1(:,n+1) + b(:,n));     
        elseif n == ng
            uh0(:,n) = 1/3*uh0(:,n) + 2/3*uh1(:,n) + 2/3*dt*(H*uh1(:,n) + L*uh1(:,n-1) + R*uh1(:,1) + b(:,n));
        else
            uh0(:,n) = 1/3*uh0(:,n) + 2/3*uh1(:,n) + 2/3*dt*(H*uh1(:,n) + L*uh1(:,n-1) + R*uh1(:,n+1) + b(:,n));
        end

    end

    t_e = t_e + dt;
    
end
tl = toc


%% evaluate error(involving L1,L2,L\infinite)

error1 = evaluate_error_Lnorm('exact_fun2',1,T_end,uh0,ng,mp,h,mid_points);
error2 = evaluate_error_Lnorm('exact_fun2',2,T_end,uh0,ng,mp,h,mid_points);
errorf = evaluate_error_absolute('exact_fun2',T_end,uh0,P,T,ng,mp);

error=[error1,error2,errorf];
plot_solution('exact_fun2',T_end,uh0,ng,mp,P,T)



