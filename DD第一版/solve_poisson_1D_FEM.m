function [phi,Tb] = solve_poisson_1D_FEM(mp,P,T,Gauss_coefficient,uh)
% solve \phi_{xx} = e/{\epsilon}*(n - nd)
% we know n(x,0), and \phi satisfy the dirichlet boundary condition,so we
% could get \phi0, and after in the io(1,2,3,4)th step in Rk, we need this
% proceding to get the right hand side f(x+h/2,t+h/2);

%% nargin
if nargin < 1
    ng = 5;
    left = 0;
    right = 0.6;
    P = linspace(left,right,ng+1);
    T = [linspace(1,ng,ng);linspace(2,ng+1,ng)];
    mp= 1;
    Gauss_coefficient = generate_1D_gaussian_quadrature_r(mp);
    a = generate_1D_matrix_r(Gauss_coefficient,mp,mp,101,0,101,0);
    uh = L2_projection_1D('initial_condition',a,ng,mp,P,T);
elseif nargin < 2
    P = linspace(0,0.6,11);
    T = [linspace(1,10,10);linspace(2,11,10)];
    Gauss_coefficient = generate_1D_gaussian_quadrature_r(mp);
    a = generate_1D_matrix_r(Gauss_coefficient,mp,mp,101,0,101,0);
    uh = L2_projection_1D('initial_condition',a,10,mp,P,T);
end

%% basis type
if mp == 1
    basis_type = 101;
elseif mp == 2
    basis_type = 102;
elseif mp == 3
    basis_type = 103;
end

%% generate the finite element information matrices
if mp == 1
    Pb = P;
    Tb = T;
else
    [Pb,Tb] = generate_1D_FEM_partition(P,T,mp);
end

%% generate reference matrices(use reference function)
% The reference cell matrix of the finite element method
% The advantage of using this is that you can count it and then store it so you don't have to count it again afterwards
M1 = generate_1D_matrix_FEM_r(Gauss_coefficient,mp,1,1,basis_type);

%% assemble matrix and right hand vector
A = assemble_1D_matrix_FEM_local(M1,P,T,Tb,mp);
b = generate_1D_vector_FEM('f2_fun',uh,Gauss_coefficient,P,T,Tb,basis_type);
% b = zeros(mp*ng+1,1);
%% boundary treatment
boundary_nodes = [1,1;1,length(T)*mp+1];
[A,b] = boundary_treatment_FEM('boundary_fun',boundary_nodes,A,b,Pb);

%% solution
phi = A\b;

%% error
% ng = length(T);
% error1 = evaluate_error_Lnorm_FEM('temp_exact_fun',1,phi,ng,mp,P,T,Tb,basis_type)
% error2 = evaluate_error_Lnorm_FEM('temp_exact_fun',2,phi,ng,mp,P,T,Tb,basis_type)
% errorinf = evaluate_error_absolute_FEM('temp_exact_fun',phi,P,T,Tb,ng,mp,basis_type)
%% plot
% Output grid midpoints and grid points
% plot_phi_FEM(phi,P,T,Tb,basis_type);
% hold on
