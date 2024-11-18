function [phi,E_value] = solve_poisson_1D_FEM(mp,ng,P,T,Pb,Tb,Gauss_coefficient,h,mid_points,uh)
% solve \phi_{xx} = e/{\epsilon}*(n - nd)
% we know n(x,0), and \phi satisfy the dirichlet boundary condition,so we
% could get \phi0, and after in the io(1,2,3,4)th step in Rk, we need this
% proceding to get the right hand side f(x+h/2,t+h/2);

%% nargin
% if nargin < 1
%     left = 0;
%     right = 0.6;
%     mp = 1;
%     ng = 5;
%     P = linspace(left,right,ng+1);
%     T=[linspace(1,ng,ng);linspace(2,ng+1,ng)];
%     Gauss_coefficient = generate_1D_gaussian_quadrature_r(mp);
% 
%     h = max(diff(P));
%     mid_points = (P(T(1,:))+P(T(2,:)))/2;
%     inv_a = generate_1D_inverse_mass_matrix(mp);
%     uh=L2_projection_1D(inv_a,ng,mp,h,mid_points);
% end

%% basis type
if mp == 1
    basis_type = 101;
else
    basis_type = 102;
end
Tb = Tb(:,2:ng+1);
% 
% %% generate the finite element information matrices
% if basis_type == 101
%     Pb = P;
%     Tb = T;
% else
%     [Pb,Tb] = generate_1D_FEM_partition(P,T);
% end

%% generate reference matrices(use reference function)
% The reference cell matrix of the finite element method
% The advantage of using this is that you can count it and then store it so you don't have to count it again afterwards
% M1 = generate_1D_matrix_FEM_r(1,1,basis_type);

%% assemble matrix and right hand vector
A = assemble_1D_matrix_FEM_local(Tb,h,basis_type);
b = generate_1D_vector_FEM(Gauss_coefficient,uh,Tb,h,mid_points,mp,ng,basis_type);
% b = zeros(mp*ng+1,1);
%% boundary treatment
boundary_nodes = [1,1;1,Tb(end,end)];
[A,b] = boundary_treatment_FEM(boundary_nodes,A,b,Pb);

%% solution
phi = A\b;
E_value = evaluate_E(phi,Tb,basis_type,ng,h);

% 周期边界
E_value = [E_value(:,ng),E_value,E_value(:,1)];

%% error
% ng = length(T);
% error1 = evaluate_error_Lnorm_FEM('temp_exact_fun',1,phi,ng,mp,P,T,Tb,basis_type)
% error2 = evaluate_error_Lnorm_FEM('temp_exact_fun',2,phi,ng,mp,P,T,Tb,basis_type)
% errorinf = evaluate_error_absolute_FEM('temp_exact_fun',phi,P,T,Tb,ng,mp,basis_type)
%% plot
% Output grid midpoints and grid points
% hold on
% plot_phi_FEM(phi,P,T,Tb,basis_type);
plot(P(1:ng),E_value(1,2:ng+1));

