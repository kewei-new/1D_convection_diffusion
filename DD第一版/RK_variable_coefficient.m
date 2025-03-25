function uh=RK_variable_coefficient(Gauss_coefficient,M1,inv_a,uh0,t,dt,mp,mt,P,T,trial_basis_type,test_basis_type)
% RK# / RK4
% result(:,i): the i th Rk iteration's solution

result = zeros(length(uh0),mt-1);
ng = length(T);
matrix_D = generate_1D_boundary_multiple(mp,trial_basis_type,0,test_basis_type,0);
xmu = 0.75;
if mt==3

    % io = 1
    [phi,Tb] = solve_poisson_1D_FEM(mp,P,T,Gauss_coefficient,uh0);
    M2 = generate_1D_variable_coefficient_matrix_local('evaluate_local_FEM',inv_a,Gauss_coefficient,P,T,trial_basis_type,0,test_basis_type,1,Tb,phi);
    [L,H,R] = generate_1D_variable_coefficient_flux_local(inv_a,matrix_D,phi,mp,ng,P,T,Tb);
    A = xmu*assemble_1D_matrix_with_variable(ng,mp,M2,L,H,R);
    b = zeros(ng*(mp+1),1);
    % b = generate_1D_vector('f1_fun',t,Gauss_coefficient,ng,mp,P,T);
    % boundary_nodes = [1;1;1];
    % [A,b] = boundary_treatment_1D(boundary_nodes,A+M1,b,ng,mp);
    result(:,1) = uh0 + dt*rhs_variable_coefficient(uh0,A+M1,b);

    % io = 2
    [phi,Tb] = solve_poisson_1D_FEM(mp,P,T,Gauss_coefficient,result(:,1));
    M2 = generate_1D_variable_coefficient_matrix_local('evaluate_local_FEM',inv_a,Gauss_coefficient,P,T,trial_basis_type,0,test_basis_type,1,Tb,phi);
    [L,H,R] = generate_1D_variable_coefficient_flux_local(inv_a,matrix_D,phi,mp,ng,P,T,Tb);
    A = xmu*assemble_1D_matrix_with_variable(ng,mp,M2,L,H,R);
    b = zeros(ng*(mp+1),1);
    % b = generate_1D_vector('f1_fun',t+dt/2,Gauss_coefficient,ng,mp,P,T);
    % boundary_nodes = [1;1;1];
    % [A,b] = boundary_treatment_1D(boundary_nodes,A+M1,b,ng,mp);
    result(:,2) = 3/4*uh0 + 1/4*result(:,1) + 1/4*dt*rhs_variable_coefficient(result(:,1),A+M1,b);

    % io = 3
    [phi,Tb] = solve_poisson_1D_FEM(mp,P,T,Gauss_coefficient,result(:,2));
    M2 = generate_1D_variable_coefficient_matrix_local('evaluate_local_FEM',inv_a,Gauss_coefficient,P,T,trial_basis_type,0,test_basis_type,1,Tb,phi);
    [L,H,R] = generate_1D_variable_coefficient_flux_local(inv_a,matrix_D,phi,mp,ng,P,T,Tb);
    A = xmu*assemble_1D_matrix_with_variable(ng,mp,M2,L,H,R);
    b = zeros(ng*(mp+1),1);
    % b = generate_1D_vector('f1_fun',t+dt,Gauss_coefficient,ng,mp,P,T);
    % boundary_nodes = [1;1;1];
    % [A,b] = boundary_treatment_1D(boundary_nodes,A+M1,b,ng,mp);
    uh = 1/3*uh0 + 2/3*result(:,2) + 2/3*dt*rhs_variable_coefficient(result(:,2),A+M1,b);
    

elseif mt==4
    % 没有完成

end
