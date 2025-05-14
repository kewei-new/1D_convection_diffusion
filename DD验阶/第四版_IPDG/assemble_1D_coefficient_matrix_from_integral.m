function M = assemble_1D_coefficient_matrix_from_integral(Gauss_coefficient,inv_a,uh,Tb,h,mp,n,xmu)

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
M = zeros(mp+1,mp+1);

if mp == 1
    basis_type = 101;
    m_FEM = 1;
else
    basis_type = 102;
    m_FEM = 2;
end

% E
uh_local = uh(Tb(:,n));
temp_value = zeros(length(Gp),m_FEM+1);
for m = 1:m_FEM+1  
    temp_value(:,m) = reference_basis_FEM(Gp,1,m,basis_type)*h^(-1); 
end
E_value = -temp_value*uh_local;

% inv_a * int (En)v_x
for alpha = 0:mp
    % test(i)
    for beta = 0:mp
        % trial(j)                
        M(alpha+1,beta+1) =xmu*inv_a(alpha+1)*Gw'*(E_value.*reference_basis(Gp,beta,0).*reference_basis(Gp,alpha,1)); 
    end
end

M = 1/h*M;