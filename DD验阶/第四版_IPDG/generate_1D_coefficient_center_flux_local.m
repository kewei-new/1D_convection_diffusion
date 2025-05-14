function [H,R,L]=generate_1D_coefficient_center_flux_local(Gauss_coefficient,inv_a,uh,E_point,Tb,mp,h,n,xmu,basis_type)

H = zeros(mp+1,mp+1);
R = zeros(mp+1,mp+1);
L = zeros(mp+1,mp+1);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

if basis_type == 101
    m_FEM = 1;
else
    m_FEM = 2;
end


% E
uh_local = uh(Tb(:,n));
temp_value = zeros(length(Gp),m_FEM+1);
for m = 1:m_FEM+1  
    temp_value(:,m) = reference_basis_FEM(Gp,1,m,basis_type)*h^(-1); 
end
E_value = -temp_value*uh_local;

for alpha = 0:mp
    % test
    for beta = 0:mp
        % trial
        % H(alpha+1,beta+1) = inv_a(alpha+1)*(Gw'*(E_value.*reference_basis(Gp,beta,0).*reference_basis(Gp,alpha,1)));
        H(alpha+1,beta+1) =-xmu*inv_a(alpha+1)*(Gw'*(E_value.*reference_basis(Gp,beta).*reference_basis_x(Gp,alpha)) ...
                           -1/2*(E_point(2,n)*reference_basis(0.5,beta)*reference_basis(0.5,alpha)-E_point(1,n)*reference_basis(-0.5,beta)*reference_basis(-0.5,alpha)));
        R(alpha+1,beta+1) = xmu*inv_a(alpha+1)*1/2*E_point(1,n+1)*reference_basis(-0.5,beta)*reference_basis(0.5,alpha);
        L(alpha+1,beta+1) =-xmu*inv_a(alpha+1)*1/2*E_point(2,n-1)*reference_basis(0.5,beta)*reference_basis(-0.5,alpha);

    end
end

H = H/h;
R = R/h;
L = L/h;