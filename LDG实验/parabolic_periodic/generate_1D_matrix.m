function [H,R,L] = generate_1D_matrix(Gauss_coefficient,mp,inv_mass,h)
% ut + ux = uxx + f
% ut + q = qx + f
% q = ux

% the second layer matrix，
% (ut,v) - (ux,v) = -(q,vx) + (q,v)e + (f,v)
% the second layer matrix
% (q,p) = -(u,px) + (u,p)e


% (q,v)e = q(1+)v(1-) - q(-1+)v(-1+)
% (u,p)e = u(1-)p(1-) - u(-1-)v(-1+)


Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

H = zeros(mp+1,mp+1,2);
R = zeros(mp+1,mp+1,2);
L = zeros(mp+1,mp+1,2);

for i = 0:mp
    for j = 0:mp
        
        H(i+1,j+1,1) = 1/h*inv_mass(i+1)*(Gw'*(-reference_basis(Gp,j).*reference_basis_x(Gp,i)) ...
                  -reference_basis(-0.5,j)*reference_basis(-0.5,i));
        R(i+1,j+1,1) = 1/h*inv_mass(i+1)*(reference_basis(-0.5,j)*reference_basis(0.5,i));
        L(i+1,j+1,1) = 0;
            
        H(i+1,j+1,2) = 1/h*inv_mass(i+1)*(Gw'*(-reference_basis(Gp,j).*reference_basis_x(Gp,i)) ...
                  +reference_basis(0.5,j)*reference_basis(0.5,i));
        R(i+1,j+1,2) = 0;
        L(i+1,j+1,2) =-1/h*inv_mass(i+1)*(reference_basis(0.5,j)*reference_basis(-0.5,i));

    end
end


end