function [M1,M2,Mn]=assemble_LDG_matrix(Gauss_coefficient,inv_mass,mp,ng,h)
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
M1 = sparse((mp+1)*ng,(mp+1)*ng);
C_p = mp/max(h);

for n = 1:ng

    if n == 1
        % M1 = -(u,vx) + u(1/2)*v(1/2)
        % b1 = -u(0)*v(-1/2)
        for alpha = 0:mp %test
            for beta = 0:mp %trial
                M1(alpha+1,beta+1) =inv_mass(alpha+1)*h(n)^(-1)*(-Gw'*(reference_basis(Gp,beta).*reference_basis_x(Gp,alpha))...
                                              +reference_basis(1/2,beta)*reference_basis(1/2,alpha));
            end
        end

    elseif n <ng
        % M1 = -(u,vx) + uj(1/2)v(1/2) - uj-1(1/2)v(-1/2)
        for alpha = 0:mp %test
            for beta = 0:mp %trial
                M1((n-1)*(mp+1)+(alpha+1),(n-1)*(mp+1)+(beta+1)) =inv_mass(alpha+1)*h(n)^(-1)*...
                    (-Gw'*(reference_basis(Gp,beta).*reference_basis_x(Gp,alpha))...
                     +reference_basis(1/2,beta)*reference_basis(1/2,alpha));
                M1((n-1)*(mp+1)+(alpha+1),(n-2)*(mp+1)+(beta+1)) =inv_mass(alpha+1)*h(n)^(-1)*...
                    (-reference_basis(1/2,beta)*reference_basis(-1/2,alpha));
            end
        end
    else
        % M1 = -(u,vx) - uj-1(1/2)v(-1/2)
        % b1 = u(pi)*v(1/2)
        for alpha = 0:mp %test
            for beta = 0:mp %trial
                M1((n-1)*(mp+1)+(alpha+1),(n-1)*(mp+1)+(beta+1)) =inv_mass(alpha+1)*h(n)^(-1)*...
                    (-Gw'*(reference_basis(Gp,beta).*reference_basis_x(Gp,alpha)));
                M1((n-1)*(mp+1)+(alpha+1),(n-2)*(mp+1)+(beta+1)) =inv_mass(alpha+1)*h(n)^(-1)*...
                    (-reference_basis(1/2,beta)*reference_basis(-1/2,alpha));
            end
        end
    end
end
%% 求(qx,v)
M2 = sparse((mp+1)*ng,(mp+1)*ng);
M3 = sparse((mp+1)*ng,(mp+1)*ng);
for n = 1:ng

    if n < ng
        % M2 = -(q,vx) + qj+1(-1/2)v(1/2) - qj(-1/2)v(-1/2)
        for alpha = 0:mp %test
            for beta = 0:mp %trial
                M2((n-1)*(mp+1)+(alpha+1),(n-1)*(mp+1)+(beta+1)) =-Gw'*(reference_basis(Gp,beta).*reference_basis_x(Gp,alpha))...
                     -reference_basis(-1/2,beta)*reference_basis(-1/2,alpha);
                M2((n-1)*(mp+1)+(alpha+1),(n)*(mp+1)+(beta+1)) =reference_basis(-1/2,beta)*reference_basis(1/2,alpha);
            end
        end
    else
        % M2 = -(q,vx) - q(-1/2)v(-1/2) + q(1/2)v(1/2)
        % M3 = -C_p*u(1/2)v(1/2)
        % b2 = C_p*u(pi)*v(1/2)
        for alpha = 0:mp %test
            for beta = 0:mp %trial
                M2((n-1)*(mp+1)+(alpha+1),(n-1)*(mp+1)+(beta+1)) =-Gw'*(reference_basis(Gp,beta).*reference_basis_x(Gp,alpha))...
                     -reference_basis(-1/2,beta)*reference_basis(-1/2,alpha)...
                     +reference_basis(1/2,beta)*reference_basis(1/2,alpha);
                M3((n-1)*(mp+1)+(alpha+1),(n-1)*(mp+1)+(beta+1)) =-C_p*reference_basis(1/2,beta)*reference_basis(1/2,alpha);
            end
        end
    end
end

Mn = M2*M1+M3;
end