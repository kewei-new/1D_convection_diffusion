function [uh,qh]=solve_1D_LDG(Gauss_coefficient,inv_mass,h,mid_points,mp,ng)
% -u_xx = cos(x) = f , let q=-u_x
% 以下为推导过程，periodic边界和dirichlet边界的flux不一样
% (q_x,v) = (f,v)
% (-u_x,r) = (q,r)
% -(q,v_x) + (q,v)e = (f,v)
% (u,r_x) - (u,r)e = (q,r)

% 取flux(u)： u-
% 取flux(q)： q+ 在第n个区间q(N) =  q(N)- - C_p*(u(N-) - u(N+))
C_p = mp/max(h);

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
f_terms = zeros((mp+1)*ng,1);


%% (q,v) = (ux,v)
% Q = M1*U + b
M1 = sparse((mp+1)*ng,(mp+1)*ng);
b1 = zeros(ng*(mp+1),1);
for n = 1:ng

    if n == 1        
        % M1 = -(u,vx) + u(1/2)*v(1/2)
        % b1 = -u(0)*v(-1/2)
        for alpha = 0:mp %test
            for beta = 0:mp %trial
                M1(alpha+1,beta+1) =inv_mass(alpha+1)*h(n)^(-1)*(-Gw'*(reference_basis(Gp,beta).*reference_basis_x(Gp,alpha))...
                                              +reference_basis(1/2,beta)*reference_basis(1/2,alpha)); 
            end
            b1(alpha+1) =-inv_mass(alpha+1)*h(n)^(-1)*LDG_bound(mid_points(1)-h(1)/2)*reference_basis(-1/2,alpha);
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
            b1((n-1)*(mp+1)+(alpha+1)) = inv_mass(alpha+1)*h(n)^(-1)*LDG_bound(mid_points(ng)+h(ng)/2)*reference_basis(1/2,alpha);
        end 
    end

end

%% 求(qx,v)
M2 = sparse((mp+1)*ng,(mp+1)*ng);
M3 = sparse((mp+1)*ng,(mp+1)*ng);
b2 = zeros(ng*(mp+1),1);
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
           b2((n-1)*(mp+1)+(alpha+1)) = C_p*LDG_bound(mid_points(ng)+h(ng)/2)*reference_basis(1/2,alpha); 
        end              
    end
end

%% 求(f,v)
for n = 1:ng
    Gp_lcoal = Gp*h(n) + mid_points(n);
    for alpha = 0:mp
        f_terms((n-1)*(mp+1)+(alpha+1),1) = Gw'*(f_fun(Gp_lcoal).*reference_basis(Gp,alpha))*h(n);
    end
end


%% 综上
% (qx,v) = M2*Q + M3*U + b2
% Q = M1*U + b1
% M2*(M1*U + b1) + M3*U + b2= (f,v)
% (M2*M1+M3)*U + M2*b1 + b2 = f_term


f_terms = f_terms - M2*b1 - b2;
uh = (M2*M1+M3)\f_terms;
qh = M1*uh + b1;

uh = reshape(uh,mp+1,ng);
qh = reshape(qh,mp+1,ng);
