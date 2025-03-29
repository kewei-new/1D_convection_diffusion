function [f_values,b2] = generate_1D_LDG_vector(Gauss_coefficient,inv_mass,ng,mp,h,mid_points,uh,M1,t)

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
alpha = mp/h;
% 
b2 = zeros(ng*(mp+1),1);
n = 1;
g_a = LDG_bound(mid_points(n)-h/2,t);
for m = 0:mp
    b2((n-1)*(mp+1)+1+m,1) =-1/h*inv_mass(m+1)*g_a*reference_basis(-1/2,m);
end
n = ng;
g_b = LDG_bound(mid_points(n)+h/2,t);
for m = 0:mp
    b2((n-1)*(mp+1)+1+m,1) = 1/h*inv_mass(m+1)*g_b*reference_basis(1/2,m);
end
% 
b1 = zeros(ng*(mp+1),1);
n = ng;
for m = 0:mp
    b1((n-1)*(mp+1)+1+m,1) = alpha*g_b*reference_basis(1/2,m);
end
% 
f_values = -b1 - M1*b2;
f_values = reshape(f_values,mp+1,ng);

% 
for n = 1:ng

    if abs(t)<1e-12
        f_values = zeros(mp+1,ng);
        break
    end

    Gp_local=Gp*h+mid_points(n);
    f2 = sin(t)*cos(Gp_local)+sin(Gp_local)*cos(t);
    n_value = zeros(length(Gp),1);

    for m = 0:mp
        n_value = n_value + uh(m+1,n)*reference_basis(Gp,m);
    end

    for m = 0:mp
        f_values(m+1,n)=f_values(m+1,n)+Gw'*((n_value -  f2).*reference_basis(Gp,m))*h;
    end

    
end
