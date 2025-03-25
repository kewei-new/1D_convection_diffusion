function f_values = generate_1D_LDG_vector(f_values,Gauss_coefficient,ng,mp,h,mid_points,uh)

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);
xe = 0.1602;
epsilon = 11.7*8.85418;

for n = 1:ng

    Gp_local=Gp*h+mid_points(n);
    xnd = nd(Gp_local);
    n_value = zeros(length(Gp),1);

    for m = 0:mp
        n_value = n_value + uh(m+1,n)*reference_basis(Gp,m);
    end

    for m = 0:mp
        f_values(m+1,n)=f_values(m+1,n)+Gw'*(xe/epsilon*(n_value - xnd).*reference_basis(Gp,m))*h;
    end

    
end
