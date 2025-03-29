function f_values = generate_1D_LDG_vector(f_values,Gauss_coefficient,ng,mp,h,mid_points)

Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

for n = 1:ng

    Gp_local=Gp*h+mid_points(n);

    for m = 0:mp
        f_values(m+1,n)=f_values(m+1,n)+Gw'*(f_fun(Gp_local).*reference_basis(Gp,m))*h;
    end
end
