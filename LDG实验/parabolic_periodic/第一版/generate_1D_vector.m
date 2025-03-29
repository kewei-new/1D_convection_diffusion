function f_values = generate_1D_vector(Gauss_coefficient,ng,mp,h,mid_points,inv_mass,t)

f_values = zeros(mp+1,ng);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

for n = 1:ng

    Gp_local=Gp*h+mid_points(n);
    f_local=f_fun(Gp_local,t);
    
    for m = 0:mp
        f_values(m+1,n)=inv_mass(m+1)*Gw'*(f_local.*reference_basis(Gp,m));
    end
end

f_values = [f_values(:,ng),f_values,f_values(:,1)];