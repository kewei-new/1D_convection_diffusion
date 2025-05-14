function b_vector = generate_1D_vector(Gauss_coefficient,t,inv_a,ng,mp,h,mid_points)

b_vector = zeros(mp+1,ng);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

for n = 1:ng

    Gp_local_n = Gp*h+ mid_points(n);

    for m = 0:mp

        b_vector(m+1,n) = inv_a(m+1)*Gw'*(f1_fun(Gp_local_n,t).*reference_basis(Gp,m));

    end

end