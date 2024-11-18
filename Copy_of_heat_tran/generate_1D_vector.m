function b = generate_1D_vector(t,Gauss_coefficient,inv_a,mid_points,mp,ng,h)

b = zeros(ng*(mp+1),1);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

for n = 1:ng

    % index = (n-1)*(mp+1);
    Gp_local = Gp*h + mid_points(n);
    for m = 0:mp
        
        % index_local = index+m+1;
        s = m+1;
        b((n-1)*(mp+1)+m+1,1) = inv_a(s,s)*Gw'*(f_fun1(Gp_local,t).*reference_basis(Gp,m,0));

    end

end
