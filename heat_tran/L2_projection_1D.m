function uh0 = L2_projection_1D(initial_fun,Gauss_coefficient,inv_a,mid_points,mp,ng,h,t)

uh0 = zeros(mp+1,ng);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

for n = 1:ng

    for m = 0:mp
        
        uh0(m+1,n) = inv_a(m+1,m+1)*Gw'*(feval(initial_fun,Gp*h+mid_points(n),t).*reference_basis(Gp,m,0));

    end
end