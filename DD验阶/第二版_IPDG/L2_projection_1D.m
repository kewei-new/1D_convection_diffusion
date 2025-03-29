function uh0 = L2_projection_1D(initial_fun,inv_a,ng,mp,h,mid_points)
% evaluate discrete initial condition

Gauss_coefficient = generate_1D_gaussian_quadrature_r(4);
Gp = Gauss_coefficient(:,1);
Gw = Gauss_coefficient(:,2);

uh0 = zeros(mp+1,ng);
for n = 1:ng
    
    Gp_local=Gp*h+mid_points(n);
    ini_value=feval(initial_fun,Gp_local);
    
    for m = 0:mp
        uh0(m+1,n)=inv_a(m+1)*Gw'*(ini_value.*reference_basis(Gp,m));
    end
end

