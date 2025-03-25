function result = evaluate_error_Lnorm_FEM(exact_fun,Lnorm,uh,ng,mp,P,T,Tb,basis_type)

Gauss_coefficient = generate_1D_gaussian_quadrature_r(4);
GP = Gauss_coefficient(:,1);
GW = Gauss_coefficient(:,2);
Gpn = length(GP);

uh = uh(Tb);

result = 0;
for n = 1:ng
    
    local_points = P(:,T(:,n));
    h_local = max(local_points) - min(local_points);
    mid_point = mean(local_points);
    Gauss_points_local = GP*h_local + mid_point;

    uh_local = uh(:,n);

    for k = 1:Gpn
        uh_a = 0;
        % 求解对应Gauss点上的逼近函数值
        for m = 1:mp+1
            uh_a = uh_a + uh_local(m)*basis_reference_FEM(GP(k),0,m,basis_type);
        end
        % Gauss积分求和
         result = result + GW(k)*abs(feval(exact_fun,Gauss_points_local(k))-uh_a)^Lnorm*h_local;
    end

end

result = result^(1/Lnorm);