function result = evaluate_error_Lnorm(exact_fun,Lnorm,t,uh,ng,mp,P,T)

Gauss_coefficient = generate_1D_gaussian_quadrature_r(4);
GP = Gauss_coefficient(:,1);
GW = Gauss_coefficient(:,2);
Gpn = length(GP);

uh = reshape(uh,mp+1,ng);

result = 0;
for n = 1:ng
    
    local_points = P(:,T(:,n));
    h_local = max(local_points) - min(local_points);
    mid_point = sum(local_points)/length(local_points);
    Gauss_points_local = GP*h_local + mid_point;

    uh_local = uh(:,n);

    for k = 1:Gpn
        uh_a = 0;
        % 求解对应Gauss点上的逼近函数值
        for m = 0:mp
            uh_a = uh_a + uh_local(m+1)*reference_basis(GP(k),m,101,0);
        end
        % Gauss积分求和
         result = result + GW(k)*abs(feval(exact_fun,Gauss_points_local(k),t)-uh_a)^Lnorm;
    end

end

result = result^(1/Lnorm);