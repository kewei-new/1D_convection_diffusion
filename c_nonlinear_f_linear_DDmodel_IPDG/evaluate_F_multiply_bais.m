function right_term = evaluate_F_multiply_bais(f_fun,Gauss_coefficient,ng,mp,P,T)
% evaluate (f,v)
% right_term: shape is (mp+1)*ng; every row is in an element

GP = Gauss_coefficient(:,1);
GW = Gauss_coefficient(:,2);
Gpn = length(GP);

right_term = zeros(mp+1,ng);

for n = 1:ng
    
    local_points = P(:,T(:,n));

    % 这里是进行仿射变换
    h_local = max(local_points) - min(local_points);
    mid_point = mean(local_points);
    Gauss_points_local = GP*h_local + mid_point;

    % 右端项，这里没有对intial_fun进行放缩还原，但前提条件是在[-0.5,0.5]上有意义
    for m = 0:mp 
        for k = 1:Gpn
                        
            right_term(m+1,n) = right_term(m+1,n) + ...
                GW(k)*feval(f_fun,Gauss_points_local(k),0)*reference_basis(GP(k),m,101,0);

        end
    end

end