function result = evaluate_error_absolute(exact_fun,t,uh,P,T,ng,mp)
% 只计算左端点
uh = reshape(uh,mp+1,ng);

result = 0;
for n = 1:ng
    
    local_points = P(:,T(:,n));
    h_local = max(local_points) - min(local_points);
    mid_point = sum(local_points)/length(local_points);
    refer_points = (min(local_points)-mid_point)/h_local;
    uh_local = uh(:,n);
    
    uh_a = 0;
    for m = 0:mp
        uh_a = uh_a + uh_local(m+1)*reference_basis(refer_points,m,101,0);
    end
    temp = abs(feval(exact_fun,min(local_points),t)-uh_a);

    result = max(result,temp);

end