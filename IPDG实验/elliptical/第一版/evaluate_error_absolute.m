function result = evaluate_error_absolute(exact_fun,t,uh,h,mid_points,ng,mp)

result = 0;
for n = 1:ng
    
    uh_local = uh(:,n);
    uh_a = 0;
    for m = 0:mp
        uh_a = uh_a + uh_local(m+1)*reference_basis(0,m);
    end
    temp = abs(feval(exact_fun,mid_points(n),t)-uh_a);

    result = max(result,temp);
    % if result == temp
    %     n
    % end
end