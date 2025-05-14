function plot_solution(exact_fun,t,uh,ng,mp,P,T)

uh = reshape(uh,mp+1,ng);

x = 1/2*(P(:,T(1,:))+P(:,T(2,:)));
uh_a_list=[];
exact_uh = [];

for n = 1:ng
    
    local_points = P(:,T(:,n));
    mid_point = sum(local_points)/length(local_points);
    uh_local = uh(:,n);
    
    uh_a = 0;
    for m = 0:mp
        uh_a = uh_a + uh_local(m+1)*reference_basis(0,m,101,0);
    end
    uh_a_list = [uh_a_list,uh_a];

    exact_uh = [exact_uh,feval(exact_fun,mid_point,t)];

end

plot(x,uh_a_list,'-',x,exact_uh,'o');
