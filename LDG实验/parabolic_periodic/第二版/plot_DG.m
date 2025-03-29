function plot_DG(exact_fun,t,uh,ng,mp,mid_points)

uh_a_list=[];
exact_uh = [];

for n = 1:ng
    
    uh_local = uh(:,n);
    uh_a = 0;
    for m = 0:mp
        uh_a = uh_a + uh_local(m+1)*reference_basis(0,m);
    end

    uh_a_list = [uh_a_list,uh_a];
    exact_uh = [exact_uh,feval(exact_fun,mid_points(n),t)];

end

plot(mid_points,uh_a_list,'b-',mid_points,exact_uh,'o-');
