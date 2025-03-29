function plot_DG(uh,ng,mp,mid_points)

uh_a_list=[];
val_name = inputname(1);
for n = 1:ng
    
    uh_local = uh(:,n);
    uh_a = 0;
    for m = 0:mp
        uh_a = uh_a + uh_local(m+1)*reference_basis(-1/2,m);
    end

    uh_a_list = [uh_a_list,uh_a];

end

plot(mid_points,uh_a_list,'b-');
% 
% hold on 

% plot(mid_points,exact_fun_x(mid_points,0));
title([val_name,'在终止时刻的图像']);
