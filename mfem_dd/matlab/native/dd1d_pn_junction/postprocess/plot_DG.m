function plot_DG(uh,ng,mp,mid_points)

uh_a_list=[];
val_name = inputname(1);
for n = 1:ng

    uh_local = uh(:,n);
    uh_a = 0;
    for m = 0:mp
        uh_a = uh_a + uh_local(m+1)*reference_basis(0,m);
    end

    uh_a_list = [uh_a_list,uh_a];

end

plot(mid_points,uh_a_list,'o-');

% hold on
%
% plot(mid_points,nd(mid_points),'o-');

title([val_name,'在终止时刻的图像']);
