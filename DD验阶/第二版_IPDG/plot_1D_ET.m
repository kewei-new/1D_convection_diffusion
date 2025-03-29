function  plot_value = plot_1D_ET(mid_points,u1h,phih,phixh,T_end)
% 画中点的值
% 每一行依次为n,E,T,phi,phix
close(gcf);

ng = size(u1h,2);
mp = size(u1h,1)-1;


% 画5个图
plot_value = zeros(3,ng);

for n = 1:ng
    
    u1_local = evaluate_uh_DG(0.5,u1h(:,n),mp);
    phi_local = evaluate_uh_DG(0,phih(:,n),mp);
    phix_local = evaluate_uh_DG(0,phixh(:,n),mp);

    plot_value(:,n) = [u1_local,phi_local,-phix_local]';

end

figure;
t = tiledlayout(2, 2); % 创建2行2列布局
t.Padding = 'compact'; % 紧凑间距

val_name = {'n','phi','-phix'};

for i = 1:3
    nexttile; % 自动填充位置
    plot(mid_points,plot_value(i,:),'-o'); 
    hold on
    plot(mid_points,exact_fun(mid_points,T_end,i),'-*')
    title([val_name(i),'在终止时刻的图像']);
end
