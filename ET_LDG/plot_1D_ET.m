function  plot_value = plot_1D_ET(mid_points,u1h,u2h,phih,phixh)
% 画中点的值
% 每一行依次为n,E,T,phi,phix
close(gcf);

ng = size(u1h,2);
mp = size(u1h,1)-1;

xm = 0.26*0.9109;
xe = 0.1602;
xk = 0.138046e-4;

% 画5个图
plot_value = zeros(5,ng);

for n = 1:ng
    
    u1_local = evaluate_uh_DG(0.5,u1h(:,n),mp);
    u2_local = evaluate_uh_DG(0.5,u2h(:,n),mp);

    n_local = u1_local/xe;
    E_local = u2_local/u1_local *xm*xe;
    Tem_local = 2./5./xk*(sqrt(abs(1.+10*E_local/3))-1.);
    phi_local = evaluate_uh_DG(0,phih(:,n),mp);
    phix_local = evaluate_uh_DG(0,phixh(:,n),mp);

    plot_value(:,n) = [n_local,E_local,Tem_local,phi_local,-phix_local]';

end

figure;
t = tiledlayout(2, 3); % 创建2行3列布局
t.Padding = 'compact'; % 紧凑间距

val_name = {'n','E','Temperature','phi','-phix'};

for i = 1:5
    nexttile; % 自动填充位置
    plot(mid_points,plot_value(i,:),'-o'); 
    title([val_name(i),'在终止时刻的图像']);
end
