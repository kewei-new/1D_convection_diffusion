function [solution,X_points] = plot_1D_DG(uh)
% 画图并输出结果
global nx hx Xm dimPk
global phiG
global lambda

NumGLP = size(phiG,1);
solution = zeros(NumGLP,nx);
X_points = zeros(NumGLP,nx);

for i = 1:nx

    X_points(:,i) = lambda(:)*hx + Xm(i);
    for d = 1:dimPk

        solution(:,i) = solution(:,i) + uh(d,i)*phiG(:,d);
        
    end
end

x_temp = reshape(X_points,[NumGLP*nx,1]);
u_temp = reshape(solution,[NumGLP*nx,1]);

% 将数据按x坐标排序
[x_sorted, sort_idx] = sort(x_temp);
u_sorted = u_temp(sort_idx);

% 绘制曲线图
figure;
hold on
plot(x_sorted, u_sorted, 'b-', 'LineWidth', 1.5);
xlabel('x');
ylabel('u');
title('Numerical Solution');
grid on;