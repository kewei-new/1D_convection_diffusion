function plot_pn_results_1D(x, doping, n, phi, E, outdir, tag)
% 中文美化版 PN 结绘图（单图输出版）

    if nargin < 7
        tag = '默认';
    end

    % ===== 字体 =====
    set(groot, 'defaultAxesFontName', 'Microsoft YaHei');

    lw = 2;
    fs = 12;

    %% ===== 1. 掺杂分布 =====
    fig1 = figure('Position', [200, 200, 700, 500]);
    plot(x, doping, 'LineWidth', lw);
    title('掺杂分布');
    xlabel('位置 x');
    ylabel('掺杂浓度');
    % legend('N_D(x)','Location','best');
    grid on; set(gca,'FontSize',fs);

    saveas(fig1, fullfile(outdir, sprintf('掺杂_%s.png', tag)));
    saveas(fig1, fullfile(outdir, sprintf('掺杂_%s.fig', tag)));

    %% ===== 2. 载流子密度 =====
    fig2 = figure('Position', [200, 200, 700, 500]);
    plot(x, n, 'LineWidth', lw);
    title('载流子密度');
    xlabel('位置 x');
    ylabel('n(x)');
    % legend('电子密度','Location','best');
    grid on; set(gca,'FontSize',fs);

    saveas(fig2, fullfile(outdir, sprintf('载流子_%s.png', tag)));
    % saveas(fig2, fullfile(outdir, sprintf('载流子_%s.fig', tag)));

    %% ===== 3. 电势 =====
    fig3 = figure('Position', [200, 200, 700, 500]);
    plot(x, phi, 'LineWidth', lw);
    title('电势分布');
    xlabel('位置 x');
    ylabel('\phi(x)');
    % legend('电势','Location','best');
    grid on; set(gca,'FontSize',fs);

    saveas(fig3, fullfile(outdir, sprintf('电势_%s.png', tag)));
    % saveas(fig3, fullfile(outdir, sprintf('电势_%s.fig', tag)));

    %% ===== 4. 电场 =====
    fig4 = figure('Position', [200, 200, 700, 500]);
    plot(x, E, 'LineWidth', lw);
    title('电场分布');
    xlabel('位置 x');
    ylabel('E(x)');
    % legend('电场强度','Location','best');
    grid on; set(gca,'FontSize',fs);

    saveas(fig4, fullfile(outdir, sprintf('电场_%s.png', tag)));
    % saveas(fig4, fullfile(outdir, sprintf('电场_%s.fig', tag)));

    fprintf('四张图已保存至 result 文件夹（标签：%s）\n', tag);
end