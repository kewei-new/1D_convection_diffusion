function plot_iv_qv_results(bias, current, Qmag, Wdep, outdir)
%PLOT_IV_QV_RESULTS Plot and save I-V, Q-V, Wdep-V curves.

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    set(groot, 'defaultAxesFontName', 'Microsoft YaHei');

    lw = 2;
    fs = 12;

    % I-V
    fig1 = figure('Position', [200, 200, 700, 500]);
    plot(bias, current, '-o', 'LineWidth', lw, 'MarkerSize', 6);
    title('电流-电压曲线');
    xlabel('偏压 V');
    ylabel('右端接触电流');
    grid on;
    set(gca, 'FontSize', fs);
    saveas(fig1, fullfile(outdir, 'IV曲线.png'));
    saveas(fig1, fullfile(outdir, 'IV曲线.fig'));

    % Q-V
    fig2 = figure('Position', [200, 200, 700, 500]);
    plot(bias, Qmag, '-o', 'LineWidth', lw, 'MarkerSize', 6);
    title('电荷-电压曲线');
    xlabel('偏压 V');
    ylabel('空间电荷指标 Q_{mag}');
    grid on;
    set(gca, 'FontSize', fs);
    saveas(fig2, fullfile(outdir, 'QV曲线.png'));
    saveas(fig2, fullfile(outdir, 'QV曲线.fig'));

    % Wdep-V
    fig3 = figure('Position', [200, 200, 700, 500]);
    plot(bias, Wdep, '-o', 'LineWidth', lw, 'MarkerSize', 6);
    title('耗尽层宽度-电压曲线');
    xlabel('偏压 V');
    ylabel('耗尽层宽度 W_{dep}');
    grid on;
    set(gca, 'FontSize', fs);
    saveas(fig3, fullfile(outdir, 'WdepV曲线.png'));
    saveas(fig3, fullfile(outdir, 'WdepV曲线.fig'));
end