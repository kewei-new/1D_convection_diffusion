function plot_cv_results(bias, Cqs, outdir)
%PLOT_CV_RESULTS Plot and save C-V curve.

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    set(groot, 'defaultAxesFontName', 'Microsoft YaHei');

    lw = 2;
    fs = 12;

    fig = figure('Position', [200, 200, 700, 500]);
    plot(bias, Cqs, '-o', 'LineWidth', lw, 'MarkerSize', 6);
    title('准静态电容-电压曲线');
    xlabel('偏压 V');
    ylabel('准静态电容 C_{qs}');
    grid on;
    set(gca, 'FontSize', fs);

    saveas(fig, fullfile(outdir, 'CV曲线.png'));
    saveas(fig, fullfile(outdir, 'CV曲线.fig'));
end