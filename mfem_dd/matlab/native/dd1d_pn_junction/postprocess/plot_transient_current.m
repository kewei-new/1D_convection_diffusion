function plot_transient_current(time, current, outdir)
%PLOT_TRANSIENT_CURRENT Plot and save transient current curve.

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    set(groot, 'defaultAxesFontName', 'Microsoft YaHei');

    lw = 2;
    fs = 12;

    fig = figure('Position', [200, 200, 700, 500]);
    plot(time, current, '-o', 'LineWidth', lw, 'MarkerSize', 4);
    title('瞬态电流响应');
    xlabel('时间 t');
    ylabel('右端接触电流');
    grid on;
    set(gca, 'FontSize', fs);

    saveas(fig, fullfile(outdir, '瞬态电流.png'));
    saveas(fig, fullfile(outdir, '瞬态电流.fig'));
end