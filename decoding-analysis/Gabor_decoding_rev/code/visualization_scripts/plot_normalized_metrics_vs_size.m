%{
Filename: plot_normalized_metrics_vs_size.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Normalized transfer metrics of the Gabor models (pooled means,
    transfer_metrics.m) as a function of the number of model cells, full space:
    PT / control, PT - control, headroom and its denominator self - control,
    simple and complex model on the same axes (colors as Fig 15). Shading is the
    percentile-bootstrap 95% CI over cell draws (folds for the whole model, the
    last point of each line). The dotted line marks the Fig 15a size (48).

Output (results/figures/05_normalized_metrics_vs_size/):
    gabor_metrics_vs_cells            the figure described above
    normalized_metrics_vs_size.csv    every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '05_normalized_metrics_vs_size';
fprintf('plot_normalized_metrics_vs_size\n');

spec = {'ratio', 'PT / control', 1; 'diff', 'PT - control', 0; 'headroom', 'Headroom', 0; 'denom', 'self - control', 0};
seed = 7000;
rows = {};
V = struct();
for j = 1:2
    m = S.models{j};  r = R.cue.(m);
    x = r.N_list(:)';  nN = numel(x);
    V.(m) = struct('x', x, 'ratio', nan(1, nN), 'ratio_ci', nan(nN, 2), 'diff', nan(1, nN), 'diff_ci', nan(nN, 2), ...
        'headroom', nan(1, nN), 'headroom_ci', nan(nN, 2), 'denom', nan(1, nN), 'denom_ci', nan(nN, 2));
    for iN = 1:nN
        D = r.n_draws(iN);  ud = 3 + (D == 1);
        seed = seed + 1;
        B = bootstrap_metrics(r.self(iN, :, 1:D, :), r.pt_rot(iN, :, 1:D, :), r.ctrl_m(iN, :, 1:D, :), ud, S, seed);
        for f = {'ratio', 'diff', 'headroom'}
            V.(m).(f{1})(iN) = B.(f{1});
            V.(m).([f{1} '_ci'])(iN, :) = B.([f{1} '_ci']);
        end
        V.(m).denom(iN) = B.self - B.ctrl;
        rows(end + 1, :) = {S.label.(m), x(iN), B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
            B.diff, B.diff_ci(1), B.diff_ci(2), B.headroom, B.headroom_ci(1), B.headroom_ci(2), B.self - B.ctrl, B.n_units}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'model', 'N', 'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'headroom', 'headroom_lo', 'headroom_hi', 'self_minus_ctrl', 'n_units_boot'});
save_table(T, S, out_dir, 'normalized_metrics_vs_size');

fig = figure('Visible', 'off', 'Position', [100 100 1300 800], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:size(spec, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, 2);
    for j = 1:2
        m = S.models{j};  v = V.(m);  c = S.c.model.(m);
        ci = v.([spec{s, 1} '_ci']);
        if ~all(isnan(ci(:)))
            fill(ax, [v.x, fliplr(v.x)], [ci(:, 1)', fliplr(ci(:, 2)')], c, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        h(j) = plot(ax, v.x, v.(spec{s, 1}), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3, 'LineWidth', 1.5);
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, S.fig15a_N, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xlim(ax, [0 200]);  box(ax, 'on');  ylabel(ax, spec{s, 2});  xlabel(ax, 'Number of model cells');
    if s == 1, legend(ax, h, {S.label.simple, S.label.complex}, 'Location', 'best'); end
end
title(tl, 'Gabor-model cue transfer vs number of cells, full space (shading: bootstrap 95% CI over cell draws; dotted line: 48 cells)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'gabor_metrics_vs_cells');
