%{
Filename: plot_pc_sweep_metrics.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Normalized transfer metrics of the Gabor models (PT / control, PT - control,
    headroom; pooled means, transfer_metrics.m) vs the number of retained PCs
    per side k (log axis), whole models. Fixing k fixes the rotation's freedom
    (k(k-1)/2 parameters) whatever the model size. Shading is the
    percentile-bootstrap 95% CI over folds; the markers at "full" give the
    full-space value of the same whole model. Colors as Fig 15.

Output (results/figures/07_pc_sweep_metrics/):
    pc_gabor_metrics_vs_k       simple and complex model
    pc_sweep_metrics.csv        every plotted point (k = NaN: full space)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '07_pc_sweep_metrics';
fprintf('plot_pc_sweep_metrics\n');
seed = 35000;

spec = {'ratio', 'PT / control', 1; 'diff', 'PT - control (accuracy)', 0; 'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0};
k_full = 75;
rows = {};
fig = figure('Visible', 'off', 'Position', [100 100 1650 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
K = struct();  Bf = struct();
for j = 1:2
    m = S.models{j};  r = R.pc.(m);  rf = R.cue.(m);
    n = numel(r.k_list);
    K.(m) = struct('k', r.k_list(:)', 'ratio', zeros(1, n), 'ratio_ci', zeros(n, 2), 'diff', zeros(1, n), ...
        'diff_ci', zeros(n, 2), 'headroom', zeros(1, n), 'headroom_ci', zeros(n, 2));
    for ik = 1:n
        seed = seed + 1;
        B = bootstrap_metrics(r.self(ik, :, :), r.pt_rot(ik, :, :), r.ctrl_m(ik, :, :), 3, S, seed);
        for f = {'ratio', 'diff', 'headroom'}
            K.(m).(f{1})(ik) = B.(f{1});
            K.(m).([f{1} '_ci'])(ik, :) = B.([f{1} '_ci']);
        end
        rows(end + 1, :) = {S.label.(m), r.k_list(ik), B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
            B.diff, B.diff_ci(1), B.diff_ci(2), B.headroom, B.headroom_ci(1), B.headroom_ci(2)}; %#ok<SAGROW>
    end
    iN = numel(rf.N_list);  seed = seed + 1;
    Bf.(m) = bootstrap_metrics(rf.self(iN, :, 1, :), rf.pt_rot(iN, :, 1, :), rf.ctrl_m(iN, :, 1, :), 4, S, seed);
    B = Bf.(m);
    rows(end + 1, :) = {S.label.(m), NaN, B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
        B.diff, B.diff_ci(1), B.diff_ci(2), B.headroom, B.headroom_ci(1), B.headroom_ci(2)}; %#ok<SAGROW>
end
for s = 1:3
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, 2);
    for j = 1:2
        m = S.models{j};  k = K.(m);  c = S.c.model.(m);  ci = k.([spec{s, 1} '_ci']);
        fill(ax, [k.k, fliplr(k.k)], [ci(:, 1)', fliplr(ci(:, 2)')], c, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        h(j) = plot(ax, k.k, k.(spec{s, 1}), 'Color', c, 'LineWidth', 1.6, 'Marker', '.', 'MarkerSize', 10);
        v = Bf.(m).(spec{s, 1});  cf = Bf.(m).([spec{s, 1} '_ci']);
        errorbar(ax, k_full * (1 + 0.04 * (j - 1.5)), v, v - cf(1), cf(2) - v, 'o', 'Color', c, 'MarkerFaceColor', c, ...
            'MarkerSize', 5, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, S.k_main, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xline(ax, 60, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');  xlim(ax, [0.9 95]);
    set(ax, 'XTick', [1 2 5 10 20 49 k_full], 'XTickLabel', {'1', '2', '5', '10', '20', '49', 'full'});
    box(ax, 'on');  xlabel(ax, 'Retained PCs per side, k');  ylabel(ax, spec{s, 2});
    if s == 1
        lg = legend(ax, h, {S.label.simple, S.label.complex}, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
    end
end
title(tl, 'Gabor-model cue transfer vs number of PCs per rendering, whole models; "full" = full space, same cells', 'FontSize', 11);
subtitle(tl, 'lines: pooled means; shading: bootstrap 95% CI over folds', 'FontSize', 9);
save_figure(fig, S, out_dir, 'pc_gabor_metrics_vs_k');

T = cell2table(rows, 'VariableNames', {'model', 'k', 'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'headroom', 'headroom_lo', 'headroom_hi'});
save_table(T, S, out_dir, 'pc_sweep_metrics');
