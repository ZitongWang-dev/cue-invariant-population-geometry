%{
Filename: plot_normalized_metric_bars.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Normalized transfer metrics of the Gabor models as bar charts, all from
    pooled means (transfer_metrics.m): PT / control, PT - control, headroom =
    (PT - control) / (self - control), and the headroom denominator
    self - control (shown because the models are near-noiseless; here the
    denominator stays large, 0.66 or more, so headroom is well defined).
    Error bars are percentile-bootstrap 95% CIs over the resampling unit: cell
    draws, or folds when a whole model is used once; the denominator panel has
    none.

Output (results/figures/04_normalized_metric_bars/):
    gabor_metrics_by_design        groups: 48 cells (Fig 15a), whole models (Fig 15b), and
                                   whole models at k = 5; bars: simple, complex
    gabor_metrics_per_pair_48      per rendering pair at 48 cells (Fig 15a design)
    normalized_metric_bars.csv     every plotted number
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '04_normalized_metric_bars';
fprintf('plot_normalized_metric_bars\n');
rows = {};
seed = 500;
labels = cellfun(@(m) S.label.(m), S.models, 'UniformOutput', false);
colors = [S.c.model.simple; S.c.model.complex];

%% A. by design: 48 cells, whole model (full space), whole model at k = 5
groups = {'48 cells (Fig 15a)', 'whole model (Fig 15b)', sprintf('whole model, k = %d', S.k_main)};
B = cell(numel(groups), 2);
for j = 1:2
    m = S.models{j};  r = R.cue.(m);  rk = R.pc.(m);
    for g = 1:2
        if g == 1, iN = find(r.N_list == S.fig15a_N); else, iN = numel(r.N_list); end
        D = r.n_draws(iN);  ud = 3 + (D == 1);                  % draws, or folds when one draw
        [B{g, j}, seed] = boot(r.self(iN, :, 1:D, :), r.pt_rot(iN, :, 1:D, :), r.ctrl_m(iN, :, 1:D, :), ud, S, seed);
        rows = add_row(rows, 'full', labels{j}, 'all pairs', r.N_list(iN), B{g, j});
    end
    ik = find(rk.k_list == S.k_main);
    [B{3, j}, seed] = boot(rk.self(ik, :, :), rk.pt_rot(ik, :, :), rk.ctrl_m(ik, :, :), 3, S, seed);
    rows = add_row(rows, sprintf('k=%d', S.k_main), labels{j}, 'all pairs', rk.N, B{3, j});
end
bar_figure(S, B, groups, labels, colors, 'Gabor-model cue transfer: 48 cells, whole models and whole models at k = 5', ...
    out_dir, 'gabor_metrics_by_design');

%% B. per rendering pair at 48 cells
pairs = {'EC-AC', 'EC-EX', 'AC-EX'};
B = cell(numel(pairs), 2);
for j = 1:2
    m = S.models{j};  r = R.cue.(m);
    iN = find(r.N_list == S.fig15a_N);  D = r.n_draws(iN);  ud = 3 + (D == 1);
    for g = 1:numel(pairs)
        P = strcmp(r.pair_group, pairs{g});
        [B{g, j}, seed] = boot(r.self(iN, P, 1:D, :), r.pt_rot(iN, P, 1:D, :), r.ctrl_m(iN, P, 1:D, :), ud, S, seed);
        rows = add_row(rows, 'full (48 cells)', labels{j}, pairs{g}, S.fig15a_N, B{g, j});
    end
end
bar_figure(S, B, pairs, labels, colors, 'Gabor-model cue transfer per rendering pair, full space, 48 cells (Fig 15a design)', ...
    out_dir, 'gabor_metrics_per_pair_48');

%% table
T = cell2table(rows, 'VariableNames', {'space', 'model', 'set', 'n_cells', ...
    'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', 'headroom', 'headroom_lo', 'headroom_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'self_minus_ctrl', 'n_units_boot'});
save_table(T, S, out_dir, 'normalized_metric_bars');

%% ===========================================================================
function [B, seed] = boot(Sa, Pa, Ca, unit_dim, S, seed)
% bootstrap_metrics with a fresh seed per call
seed = seed + 1;
B = bootstrap_metrics(Sa, Pa, Ca, unit_dim, S, seed);
end

function rows = add_row(rows, space, model, set, n, B)
rows(end + 1, :) = {space, model, set, n, B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
    B.headroom, B.headroom_ci(1), B.headroom_ci(2), B.diff, B.diff_ci(1), B.diff_ci(2), B.self - B.ctrl, B.n_units};
end

function bar_figure(S, B, groups, series_names, colors, ttl, out_dir, name)
% B{group, series}; four panels: PT/control, PT - control, headroom, self - control
spec = {'ratio', 'PT / control', 1; 'diff', 'PT - control (accuracy)', 0; ...
        'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0; 'denom', 'self - control (headroom denominator)', 0};
fig = figure('Visible', 'off', 'Position', [100 100 1800 480], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
ns = numel(series_names);  w = 0.7 / ns;
for s = 1:size(spec, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    hs = gobjects(1, ns);
    for g = 1:numel(groups)
        for j = 1:ns
            x = g - 0.35 + w * (j - 0.5);
            if strcmp(spec{s, 1}, 'denom')
                v = B{g, j}.self - B{g, j}.ctrl;  ci = [NaN NaN];
            else
                v = B{g, j}.(spec{s, 1});  ci = B{g, j}.([spec{s, 1} '_ci']);
            end
            hs(j) = bar(ax, x, v, w * 0.9, 'FaceColor', colors(j, :), 'EdgeColor', 'k');
            if ~isnan(ci(1)), errorbar(ax, x, v, v - ci(1), ci(2) - v, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 4); end
        end
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3]);
    set(ax, 'XTick', 1:numel(groups), 'XTickLabel', groups, 'FontSize', 9);
    ylabel(ax, spec{s, 2});  box(ax, 'on');  xlim(ax, [0.5, numel(groups) + 0.5]);
    if s == 1
        lg = legend(ax, hs, series_names, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
    end
end
title(tl, ttl, 'FontSize', 11);
subtitle(tl, 'bars: pooled means; error bars: bootstrap 95% CI over cell draws (folds when a whole model is used once)', 'FontSize', 9);
save_figure(fig, S, out_dir, name);
end
