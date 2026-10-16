%{
Filename: plot_normalized_metric_bars.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Normalized transfer metrics of the CNN layers as bar charts, all from pooled
    means (transfer_metrics.m): PT / control, PT - control, headroom =
    (PT - control) / (self - control), and the headroom denominator
    self - control. The last panel matters for the CNNs: self-decoding is near
    1 and the corrected control is high, so the denominator can be small and
    headroom then swings with small changes of PT or control.
    Error bars are percentile-bootstrap 95% CIs over the resampling unit: unit
    draws, or folds when a whole layer is used once (and for the hierarchy,
    whose 5 draws differ by one unit); the denominator panel has none.

Output (results/figures/04_normalized_metric_bars/):
    cue_metrics_full_vs_k5           whole layers: full unit space (solid) vs top 5 PCs (light)
    cue_metrics_fig3d_size_per_pair  per rendering pair at the Fig 3d size (96 / 64 units)
    hierarchy_metrics_full_vs_k5     AlexNet pool2 <-> pool5, per direction and both pooled
    normalized_metric_bars.csv       every plotted number
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '04_normalized_metric_bars';
fprintf('plot_normalized_metric_bars\n');
rows = {};
seed = 300;
labels = cellfun(@(k) R.cue.(k).population, S.keys, 'UniformOutput', false);
colors = cell2mat(cellfun(@(k) S.c.layer.(k), S.keys(:), 'UniformOutput', false));

%% A. whole layers: full space vs k = 5 (same units)
[Bf, Bk] = deal(cell(1, numel(S.keys)));
for i = 1:numel(S.keys)
    r  = R.cue.(S.keys{i});  iN = numel(r.N_list);                  % whole layer (one draw: folds)
    [Bf{i}, seed] = boot(r.self(iN, :, 1, :), r.pt_rot(iN, :, 1, :), r.ctrl_m(iN, :, 1, :), 4, S, seed);
    rows = add_row(rows, 'cue', 'full', labels{i}, 'all pairs', r.N_list(iN), Bf{i});
    rk = R.pc_cue.(S.keys{i});  ik = find(rk.k_list == S.k_main);
    [Bk{i}, seed] = boot(rk.self(ik, :, :), rk.pt_rot(ik, :, :), rk.ctrl_m(ik, :, :), 3, S, seed);
    rows = add_row(rows, 'cue', sprintf('k=%d', S.k_main), labels{i}, 'all pairs', rk.N, Bk{i});
end
two_space_figure(S, labels, Bf, Bk, colors, ...
    'CNN cue transfer, whole layers: full unit space vs top 5 PCs per rendering (same units)', out_dir, 'cue_metrics_full_vs_k5');

%% B. per rendering pair at the Fig 3d size
groups = {'EC-AC', 'EC-EX', 'AC-EX'};
B = cell(numel(groups), numel(S.keys));
for i = 1:numel(S.keys)
    r  = R.cue.(S.keys{i});
    iN = find(r.N_list == S.fig3d_N.(r.net));  D = r.n_draws(iN);
    ud = 3 + (D == 1);                                              % draws, or folds when one draw
    for g = 1:numel(groups)
        P = strcmp(r.pair_group, groups{g});
        [B{g, i}, seed] = boot(r.self(iN, P, 1:D, :), r.pt_rot(iN, P, 1:D, :), r.ctrl_m(iN, P, 1:D, :), ud, S, seed);
        rows = add_row(rows, 'cue', 'full (Fig 3d size)', labels{i}, groups{g}, r.N_list(iN), B{g, i});
    end
end
per_group_figure(S, B, groups, labels, colors, ...
    'CNN cue transfer per rendering pair, full unit space at the Fig 3d size (AlexNet 96, VGG 64 units)', ...
    out_dir, 'cue_metrics_fig3d_size_per_pair');

%% C. hierarchy: per direction and pooled, full space vs k = 5
r = R.hier;  rk = R.pc_hier;  ik = find(rk.k_list == S.k_main);
sets  = {1, 2, [1 2]};
names = {'pool2-to-pool5', 'pool5-to-pool2', 'both directions'};
[Bf, Bk] = deal(cell(1, 3));
for s = 1:3
    di = sets{s};
    [Bf{s}, seed] = boot(r.self(di, :, :, :), r.pt_rot(di, :, :, :), r.ctrl_m(di, :, :, :), 4, S, seed);   % folds
    [Bk{s}, seed] = boot(rk.self(di, :, ik, :), rk.pt_rot(di, :, ik, :), rk.ctrl_m(di, :, ik, :), 4, S, seed);
    rows = add_row(rows, 'hierarchy', 'full', 'AlexNet pool2 <-> pool5', names{s}, r.N_match, Bf{s});
    rows = add_row(rows, 'hierarchy', sprintf('k=%d', S.k_main), 'AlexNet pool2 <-> pool5', names{s}, min(rk.N_layer), Bk{s});
end
two_space_figure(S, names, Bf, Bk, [S.c.dir; 0.45 0.45 0.45], ...
    ['AlexNet pool2 ' char(8596) ' pool5: full unit space (pool2 subsampled to 249) vs top 5 PCs (whole layers)'], ...
    out_dir, 'hierarchy_metrics_full_vs_k5');

%% table
T = cell2table(rows, 'VariableNames', {'comparison', 'space', 'layer', 'set', 'n_units', ...
    'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', 'headroom', 'headroom_lo', 'headroom_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'self_minus_ctrl', 'n_units_boot'});
save_table(T, S, out_dir, 'normalized_metric_bars');

%% ===========================================================================
function [B, seed] = boot(Sa, Pa, Ca, unit_dim, S, seed)
% bootstrap_metrics with a fresh seed per call
seed = seed + 1;
B = bootstrap_metrics(Sa, Pa, Ca, unit_dim, S, seed);
end

function rows = add_row(rows, comparison, space, layer, set, n, B)
rows(end + 1, :) = {comparison, space, layer, set, n, B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
    B.headroom, B.headroom_ci(1), B.headroom_ci(2), B.diff, B.diff_ci(1), B.diff_ci(2), B.self - B.ctrl, B.n_units};
end

function spec = metric_spec()
% field, axis label, reference line
spec = {'ratio', 'PT / control', 1; 'diff', 'PT - control (accuracy)', 0; ...
        'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0; 'denom', 'self - control (headroom denominator)', 0};
end

function [v, lo, hi] = metric_value(B, f)
if strcmp(f, 'denom')
    v = B.self - B.ctrl;  lo = NaN;  hi = NaN;
else
    v = B.(f);  ci = B.([f '_ci']);  lo = ci(1);  hi = ci(2);
end
end

function two_space_figure(S, labels, Bf, Bk, colors, ttl, out_dir, name)
% four panels; per group: full space (solid) and k = 5 (light, dashed edge)
spec = metric_spec();
fig = figure('Visible', 'off', 'Position', [100 100 1900 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:size(spec, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    for i = 1:numel(labels)
        for j = 1:2
            if j == 1, B = Bf{i}; else, B = Bk{i}; end
            x = i + (j - 1.5) * 0.36;
            [v, lo, hi] = metric_value(B, spec{s, 1});
            b = bar(ax, x, v, 0.34, 'FaceColor', colors(i, :), 'EdgeColor', 'k');
            if j == 2, b.FaceAlpha = 0.35; b.LineStyle = '--'; end
            if ~isnan(lo), errorbar(ax, x, v, v - lo, hi - v, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 3); end
        end
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3]);
    set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 30, 'FontSize', 9);
    ylabel(ax, spec{s, 2});  box(ax, 'on');  xlim(ax, [0.4, numel(labels) + 0.6]);
    if s == 1
        hf = bar(ax, nan, nan, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
        hk = bar(ax, nan, nan, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.35, 'EdgeColor', 'k', 'LineStyle', '--');
        lg = legend(ax, [hf hk], {'full unit space', sprintf('top %d PCs', S.k_main)}, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
    end
end
title(tl, ttl, 'FontSize', 11);
subtitle(tl, 'bars: pooled means; error bars: bootstrap 95% CI over folds (whole layers) or unit draws', 'FontSize', 9);
save_figure(fig, S, out_dir, name);
end

function per_group_figure(S, B, groups, series_names, colors, ttl, out_dir, name)
% B{group, series}; series = layers (AlexNet solid, VGG light with dashed edge)
spec = metric_spec();
fig = figure('Visible', 'off', 'Position', [100 100 1900 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
ns = numel(series_names);  w = 0.84 / ns;
for s = 1:size(spec, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    hs = gobjects(1, ns);
    for g = 1:numel(groups)
        for j = 1:ns
            x = g - 0.42 + w * (j - 0.5);
            [v, lo, hi] = metric_value(B{g, j}, spec{s, 1});
            hs(j) = bar(ax, x, v, w * 0.9, 'FaceColor', colors(j, :), 'EdgeColor', 'k');
            if startsWith(series_names{j}, 'VGG'), hs(j).FaceAlpha = 0.45; hs(j).LineStyle = '--'; end
            if ~isnan(lo), errorbar(ax, x, v, v - lo, hi - v, 'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 2); end
        end
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3]);
    set(ax, 'XTick', 1:numel(groups), 'XTickLabel', groups, 'FontSize', 10);
    ylabel(ax, spec{s, 2});  box(ax, 'on');  xlim(ax, [0.5, numel(groups) + 0.5]);
    if s == 1
        lg = legend(ax, hs, series_names, 'Orientation', 'horizontal', 'NumColumns', 8);
        lg.Layout.Tile = 'south';
    end
end
title(tl, ttl, 'FontSize', 11);
subtitle(tl, 'bars: pooled means; error bars: bootstrap 95% CI over unit draws (folds for the pool1 layers, used whole)', 'FontSize', 9);
save_figure(fig, S, out_dir, name);
end
