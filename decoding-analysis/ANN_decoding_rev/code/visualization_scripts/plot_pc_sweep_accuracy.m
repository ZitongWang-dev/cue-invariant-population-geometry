%{
Filename: plot_pc_sweep_accuracy.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    CNN decoding accuracy vs the number of retained PCs per side k (log axis,
    k = 1-10, 15-45, 49), whole layers: self-decoding, PT (rotation only),
    corrected control and non-PT (source trials projected onto the target's
    PCs; cue transfer only). Mean +- SEM over all repeat estimates; dotted line
    at k = 5. The full-space whole-layer values are drawn as markers at "full".

Output (results/figures/06_pc_sweep_accuracy/):
    pc_cue_accuracy_vs_k         one panel per layer (AlexNet top row, VGG-19 bottom row)
    pc_hierarchy_accuracy_vs_k   AlexNet pool2 <-> pool5, one panel per direction
    pc_sweep_accuracy.csv        every plotted point (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '06_pc_sweep_accuracy';
fprintf('plot_pc_sweep_accuracy\n');

F  = {'self', 'pt_rot', 'ctrl_m', 'nopt'};
LN = {'self-decoding', 'PT (rotation only)', 'PT control (corrected)', 'non-PT (target PCs)'};
LS = {'-', '--', ':', '-.'};  LW = [1.5 1.5 1.8 1];
k_full = 75;                                      % x position of the full-space markers
rows = {};

%% cue transfer, one panel per layer
fig = figure('Visible', 'off', 'Position', [100 100 2000 800], 'Color', 'w');
tl  = tiledlayout(fig, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
tiles = [1 2 3 6 7 8 9 10];                       % AlexNet in row 1, VGG-19 in row 2
for i = 1:numel(S.keys)
    key = S.keys{i};  r = R.pc_cue.(key);  rf = R.cue.(key);  c = S.c.layer.(key);
    ax = nexttile(tl, tiles(i));  hold(ax, 'on');
    h = gobjects(1, numel(F));
    k = r.k_list(:)';
    iN = numel(rf.N_list);                        % whole layer, full space
    for f = 1:numel(F)
        [mu, se] = deal(nan(size(k)));
        for ik = 1:numel(k)
            [mu(ik), se(ik)] = mean_sem(r.(F{f})(ik, :, :));
            rows(end + 1, :) = {'cue', r.population, k(ik), LN{f}, mu(ik), se(ik)}; %#ok<SAGROW>
        end
        h(f) = errorbar(ax, k, mu, se, 'Color', c, 'LineStyle', LS{f}, 'LineWidth', LW(f), 'Marker', '.', 'MarkerSize', 10, 'CapSize', 2);
        [mf, sf] = mean_sem(rf.(F{f})(iN, :, 1, :));
        errorbar(ax, k_full, mf, sf, 'o', 'Color', c, 'MarkerSize', 4, 'HandleVisibility', 'off');
        rows(end + 1, :) = {'cue', r.population, NaN, [LN{f} ', full space'], mf, sf}; %#ok<SAGROW>
    end
    finish_k_axis(ax, S, k_full);
    ylim(ax, [0 1.05]);  ylabel(ax, 'Decoding accuracy');
    title(ax, sprintf('%s (%d units)', r.population, r.N));
    if i == 1, lg = legend(ax, h, LN); lg.Layout.Tile = 4; end
end
title(tl, sprintf(['CNN cue transfer vs number of PCs, whole layers (mean +- SEM over repeat estimates; dotted line: k = %d; ' ...
    'circles at "full": full unit space, same units)'], S.k_main), 'FontSize', 11);
save_figure(fig, S, out_dir, 'pc_cue_accuracy_vs_k');

%% hierarchy, one panel per direction
r = R.pc_hier;  rf = R.hier;
fig = figure('Visible', 'off', 'Position', [100 100 1400 480], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for di = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    c = S.c.dir(di, :);  k = r.k_list(:)';
    h = gobjects(1, 3);
    for f = 1:3
        [mu, se] = deal(nan(size(k)));
        for ik = 1:numel(k)
            [mu(ik), se(ik)] = mean_sem(r.(F{f})(di, :, ik, :));
            rows(end + 1, :) = {'hierarchy', r.dir_labels{di}, k(ik), LN{f}, mu(ik), se(ik)}; %#ok<SAGROW>
        end
        h(f) = errorbar(ax, k, mu, se, 'Color', c, 'LineStyle', LS{f}, 'LineWidth', LW(f), 'Marker', '.', 'MarkerSize', 10, 'CapSize', 2);
        [mf, sf] = mean_sem(rf.(F{f})(di, :, :, :));
        errorbar(ax, k_full, mf, sf, 'o', 'Color', c, 'MarkerSize', 4, 'HandleVisibility', 'off');
        rows(end + 1, :) = {'hierarchy', r.dir_labels{di}, NaN, [LN{f} ', full space'], mf, sf}; %#ok<SAGROW>
    end
    finish_k_axis(ax, S, k_full);
    ylim(ax, [0 1.05]);  ylabel(ax, 'Decoding accuracy');
    title(ax, r.dir_labels{di});
    if di == 1, legend(ax, h, LN(1:3), 'Location', 'northwest'); end
end
title(tl, sprintf('AlexNet pool2 %s pool5 vs number of PCs (whole layers: pool2 %d, pool5 %d units; circles at "full": full space, pool2 subsampled to %d)', ...
    char(8596), r.N_layer(1), r.N_layer(2), rf.N_match), 'FontSize', 11);
save_figure(fig, S, out_dir, 'pc_hierarchy_accuracy_vs_k');

save_table(cell2table(rows, 'VariableNames', {'comparison', 'panel', 'k', 'series', 'mean', 'sem'}), S, out_dir, 'pc_sweep_accuracy');

%% ===========================================================================
function finish_k_axis(ax, S, k_full)
yline(ax, S.chance, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xline(ax, S.k_main, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xline(ax, 60, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
set(ax, 'XScale', 'log');  xlim(ax, [0.9 95]);
set(ax, 'XTick', [1 2 5 10 20 49 k_full], 'XTickLabel', {'1', '2', '5', '10', '20', '49', 'full'});
box(ax, 'on');  xlabel(ax, 'Retained PCs per side, k');
end
