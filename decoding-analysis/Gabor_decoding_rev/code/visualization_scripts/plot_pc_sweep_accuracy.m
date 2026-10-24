%{
Filename: plot_pc_sweep_accuracy.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Gabor-model decoding accuracy vs the number of retained PCs per side k (log
    axis, k = 1-10, 15-45, 49; 48 for the complex model), whole models:
    self-decoding, PT (rotation only), corrected control and non-PT (source
    trials projected onto the target's PCs). Mean +- SEM over all repeat
    estimates; dotted line at k = 5. The full-space whole-model values are drawn
    as circles at "full".

Output (results/figures/06_pc_sweep_accuracy/):
    pc_gabor_accuracy_vs_k      one panel per model
    pc_sweep_accuracy.csv       every plotted point (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '06_pc_sweep_accuracy';
fprintf('plot_pc_sweep_accuracy\n');

F  = {'self', 'pt_rot', 'ctrl_m', 'nopt'};
LN = {'self-decoding', 'PT (rotation only)', 'PT control (corrected)', 'non-PT (target PCs)'};
LS = {'-', '--', ':', '-.'};  LW = [1.5 1.5 1.8 1];
k_full = 75;                                      % x position of the full-space markers
rows = {};

fig = figure('Visible', 'off', 'Position', [100 100 1400 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:2
    m = S.models{j};  r = R.pc.(m);  rf = R.cue.(m);  c = S.c.model.(m);
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(F));
    k = r.k_list(:)';
    iN = numel(rf.N_list);                        % whole model, full space
    for f = 1:numel(F)
        [mu, se] = deal(nan(size(k)));
        for ik = 1:numel(k)
            [mu(ik), se(ik)] = mean_sem(r.(F{f})(ik, :, :));
            rows(end + 1, :) = {S.label.(m), k(ik), LN{f}, mu(ik), se(ik)}; %#ok<SAGROW>
        end
        h(f) = errorbar(ax, k, mu, se, 'Color', c, 'LineStyle', LS{f}, 'LineWidth', LW(f), 'Marker', '.', 'MarkerSize', 10, 'CapSize', 2);
        [mf, sf] = mean_sem(rf.(F{f})(iN, :, 1, :));
        errorbar(ax, k_full, mf, sf, 'o', 'Color', c, 'MarkerSize', 4, 'HandleVisibility', 'off');
        rows(end + 1, :) = {S.label.(m), NaN, [LN{f} ', full space'], mf, sf}; %#ok<SAGROW>
    end
    yline(ax, S.chance, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, S.k_main, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xline(ax, 60, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');  xlim(ax, [0.9 95]);
    set(ax, 'XTick', [1 2 5 10 20 49 k_full], 'XTickLabel', {'1', '2', '5', '10', '20', '49', 'full'});
    ylim(ax, [0 1.05]);  box(ax, 'on');  xlabel(ax, 'Retained PCs per side, k');  ylabel(ax, 'Decoding accuracy');
    title(ax, sprintf('%s (%d cells)', S.label.(m), r.N));
    if j == 1, legend(ax, h, LN, 'Location', 'east'); end
end
title(tl, sprintf(['Gabor-model cue transfer vs number of PCs, whole models (mean +- SEM over repeat estimates; ' ...
    'dotted line: k = %d; circles at "full": full space, same cells)'], S.k_main), 'FontSize', 11);
save_figure(fig, S, out_dir, 'pc_gabor_accuracy_vs_k');
save_table(cell2table(rows, 'VariableNames', {'model', 'k', 'series', 'mean', 'sem'}), S, out_dir, 'pc_sweep_accuracy');
