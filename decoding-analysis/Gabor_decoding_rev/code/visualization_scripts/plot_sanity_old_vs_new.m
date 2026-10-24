%{
Filename: plot_sanity_old_vs_new.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Sanity check of the Gabor-model rerun against the published pipeline
    (Gabor_decoding/). Old value (x) vs new value (y) for self-decoding, non-PT
    and PT (rotation only), which should agree up to sampling noise, and for the
    control, where the old control (chance by construction) is set against the
    corrected one.
      Fig 15         simple at 48 cells (published "downsampled"), simple at 192
                     and complex at 48, every ordered rendering pair
      point 9 sweep  simple model at 10, 20, ..., 190, 192 cells, every pair
                     (procrustes_Gabor_resps_decoding_basic_incremental_pr.m)
    Old result columns: 1 self, 2 non-PT, 6 PT rotation-only, 4 old control.

Output (results/figures/01_sanity_old_vs_new/):
    old_vs_new_all_points.png/.fig   two panels: Fig 15, point 9 sweep
    old_vs_new_points.csv            every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '01_sanity_old_vs_new';
fprintf('plot_sanity_old_vs_new\n');

rows = {};   % analysis, model, item, metric, old, new
% ---- Fig 15: published basic runs ----
designs = {'simple', 48, 'even and odd combined downsampled'; 'simple', 192, 'even and odd combined'; 'complex', 48, 'complex'};
for i = 1:size(designs, 1)
    r  = R.cue.(designs{i, 1});
    iN = find(r.N_list == designs{i, 2});
    for p = 1:numel(r.pair_labels)
        stem = strrep(r.pair_labels{p}, '->', '');
        L = load(fullfile(S.old_root, designs{i, 3}, [stem '_results.mat']));  fn = fieldnames(L);
        old = mean(L.(fn{1})(:, [1 2 6 4]), 1);
        rows = add_rows(rows, 'Fig 15', sprintf('%s %d', S.label.(designs{i, 1}), designs{i, 2}), r.pair_labels{p}, ...
            old, new_means(r, iN, p));
    end
end
% ---- rebuttal point 9: incremental sweep of the simple model ----
r = R.cue.simple;
for p = 1:numel(r.pair_labels)
    stem = strrep(r.pair_labels{p}, '->', '');
    L = load(fullfile(S.old_root, 'procrustes_Gabor_incremental_pr_results', 'even and odd combined', [stem '_results.mat']));
    fn = fieldnames(L);  C = L.(fn{1});
    for j = 1:numel(C)
        iN = find(r.N_list == C{j}.neuron_num);
        if isempty(iN), continue; end
        rows = add_rows(rows, 'point 9 sweep', S.label.simple, sprintf('%s N=%d', r.pair_labels{p}, C{j}.neuron_num), ...
            mean(C{j}.acc(:, [1 2 6 4]), 1), new_means(r, iN, p));
    end
end
T = cell2table(rows, 'VariableNames', {'analysis', 'model', 'item', 'metric', 'old', 'new'});
save_table(T, S, out_dir, 'old_vs_new_points');

%% figure: one panel per analysis
analyses = {'Fig 15', 'point 9 sweep'};
titles   = {'Fig 15 (simple 48 and 192, complex 48)', 'Rebuttal point 9 sweep (simple, 10-192 cells)'};
metrics  = {'self', 'non-PT', 'PT (rotation only)', 'control'};
mk       = {'o', '^', 's', 'x'};
col      = [0.25 0.25 0.25; 0.30 0.60 0.85; 0.85 0.33 0.10; 0.80 0.10 0.10];
fig = figure('Visible', 'off', 'Position', [100 100 1050 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ia = 1:numel(analyses)
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(metrics));
    for k = 1:numel(metrics)
        sel = strcmp(T.analysis, analyses{ia}) & strcmp(T.metric, metrics{k});
        h(k) = scatter(ax, T.old(sel), T.new(sel), 24, col(k, :), mk{k}, 'LineWidth', 1);
    end
    lim = [0 1.02];
    plot(ax, lim, lim, 'k-', 'HandleVisibility', 'off');
    yline(ax, S.chance, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xline(ax, S.chance, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    axis(ax, 'square');  xlim(ax, lim);  ylim(ax, lim);  box(ax, 'on');  grid(ax, 'on');
    xlabel(ax, 'Old pipeline (Gabor\_decoding)');  ylabel(ax, 'Rerun (Gabor\_decoding\_rev)');
    sel = strcmp(T.analysis, analyses{ia}) & ~strcmp(T.metric, 'control');
    dv  = T.new(sel) - T.old(sel);
    title(ax, {titles{ia}, sprintf('self / non-PT / PT: mean |new - old| = %.3f, max %.3f', mean(abs(dv)), max(abs(dv)))}, ...
        'FontWeight', 'normal', 'FontSize', 9);
    if ia == 1
        legend(ax, h, {'self-decoding', 'non-PT', 'PT (rotation only)', ['control: old ' char(8594) ' corrected']}, 'Location', 'northwest');
    end
end
title(tl, 'Gabor-model rerun vs published pipeline (6 ordered rendering pairs)');
save_figure(fig, S, out_dir, 'old_vs_new_all_points');

%% ===========================================================================
function v = new_means(r, iN, p)
% [self, non-PT, PT rotation-only, corrected control] means of one rerun design point
D = r.n_draws(iN);
v = [mean(r.self(iN, p, 1:D, :), 'all'), mean(r.nopt(iN, p, 1:D, :), 'all'), ...
     mean(r.pt_rot(iN, p, 1:D, :), 'all'), mean(r.ctrl_m(iN, p, 1:D, :), 'all')];
end

function rows = add_rows(rows, analysis, model, item, old, new)
names = {'self', 'non-PT', 'PT (rotation only)', 'control'};
for k = 1:4
    rows(end + 1, :) = {analysis, model, item, names{k}, old(k), new(k)}; %#ok<AGROW>
end
end
