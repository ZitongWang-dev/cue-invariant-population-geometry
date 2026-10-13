%{
Filename: plot_sanity_old_vs_new.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Sanity check of the rerun against the published pipeline. Every bar of the
    old Figs 3b (cue transfer), 4a (halves) and 4b/12 (across areas, all nine
    combinations) is plotted as old value (x) vs new value (y), for self-decoding,
    non-PT and PT (rotation only), which should agree up to sampling noise, and
    for the control, where the old control (chance by construction) is set
    against the corrected one.

    Old values: neuronal_decoding/results/decoding_outputs/
        procrustes_decoding_basic_results                  (cols 1 self, 2 non-PT, 6 PT, 4 control)
        procrustes_decoding_half_and_half_populations_results
        procrustes_decoding_across_area_results
    New values: the same design points of the rerun (cue N = 48 / 109, halves
    24 / 54, across areas matched to the V1 count), one point per ordered pair,
    rendering or combination.

Output (results/figures/01_sanity_old_vs_new/):
    old_vs_new_all_bars.png/.fig   three panels: cue transfer, halves, across areas
    old_vs_new_bars.csv            every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '01_sanity_old_vs_new';
fprintf('plot_sanity_old_vs_new\n');

%% collect old and new values of every bar
rows = {};   % analysis, population/direction, item, metric, old, new
for im = 1:2
    m = S.monkeys{im};
    for a = 1:2
        pop = sprintf('%s_V%d', m, a);
        % ---- cue transfer, Fig 3b point ----
        r  = R.cue.(pop);
        iN = find(r.N_list == S.anchor_N.(m));
        for p = 1:numel(r.pair_labels)
            stem = strrep(r.pair_labels{p}, '->', '');
            old = load_old(fullfile(S.old_root, 'procrustes_decoding_basic_results', m, sprintf('V%d', a), [stem '_results.mat']));
            new = new_means(r, iN, p, r.n_draws(iN));
            rows = add_rows(rows, 'cue transfer', pop, r.pair_labels{p}, old, new);
        end
        % ---- halves, Fig 4a point ----
        r  = R.halves.(pop);
        iH = find(r.h_list == S.anchor_h.(m));
        for ri = 1:numel(r.renderings)
            old = load_old(fullfile(S.old_root, 'procrustes_decoding_half_and_half_populations_results', m, ...
                sprintf('V%d', a), [r.renderings{ri} '_results.mat']));
            new = new_means(r, iH, ri, size(r.self, 3));
            rows = add_rows(rows, 'halves', pop, upper(r.renderings{ri}), old, new);
        end
    end
    % ---- across areas, Figs 4b / 12 ----
    r = R.across.(m);
    old_dir = {'V1 decoding V2', 'V2 decoding V1'};     % old folder of direction 1 (train V1) and 2 (train V2)
    src_a   = {'v2', 'v1'};  tgt_a = {'v1', 'v2'};
    for di = 1:2
        for c = 1:numel(r.combo_source_cue)
            stem = sprintf('%s%s%s%s', src_a{di}, r.combo_source_cue{c}, tgt_a{di}, r.combo_target_cue{c});
            old = load_old(fullfile(S.old_root, 'procrustes_decoding_across_area_results', m, old_dir{di}, [stem '_results.mat']));
            new = [mean(r.self(di, c, :, :), 'all'), mean(r.nopt(di, c, :, :), 'all'), ...
                   mean(r.pt_rot(di, c, :, :), 'all'), mean(r.ctrl_m(di, c, :, :), 'all')];
            rows = add_rows(rows, 'across areas', sprintf('%s %s', m, r.dir_labels{di}(1:8)), ...
                sprintf('train %s, test %s', upper(r.combo_target_cue{c}), upper(r.combo_source_cue{c})), old, new);
        end
    end
end
T = cell2table(rows, 'VariableNames', {'analysis', 'population', 'item', 'metric', 'old', 'new'});
save_table(T, S, out_dir, 'old_vs_new_bars');

%% figure: one panel per analysis
analyses = {'cue transfer', 'halves', 'across areas'};
metrics  = {'self', 'non-PT', 'PT (rotation only)', 'control'};
mk       = {'o', '^', 's', 'x'};
col      = [0.25 0.25 0.25; 0.30 0.60 0.85; 0.85 0.33 0.10; 0.80 0.10 0.10];
fig = figure('Visible', 'off', 'Position', [100 100 1500 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for ia = 1:numel(analyses)
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(metrics));
    for k = 1:numel(metrics)
        sel = strcmp(T.analysis, analyses{ia}) & strcmp(T.metric, metrics{k});
        h(k) = scatter(ax, T.old(sel), T.new(sel), 36, col(k, :), mk{k}, 'LineWidth', 1.2);
    end
    lim = [0, ceil(max([T.old; T.new]) * 10) / 10];
    plot(ax, lim, lim, 'k-', 'HandleVisibility', 'off');
    yline(ax, S.chance, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xline(ax, S.chance, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    axis(ax, 'square');  xlim(ax, lim);  ylim(ax, lim);  box(ax, 'on');  grid(ax, 'on');
    xlabel(ax, 'Old pipeline (published figure)');  ylabel(ax, 'Rerun (neuron\_decoding\_rev)');
    % agreement of the three quantities that should not have changed
    sel = strcmp(T.analysis, analyses{ia}) & ~strcmp(T.metric, 'control');
    dv  = T.new(sel) - T.old(sel);
    title(ax, sprintf('%s: self / non-PT / PT  mean |new - old| = %.3f, max %.3f', ...
        analyses{ia}, mean(abs(dv)), max(abs(dv))), 'FontWeight', 'normal', 'FontSize', 9);
    if ia == 1
        legend(ax, h, {'self-decoding', 'non-PT', 'PT (rotation only)', ['control: old ' char(8594) ' corrected']}, 'Location', 'northwest');
    end
end
title(tl, 'Rerun vs published pipeline, every bar of Figs 3b, 4a, 4b/12 (both monkeys, both areas)');
save_figure(fig, S, out_dir, 'old_vs_new_all_bars');

%% ===========================================================================
function v = load_old(f)
% [self, non-PT, PT rotation-only, old control] means of an old result file
L = load(f);  fn = fieldnames(L);  X = L.(fn{1});
if iscell(X), X = X{1}; end
v = mean(X(:, [1 2 6 4]), 1);
end

function v = new_means(r, i1, i2, nd)
% [self, non-PT, PT rotation-only, corrected control] means of one rerun design point
v = [mean(r.self(i1, i2, 1:nd, :), 'all', 'omitnan'), mean(r.nopt(i1, i2, 1:nd, :), 'all', 'omitnan'), ...
     mean(r.pt_rot(i1, i2, 1:nd, :), 'all', 'omitnan'), mean(r.ctrl_m(i1, i2, 1:nd, :), 'all', 'omitnan')];
end

function rows = add_rows(rows, analysis, pop, item, old, new)
names = {'self', 'non-PT', 'PT (rotation only)', 'control'};
for k = 1:4
    rows(end + 1, :) = {analysis, pop, item, names{k}, old(k), new(k)}; %#ok<AGROW>
end
end
