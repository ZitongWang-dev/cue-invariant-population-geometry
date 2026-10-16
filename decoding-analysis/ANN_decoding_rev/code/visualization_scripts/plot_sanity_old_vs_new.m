%{
Filename: plot_sanity_old_vs_new.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Sanity check of the CNN rerun against the published pipeline
    (ANN_decoding/). Old value (x) vs new value (y) for self-decoding, non-PT
    and PT (rotation only), which should agree up to sampling noise, and for the
    control, where the old control (chance by construction) is set against the
    corrected one.
      Fig 3d          every layer x ordered pair at 96 (AlexNet) / 64 (VGG) units
      Fig 4d          every layer x ordered pair at the unit counts both grids share
                      (10, 20, 40, 80, 160, 320 and the whole layer)
      hierarchy       AlexNet pool2 <-> pool5, both directions x 9 combinations
    Old result columns: 1 self, 2 non-PT, 6 PT rotation-only, 4 old control.

Output (results/figures/01_sanity_old_vs_new/):
    old_vs_new_all_points.png/.fig   three panels: Fig 3d, Fig 4d, hierarchy
    old_vs_new_points.csv            every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '01_sanity_old_vs_new';
fprintf('plot_sanity_old_vs_new\n');

rows = {};   % analysis, layer / direction, item, metric, old, new
for i = 1:numel(S.keys)
    key = S.keys{i};  net = S.layers{i, 1};  L = S.layers{i, 2};
    r = R.cue.(key);
    for p = 1:numel(r.pair_labels)
        stem = strrep(r.pair_labels{p}, '->', '');
        % ---- Fig 3d: basic run at the published size ----
        old = load_old(fullfile(S.old_root, 'procrustes_ann_decoding_basic_results', net, sprintf('pool%d', L), [stem '_results.mat']), 1);
        iN  = find(r.N_list == S.fig3d_N.(net));
        rows = add_rows(rows, 'Fig 3d', r.population, r.pair_labels{p}, old, new_means(r, iN, p));
        % ---- Fig 4d: incremental run at the shared unit counts ----
        f = fullfile(S.old_root, 'procrustes_ann_decoding_basic_incremental_results', net, sprintf('pool%d', L), [stem '_results.mat']);
        old_N = old_grid(r.N_total);
        for iN = 1:numel(r.N_list)
            io = find(old_N == r.N_list(iN));
            if isempty(io), continue; end
            rows = add_rows(rows, 'Fig 4d', r.population, sprintf('%s N=%d', r.pair_labels{p}, r.N_list(iN)), ...
                load_old(f, io), new_means(r, iN, p));
        end
    end
end
% ---- hierarchy ----
r = R.hier;
old_dir = {'V1 decoding V2', 'V2 decoding V1'};    % old folder of direction 1 (train pool2) and 2 (train pool5)
src_a   = {'v2', 'v1'};  tgt_a = {'v1', 'v2'};     % old names: v1 = lower layer (pool2), v2 = higher (pool5)
for di = 1:2
    for c = 1:numel(r.combo_source_cue)
        stem = sprintf('%s%s%s%s', src_a{di}, r.combo_source_cue{c}, tgt_a{di}, r.combo_target_cue{c});
        old = load_old(fullfile(S.old_root, 'procrustes_ann_decoding_along_hierarchy_results', 'alexnet', ...
            'alexnet POOL 2&5', old_dir{di}, [stem '_results.mat']), 1);
        new = [mean(r.self(di, c, :, :), 'all'), mean(r.nopt(di, c, :, :), 'all'), ...
               mean(r.pt_rot(di, c, :, :), 'all'), mean(r.ctrl_m(di, c, :, :), 'all')];
        rows = add_rows(rows, 'hierarchy', r.dir_labels{di}(1:14), ...
            sprintf('train %s, test %s', upper(r.combo_target_cue{c}), upper(r.combo_source_cue{c})), old, new);
    end
end
T = cell2table(rows, 'VariableNames', {'analysis', 'layer', 'item', 'metric', 'old', 'new'});
save_table(T, S, out_dir, 'old_vs_new_points');

%% figure: one panel per analysis
analyses = {'Fig 3d', 'Fig 4d', 'hierarchy'};
titles   = {'Fig 3d (96 / 64 units)', 'Fig 4d (shared unit counts)', ['AlexNet pool2 ' char(8596) ' pool5']};
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
        h(k) = scatter(ax, T.old(sel), T.new(sel), 24, col(k, :), mk{k}, 'LineWidth', 1);
    end
    lim = [0 1.02];
    plot(ax, lim, lim, 'k-', 'HandleVisibility', 'off');
    yline(ax, S.chance, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xline(ax, S.chance, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    axis(ax, 'square');  xlim(ax, lim);  ylim(ax, lim);  box(ax, 'on');  grid(ax, 'on');
    xlabel(ax, 'Old pipeline (ANN\_decoding)');  ylabel(ax, 'Rerun (ANN\_decoding\_rev)');
    sel = strcmp(T.analysis, analyses{ia}) & ~strcmp(T.metric, 'control');
    dv  = T.new(sel) - T.old(sel);
    title(ax, {titles{ia}, sprintf('self / non-PT / PT: mean |new - old| = %.3f, max %.3f', mean(abs(dv)), max(abs(dv)))}, ...
        'FontWeight', 'normal', 'FontSize', 9);
    if ia == 1
        legend(ax, h, {'self-decoding', 'non-PT', 'PT (rotation only)', ['control: old ' char(8594) ' corrected']}, 'Location', 'northwest');
    end
end
title(tl, 'CNN rerun vs published pipeline (all 8 layers, 6 ordered rendering pairs)');
save_figure(fig, S, out_dir, 'old_vs_new_all_points');

%% ===========================================================================
function v = load_old(f, io)
% [self, non-PT, PT rotation-only, old control] means of an old result file
% (io = position in the old unit grid when the file holds a cell per size)
L = load(f);  fn = fieldnames(L);  X = L.(fn{1});
if iscell(X), X = X{io}; end
v = mean(X(:, [1 2 6 4]), 1);
end

function N = old_grid(N_total)
% unit grid of procrustes_ann_decoding_basic_incremental.m: 5 * 2^n below the layer size, plus the layer
N = 5 * 2.^(0:10);
N = unique([N(N < N_total), N_total]);
end

function v = new_means(r, iN, p)
% [self, non-PT, PT rotation-only, corrected control] means of one rerun design point
D = r.n_draws(iN);
v = [mean(r.self(iN, p, 1:D, :), 'all'), mean(r.nopt(iN, p, 1:D, :), 'all'), ...
     mean(r.pt_rot(iN, p, 1:D, :), 'all'), mean(r.ctrl_m(iN, p, 1:D, :), 'all')];
end

function rows = add_rows(rows, analysis, layer, item, old, new)
names = {'self', 'non-PT', 'PT (rotation only)', 'control'};
for k = 1:4
    rows(end + 1, :) = {analysis, layer, item, names{k}, old(k), new(k)}; %#ok<AGROW>
end
end
