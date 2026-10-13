%{
Filename: plot_normalized_metric_bars.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Normalized transfer metrics as bar charts, all from pooled means
    (transfer_metrics.m):
        PT / control,   headroom = (PT - control) / (self - control),   PT - control.
    Each figure has one panel per metric. Error bars are percentile-bootstrap 95%
    CIs over the resampling unit (neuron draws; folds when the whole population
    is used once).

    Note on headroom: its denominator is the TARGET's self-decoding. Pooling both
    directions of a comparison (as for the cue pairs) makes it symmetric. Per
    direction it can exceed 1 when the tested (source) population is more
    decodable than the target, e.g. KO V2-to-V1 across areas.

Output (results/figures/04_normalized_metric_bars/):
    cue_metrics_full_vs_k5        cue transfer, whole populations: full neuron space vs top 5 PCs
    cue_metrics_per_pair_full     per rendering pair, Fig 3b design (V2 subsampled to the V1 count)
    cue_metrics_per_pair_k5       per rendering pair, top 5 PCs, whole populations
    halves_metrics_full_vs_k5     halves of floor(N/2) neurons: full neuron space vs top 5 PCs
    across_metrics_full_vs_k5     across areas, per direction and both directions pooled
    cue_vs_halves_matched_size    cue transfer and halves at the same size (FR 24, KO 54)
    normalized_metric_bars.csv    every plotted number
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '04_normalized_metric_bars';
fprintf('plot_normalized_metric_bars\n');
rows = {};
seed = 100;
pops = S.pops;
popname = @(p) strrep(p, '_', ' ');

%% A. cue transfer: full space vs k = 5, whole populations
[Bf, Bk] = deal(cell(1, numel(pops)));
for ip = 1:numel(pops)
    r  = R.cue.(pops{ip});  iN = numel(r.N_list);          % whole population (one draw)
    [Bf{ip}, seed] = boot(r.self(iN, :, 1, :), r.pt_rot(iN, :, 1, :), r.ctrl_m(iN, :, 1, :), 4, S, seed);
    rows = add_row(rows, 'cue', 'full', popname(pops{ip}), 'all pairs', r.N_list(iN), Bf{ip});
    rk = R.pc_cue.(pops{ip});  ik = find(rk.k_list == S.k_main);
    [Bk{ip}, seed] = boot(rk.self(ik, :, :), rk.pt_rot(ik, :, :), rk.ctrl_m(ik, :, :), 3, S, seed);
    rows = add_row(rows, 'cue', sprintf('k=%d', S.k_main), popname(pops{ip}), 'all pairs', rk.N, Bk{ip});
end
two_space_figure(S, pops, Bf, Bk, 'Cue transfer, whole populations: full neuron space vs top 5 PCs per cue', ...
    out_dir, 'cue_metrics_full_vs_k5');

%% B. cue transfer per pair: Fig 3b design and k = 5
groups = {'EC-AC', 'EC-EX', 'AC-EX'};
for space = {'full', 'k5'}
    B = cell(numel(groups), numel(pops));
    for ip = 1:numel(pops)
        m = pops{ip}(1:2);
        for g = 1:numel(groups)
            if strcmp(space{1}, 'full')
                r  = R.cue.(pops{ip});  iN = find(r.N_list == S.anchor_N.(m));  D = r.n_draws(iN);
                P  = strcmp(r.pair_group, groups{g});
                ud = 3 + (D == 1);                               % draws, or folds when one draw
                [B{g, ip}, seed] = boot(r.self(iN, P, 1:D, :), r.pt_rot(iN, P, 1:D, :), r.ctrl_m(iN, P, 1:D, :), ud, S, seed);
                rows = add_row(rows, 'cue', 'full (Fig 3b size)', popname(pops{ip}), groups{g}, S.anchor_N.(m), B{g, ip});
            else
                r  = R.pc_cue.(pops{ip});  ik = find(r.k_list == S.k_main);
                P  = strcmp(r.pair_group, groups{g});
                [B{g, ip}, seed] = boot(r.self(ik, P, :), r.pt_rot(ik, P, :), r.ctrl_m(ik, P, :), 3, S, seed);
                rows = add_row(rows, 'cue', sprintf('k=%d', S.k_main), popname(pops{ip}), groups{g}, r.N, B{g, ip});
            end
        end
    end
    if strcmp(space{1}, 'full')
        ttl = 'Cue transfer per rendering pair, full neuron space (V1 whole, V2 subsampled to the V1 count)';
    else
        ttl = sprintf('Cue transfer per rendering pair, top %d PCs per cue (whole populations)', S.k_main);
    end
    per_group_figure(S, B, groups, cellfun(popname, pops, 'UniformOutput', false), ttl, ...
        out_dir, ['cue_metrics_per_pair_' space{1}]);
end

%% C. halves: largest halves, full space vs k = 5
[Bf, Bk] = deal(cell(1, numel(pops)));
for ip = 1:numel(pops)
    r  = R.halves.(pops{ip});  iH = numel(r.h_list);       % floor(N/2), as the PC sweep
    [Bf{ip}, seed] = boot(r.self(iH, :, :, :), r.pt_rot(iH, :, :, :), r.ctrl_m(iH, :, :, :), 3, S, seed);
    rows = add_row(rows, 'halves', 'full', popname(pops{ip}), 'all renderings', r.h_list(iH), Bf{ip});
    rk = R.pc_halves.(pops{ip});  ik = find(rk.k_list == S.k_main);
    [Bk{ip}, seed] = boot(rk.self(ik, :, :, :), rk.pt_rot(ik, :, :, :), rk.ctrl_m(ik, :, :, :), 3, S, seed);
    rows = add_row(rows, 'halves', sprintf('k=%d', S.k_main), popname(pops{ip}), 'all renderings', rk.half_size, Bk{ip});
end
two_space_figure(S, pops, Bf, Bk, 'Within-area halves of floor(N/2) neurons: full neuron space vs top 5 PCs per half', ...
    out_dir, 'halves_metrics_full_vs_k5');

%% D. across areas: per direction and pooled, full space vs k = 5
labels = {};  [Bf, Bk] = deal({});  colors = [];
sets  = {1, 2, [1 2]};  names = {'V1-to-V2', 'V2-to-V1', 'both directions'};
cols  = [S.c.V1; S.c.V2; 0.45 0.45 0.45];
for im = 1:2
    m  = S.monkeys{im};
    r  = R.across.(m);  rk = R.pc_across.(m);  ik = find(rk.k_list == S.k_main);
    for s = 1:3
        di = sets{s};
        % full space [dirs x combos x draws x folds]: unit = draws (dim 3)
        % PC sweep at one k [dirs x combos x 1 x folds]: unit = folds (dim 4)
        [Bf{end + 1}, seed] = boot(r.self(di, :, :, :), r.pt_rot(di, :, :, :), r.ctrl_m(di, :, :, :), 3, S, seed); %#ok<SAGROW>
        [Bk{end + 1}, seed] = boot(rk.self(di, :, ik, :), rk.pt_rot(di, :, ik, :), rk.ctrl_m(di, :, ik, :), 4, S, seed); %#ok<SAGROW>
        labels{end + 1} = sprintf('%s %s', m, names{s}); %#ok<SAGROW>
        colors(end + 1, :) = cols(s, :); %#ok<SAGROW>
        rows = add_row(rows, 'across', 'full', m, names{s}, r.N_match, Bf{end});
        rows = add_row(rows, 'across', sprintf('k=%d', S.k_main), m, names{s}, min(rk.N_area), Bk{end});
    end
end
two_space_figure(S, labels, Bf, Bk, ['Across areas: full neuron space (V2 subsampled to the V1 count) vs top 5 PCs ' ...
    '(whole populations); headroom per direction uses the target area''s self-decoding'], out_dir, ...
    'across_metrics_full_vs_k5', colors);

%% E. cue transfer vs halves at the same size
B = cell(2, numel(pops));
for ip = 1:numel(pops)
    m = pops{ip}(1:2);  n = S.anchor_h.(m);
    r = R.cue.(pops{ip});  iN = find(r.N_list == n);  D = r.n_draws(iN);
    [B{1, ip}, seed] = boot(r.self(iN, :, 1:D, :), r.pt_rot(iN, :, 1:D, :), r.ctrl_m(iN, :, 1:D, :), 3, S, seed);
    rows = add_row(rows, 'cue', 'full (matched to halves)', popname(pops{ip}), 'all pairs', n, B{1, ip});
    r = R.halves.(pops{ip});  iH = find(r.h_list == n);
    [B{2, ip}, seed] = boot(r.self(iH, :, :, :), r.pt_rot(iH, :, :, :), r.ctrl_m(iH, :, :, :), 3, S, seed);
    rows = add_row(rows, 'halves', 'full (matched to cue)', popname(pops{ip}), 'all renderings', n, B{2, ip});
end
per_group_figure(S, B, {'cue transfer', 'within-area halves'}, cellfun(popname, pops, 'UniformOutput', false), ...
    sprintf(['Same number of neurons (FR %d, KO %d), full neuron space: cue transfer vs two halves of one area ' ...
    '(the same-geometry reference)'], S.anchor_h.FR, S.anchor_h.KO), out_dir, 'cue_vs_halves_matched_size');

%% table
T = cell2table(rows, 'VariableNames', {'comparison', 'space', 'population', 'set', 'n_or_size', ...
    'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', 'headroom', 'headroom_lo', 'headroom_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'n_units'});
save_table(T, S, out_dir, 'normalized_metric_bars');

%% ===========================================================================
function [B, seed] = boot(Sa, Pa, Ca, unit_dim, S, seed)
% bootstrap_metrics with a fresh seed per call
seed = seed + 1;
B = bootstrap_metrics(Sa, Pa, Ca, unit_dim, S, seed);
end

function rows = add_row(rows, comparison, space, pop, set, n, B)
rows(end + 1, :) = {comparison, space, pop, set, n, B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
    B.headroom, B.headroom_ci(1), B.headroom_ci(2), B.diff, B.diff_ci(1), B.diff_ci(2), B.n_units};
end

function two_space_figure(S, labels, Bf, Bk, ttl, out_dir, name, colors)
% three panels (PT/control, headroom, PT - control); per group: full space (solid) and k = 5 (light, dashed edge)
if nargin < 8
    colors = zeros(numel(labels), 3);
    for i = 1:numel(labels)
        if contains(labels{i}, 'V1'), colors(i, :) = S.c.V1; else, colors(i, :) = S.c.V2; end
    end
end
labels = strrep(labels, '_', ' ');
fig = figure('Visible', 'off', 'Position', [100 100 1500 480], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {'ratio', 'PT / control', 1; 'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0; 'diff', 'PT - control (accuracy)', 0};
for s = 1:3
    ax = nexttile(tl);  hold(ax, 'on');
    for i = 1:numel(labels)
        for j = 1:2
            if j == 1, B = Bf{i}; else, B = Bk{i}; end
            x = i + (j - 1.5) * 0.36;
            v = B.(spec{s, 1});  ci = B.([spec{s, 1} '_ci']);
            b = bar(ax, x, v, 0.34, 'FaceColor', colors(i, :), 'EdgeColor', 'k');
            if j == 2, b.FaceAlpha = 0.35; b.LineStyle = '--'; end
            errorbar(ax, x, v, v - ci(1), ci(2) - v, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 4);
        end
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3]);
    set(ax, 'XTick', 1:numel(labels), 'XTickLabel', labels, 'XTickLabelRotation', 20, 'FontSize', 9);
    ylabel(ax, spec{s, 2});  box(ax, 'on');  xlim(ax, [0.4, numel(labels) + 0.6]);
    if s == 1
        hf = bar(ax, nan, nan, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'k');
        hk = bar(ax, nan, nan, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.35, 'EdgeColor', 'k', 'LineStyle', '--');
        lg = legend(ax, [hf hk], {'full neuron space', sprintf('top %d PCs', S.k_main)}, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
    end
end
title(tl, ttl, 'FontSize', 11);
subtitle(tl, 'bars: pooled means; error bars: bootstrap 95% CI over neuron draws (folds when a whole population is used once)', 'FontSize', 9);
save_figure(fig, S, out_dir, name);
end

function per_group_figure(S, B, groups, series_names, ttl, out_dir, name)
% B{group, series}; series = populations (KO solid, FR light), V1 blue / V2 orange
fig = figure('Visible', 'off', 'Position', [100 100 1500 480], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {'ratio', 'PT / control', 1; 'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0; 'diff', 'PT - control (accuracy)', 0};
ns = numel(series_names);  w = 0.8 / ns;
for s = 1:3
    ax = nexttile(tl);  hold(ax, 'on');
    hs = gobjects(1, ns);
    for g = 1:numel(groups)
        for j = 1:ns
            x  = g - 0.4 + w * (j - 0.5);
            v  = B{g, j}.(spec{s, 1});  ci = B{g, j}.([spec{s, 1} '_ci']);
            if contains(series_names{j}, 'V1'), c = S.c.V1; else, c = S.c.V2; end
            hs(j) = bar(ax, x, v, w * 0.9, 'FaceColor', c, 'EdgeColor', 'k');
            if startsWith(series_names{j}, 'FR'), hs(j).FaceAlpha = 0.4; end
            errorbar(ax, x, v, v - ci(1), ci(2) - v, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 3);
        end
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3]);
    set(ax, 'XTick', 1:numel(groups), 'XTickLabel', groups, 'FontSize', 10);
    ylabel(ax, spec{s, 2});  box(ax, 'on');  xlim(ax, [0.5, numel(groups) + 0.5]);
    if s == 1
        lg = legend(ax, hs, series_names, 'Orientation', 'horizontal');
        lg.Layout.Tile = 'south';
    end
end
title(tl, ttl, 'FontSize', 11);
subtitle(tl, 'bars: pooled means; error bars: bootstrap 95% CI over neuron draws (folds when a whole population is used once)', 'FontSize', 9);
save_figure(fig, S, out_dir, name);
end
