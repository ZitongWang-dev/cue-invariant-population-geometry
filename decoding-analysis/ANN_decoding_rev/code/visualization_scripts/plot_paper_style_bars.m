%{
Filename: plot_paper_style_bars.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    The CNN bar panels redrawn in the style of the published figures, now with
    the CORRECTED control as an extra set of bars. Bar color = pooling layer
    (as the published Fig 3d), bar edge = condition (plain = self, dashed =
    non-PT, thick solid = PT rotation only, dotted = control). Each figure is
    drawn twice: in the full unit space (_full) and with each side reduced to
    its top 5 PCs (_k5, whole layers). Error bars are the SEM over all repeat
    estimates, as in the old scripts (they understate the uncertainty).

Output (results/figures/02_paper_style_bars/):
    fig3d_cnn_cue_transfer_{full,k5}   AlexNet (top) and VGG-19 (bottom), per rendering pair;
                                       full space at 96 (AlexNet) / 64 (VGG) units as published,
                                       k = 5 on whole layers (non-PT = source on the target's PCs)
    hierarchy_by_training_{full,k5}    AlexNet pool2 <-> pool5, per training rendering, averaged
                                       over the 3 test renderings (the Fig 4b layout)
    hierarchy_all_combos_{full,k5}     the same for all 9 train x test combinations (Fig 12 layout)
    paper_style_bars.csv               every plotted bar (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '02_paper_style_bars';
fprintf('plot_paper_style_bars\n');

rows = {};
for space = {'full', 'k5'}
    rows = [rows; cue_bars(R, S, space{1}, out_dir)];               %#ok<AGROW>
    rows = [rows; hierarchy_bars(R, S, space{1}, false, out_dir)];  %#ok<AGROW>
    rows = [rows; hierarchy_bars(R, S, space{1}, true, out_dir)];   %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'figure', 'panel', 'group', 'series', 'mean', 'sem'});
save_table(T, S, out_dir, 'paper_style_bars');

%% ===========================================================================
function rows = cue_bars(R, S, space, out_dir)
groups = {'EC-AC', 'EC-EX', 'AC-EX'};
F  = {'self', 'nopt', 'pt_rot', 'ctrl_m'};
E4 = {'self', 'nopt', 'pt', 'ctrl'};
if strcmp(space, 'full')
    N4 = {'self', 'non-PT', 'PT (rotation only)', 'PT control (corrected)'};
else
    N4 = {'self', 'non-PT (target PCs)', 'PT (rotation only)', 'PT control (corrected)'};
end
fig = figure('Visible', 'off', 'Position', [100 100 1700 1000], 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
name = ['fig3d_cnn_cue_transfer_' space];
for in = 1:2
    net  = S.nets{in};
    keys = S.keys(strcmp(S.layers(:, 1), net));
    series = struct('name', {}, 'color', {}, 'edge', {}, 'alpha', {});
    for f = 1:numel(F)
        for j = 1:numel(keys)
            series(end + 1) = struct('name', sprintf('%s %s', strrep(keys{j}(numel(net) + 2:end), 'pool', 'pool '), N4{f}), ...
                'color', S.c.layer.(keys{j}), 'edge', E4{f}, 'alpha', []); %#ok<AGROW>
        end
    end
    [M, E] = deal(zeros(numel(groups), numel(series)));
    for g = 1:numel(groups)
        for f = 1:numel(F)
            for j = 1:numel(keys)
                [M(g, (f - 1) * numel(keys) + j), E(g, (f - 1) * numel(keys) + j)] = mean_sem(pick_cue(R, S, space, keys{j}, F{f}, groups{g}));
            end
        end
    end
    ax = nexttile(tl);
    plot_grouped_bars(ax, M, E, groups, series, S, 1.1);
    xlabel(ax, 'Stimulus pair');
    if strcmp(space, 'full')
        title(ax, sprintf('%s cue-transfer decoding, full unit space (%d units per layer, as published)', ...
            S.net_label.(net), S.fig3d_N.(net)));
    else
        title(ax, sprintf('%s cue-transfer decoding, top %d PCs per rendering (whole layers)', S.net_label.(net), S.k_main));
    end
    legend(ax, [{series.name}, {'chance (2%)'}], 'Location', 'eastoutside', 'FontSize', 8);
    for g = 1:numel(groups)
        for i = 1:numel(series)
            rows(end + 1, :) = {name, S.net_label.(net), groups{g}, series(i).name, M(g, i), E(g, i)}; %#ok<AGROW>
        end
    end
end
save_figure(fig, S, out_dir, name);
end

function x = pick_cue(R, S, space, key, field, group)
if strcmp(space, 'full')
    r  = R.cue.(key);
    iN = find(r.N_list == S.fig3d_N.(r.net));
    P  = strcmp(r.pair_group, group);
    x  = r.(field)(iN, P, 1:r.n_draws(iN), :);
else
    r  = R.pc_cue.(key);
    ik = find(r.k_list == S.k_main);
    P  = strcmp(r.pair_group, group);
    x  = r.(field)(ik, P, :);
end
end

%% ===========================================================================
function rows = hierarchy_bars(R, S, space, all_combos, out_dir)
cues = {'ac', 'ex', 'ec'};                          % group order of the published Fig 4b
if strcmp(space, 'full'), r = R.hier; else, r = R.pc_hier; end
if all_combos
    [tc, sc] = ndgrid(1:3, 1:3);  tc = tc';  sc = sc';   % training rendering slowest
    groups = arrayfun(@(i) [upper(cues{tc(i)}) ' ' char(8594) ' ' upper(cues{sc(i)})], 1:9, 'UniformOutput', false);
else
    groups = cellfun(@(c) [upper(c) ' ' char(8594) ' AC|EX|EC'], cues, 'UniformOutput', false);
end
if strcmp(space, 'full')
    F = {'self', 'nopt', 'pt_rot', 'ctrl_m'};  E4 = {'self', 'nopt', 'pt', 'ctrl'};
    N4 = {'self', 'non-PT (different units)', 'PT (rotation only)', 'PT control (corrected)'};
else
    F = {'self', 'pt_rot', 'ctrl_m'};  E4 = {'self', 'pt', 'ctrl'};
    N4 = {'self', 'PT (rotation only)', 'PT control (corrected)'};
end
dir_name = {'pool2-to-pool5', 'pool5-to-pool2'};
self_name = {'pool2 self-decoding', 'pool5 self-decoding'};       % target layer of each direction
series = struct('name', {}, 'color', {}, 'edge', {}, 'alpha', {});
for f = 1:numel(F)
    for di = 1:2
        if strcmp(E4{f}, 'self'), nm = self_name{di}; else, nm = sprintf('%s %s', dir_name{di}, N4{f}); end
        series(end + 1) = struct('name', nm, 'color', S.c.dir(di, :), 'edge', E4{f}, 'alpha', []); %#ok<AGROW>
    end
end
[M, E] = deal(zeros(numel(groups), numel(series)));
for di = 1:2
    for g = 1:numel(groups)
        if all_combos
            C = strcmp(r.combo_target_cue, cues{tc(g)}) & strcmp(r.combo_source_cue, cues{sc(g)});
        else
            C = strcmp(r.combo_target_cue, cues{g});
        end
        for f = 1:numel(F)
            if strcmp(space, 'full')
                x = r.(F{f})(di, C, :, :);
            else
                x = r.(F{f})(di, C, r.k_list == S.k_main, :);
            end
            [M(g, 2 * f - 2 + di), E(g, 2 * f - 2 + di)] = mean_sem(x);
        end
    end
end
fig = figure('Visible', 'off', 'Position', [100 100 1500 560], 'Color', 'w');
if all_combos, fig.Position(3) = 1800; end
ax = axes(fig);
plot_grouped_bars(ax, M, E, groups, series, S, 1.1);
if strcmp(space, 'full')
    info = sprintf('full unit space (pool2 subsampled to the %d units of pool5)', r.N_match);
else
    info = sprintf('top %d PCs per layer (whole layers: pool2 %d, pool5 %d units)', S.k_main, r.N_layer(1), r.N_layer(2));
end
title(ax, sprintf('AlexNet layer-to-layer decoding, %s; "A-to-B" = decoder trained on A, tested on B', info));
if all_combos
    xlabel(ax, ['Training rendering ' char(8594) ' test rendering']);
    name = ['hierarchy_all_combos_' space];
else
    xlabel(ax, ['Training rendering ' char(8594) ' test renderings']);
    name = ['hierarchy_by_training_' space];
end
legend(ax, [{series.name}, {'chance (2%)'}], 'Location', 'eastoutside', 'FontSize', 8);
save_figure(fig, S, out_dir, name);
rows = {};
for g = 1:numel(groups)
    for i = 1:numel(series)
        rows(end + 1, :) = {name, 'AlexNet pool2 <-> pool5', groups{g}, series(i).name, M(g, i), E(g, i)}; %#ok<AGROW>
    end
end
end
