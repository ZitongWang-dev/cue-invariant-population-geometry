%{
Filename: plot_paper_style_bars.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    The bar panels of the main figures redrawn in the style of the published
    figures (and of the old visualization scripts), now with the CORRECTED
    control as an extra pair of bars. KO on top, FR below, as in the paper.
    Each figure is drawn twice: in the full neuron space (_full, the published
    design) and with each side reduced to its top 5 PCs (_k5, whole populations;
    halves of floor(N/2); no non-PT bars for halves and across areas, whose two
    sides are different neurons).
    Error bars are the SEM over all repeat estimates (draws x ordered pairs x
    folds), as in the old scripts; they understate the uncertainty because those
    estimates are not independent.

Output (results/figures/02_paper_style_bars/):
    fig3b_cue_transfer_{full,k5}              per rendering pair; full space: V1 whole,
                                              V2 subsampled to the V1 count (FR 48, KO 109)
    fig4a_halves_{full,k5}                    per rendering; full space: halves of 24 (FR) / 54 (KO)
    fig4b_across_area_{full,k5}               per training rendering, averaged over the 3 test renderings
    fig12_across_area_all_combos_{full,k5}    every train x test rendering combination
    paper_style_bars.csv                      every plotted bar (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '02_paper_style_bars';
fprintf('plot_paper_style_bars\n');

rows = {};
for space = {'full', 'k5'}
    rows = [rows; cue_bars(R, S, space{1}, out_dir)];               %#ok<AGROW>
    rows = [rows; halves_bars(R, S, space{1}, out_dir)];            %#ok<AGROW>
    rows = [rows; across_bars(R, S, space{1}, false, out_dir)];     %#ok<AGROW>
    rows = [rows; across_bars(R, S, space{1}, true, out_dir)];      %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'figure', 'monkey', 'group', 'series', 'mean', 'sem'});
save_table(T, S, out_dir, 'paper_style_bars');

%% ===========================================================================
function rows = cue_bars(R, S, space, out_dir)
groups = {'EC-AC', 'EC-EX', 'AC-EX'};
F = {'self', 'nopt', 'pt_rot', 'ctrl_m'};
series = make_series(S, {'self', 'nopt', 'pt', 'ctrl'}, ...
    {'self-decoding', 'non-PT-decoding', 'PT-decoding (rotation only)', 'PT control (corrected)'}, {'V1', 'V2'});
[Ms, Es, ttl] = deal(cell(1, 2));
for im = 1:2
    m = S.monkeys{im};
    [M, E] = deal(zeros(numel(groups), numel(series)));
    for a = 1:2
        pop = sprintf('%s_V%d', m, a);
        for g = 1:numel(groups)
            for f = 1:numel(F)
                [M(g, 2 * f - 2 + a), E(g, 2 * f - 2 + a)] = mean_sem(pick_cue(R, S, space, pop, F{f}, groups{g}));
            end
        end
    end
    Ms{im} = M;  Es{im} = E;
    if strcmp(space, 'full')
        ttl{im} = sprintf('%s cue-transfer decoding, full neuron space (V1 %d neurons, V2 subsampled to %d)', ...
            m, S.anchor_N.(m), S.anchor_N.(m));
    else
        ttl{im} = sprintf('%s cue-transfer decoding, top %d PCs per cue (whole populations: V1 %d, V2 %d neurons)', ...
            m, S.k_main, R.pc_cue.([m '_V1']).N, R.pc_cue.([m '_V2']).N);
    end
end
rows = draw_and_save(S, Ms, Es, groups, series, ttl, 'Stimulus pair', [], ['fig3b_cue_transfer_' space], out_dir);
end

%% ===========================================================================
function rows = halves_bars(R, S, space, out_dir)
groups = {'AC', 'EX', 'EC'};
if strcmp(space, 'full')
    F = {'self', 'nopt', 'pt_rot', 'ctrl_m'};  E4 = {'self', 'nopt', 'pt', 'ctrl'};
    N4 = {'self-decoding', 'non-PT-decoding (different neurons)', 'PT-decoding (rotation only)', 'PT control (corrected)'};
else
    F = {'self', 'pt_rot', 'ctrl_m'};  E4 = {'self', 'pt', 'ctrl'};
    N4 = {'self-decoding', 'PT-decoding (rotation only)', 'PT control (corrected)'};
end
series = make_series(S, E4, N4, {'V1', 'V2'});
[Ms, Es, ttl] = deal(cell(1, 2));
for im = 1:2
    m = S.monkeys{im};
    [M, E] = deal(zeros(numel(groups), numel(series)));
    hs = zeros(1, 2);
    for a = 1:2
        pop = sprintf('%s_V%d', m, a);
        for g = 1:numel(groups)
            for f = 1:numel(F)
                [x, hs(a)] = pick_halves(R, S, space, pop, F{f}, groups{g});
                [M(g, 2 * f - 2 + a), E(g, 2 * f - 2 + a)] = mean_sem(x);
            end
        end
    end
    Ms{im} = M;  Es{im} = E;
    if strcmp(space, 'full')
        ttl{im} = sprintf('%s within-area halves, full neuron space (two halves of %d neurons in each area)', m, hs(1));
    else
        ttl{im} = sprintf('%s within-area halves, top %d PCs per half (halves of V1 %d, V2 %d neurons)', m, S.k_main, hs(1), hs(2));
    end
end
rows = draw_and_save(S, Ms, Es, groups, series, ttl, 'Rendering', [], ['fig4a_halves_' space], out_dir);
end

%% ===========================================================================
function rows = across_bars(R, S, space, all_combos, out_dir)
cues = {'ac', 'ex', 'ec'};                          % group order of the published figure
if all_combos
    [tc, sc] = ndgrid(1:3, 1:3);  tc = tc';  sc = sc';   % training rendering slowest
    groups = arrayfun(@(i) [upper(cues{tc(i)}) ' ' char(8594) ' ' upper(cues{sc(i)})], 1:9, 'UniformOutput', false);
else
    groups = cellfun(@(c) [upper(c) ' ' char(8594) ' AC|EX|EC'], cues, 'UniformOutput', false);
end
if strcmp(space, 'full')
    F = {'self', 'nopt', 'pt_rot', 'ctrl_m'};  E4 = {'self', 'nopt', 'pt', 'ctrl'};
    N4 = {'self', 'non-PT (different neurons)', 'PT (rotation only)', 'PT control (corrected)'};
else
    F = {'self', 'pt_rot', 'ctrl_m'};  E4 = {'self', 'pt', 'ctrl'};
    N4 = {'self', 'PT (rotation only)', 'PT control (corrected)'};
end
% blue = V1-to-V2 (decoder trained on V1, tested on V2), orange = V2-to-V1; self = target area
series = make_series(S, E4, N4, {'V1', 'V2'});
dir_name = {'V1-to-V2', 'V2-to-V1'};
for i = 1:numel(series)
    a = 2 - mod(i, 2);                              % odd series = direction 1 (blue), even = 2 (orange)
    if strcmp(series(i).edge, 'self')
        series(i).name = sprintf('V%d self-decoding', a);
    else
        series(i).name = sprintf('%s %s', dir_name{a}, N4{ceil(i / 2)});
    end
end
[Ms, Es, ttl] = deal(cell(1, 2));
for im = 1:2
    m = S.monkeys{im};
    [M, E] = deal(zeros(numel(groups), numel(series)));
    for di = 1:2
        for g = 1:numel(groups)
            if all_combos, tcue = cues{tc(g)}; scue = cues{sc(g)}; else, tcue = cues{g}; scue = ''; end
            for f = 1:numel(F)
                [x, info] = pick_across(R, S, space, m, di, F{f}, tcue, scue);
                [M(g, 2 * f - 2 + di), E(g, 2 * f - 2 + di)] = mean_sem(x);
            end
        end
    end
    Ms{im} = M;  Es{im} = E;
    ttl{im} = sprintf('%s across-area decoding, %s', m, info);
end
if all_combos
    rows = draw_and_save(S, Ms, Es, groups, series, ttl, ['Training rendering ' char(8594) ' test rendering'], 1700, ...
        ['fig12_across_area_all_combos_' space], out_dir);
else
    rows = draw_and_save(S, Ms, Es, groups, series, ttl, ['Training rendering ' char(8594) ' test renderings'], [], ...
        ['fig4b_across_area_' space], out_dir);
end
end

%% ===========================================================================
function x = pick_cue(R, S, space, pop, field, group)
m = pop(1:2);
if strcmp(space, 'full')
    r  = R.cue.(pop);
    iN = find(r.N_list == S.anchor_N.(m));
    P  = strcmp(r.pair_group, group);
    x  = r.(field)(iN, P, 1:r.n_draws(iN), :);
else
    r  = R.pc_cue.(pop);
    ik = find(r.k_list == S.k_main);
    P  = strcmp(r.pair_group, group);
    x  = r.(field)(ik, P, :);
end
end

function [x, h] = pick_halves(R, S, space, pop, field, group)
m = pop(1:2);
if strcmp(space, 'full')
    r  = R.halves.(pop);
    iH = find(r.h_list == S.anchor_h.(m));
    ri = strcmp(r.renderings, lower(group));
    x  = r.(field)(iH, ri, :, :);
    h  = r.h_list(iH);
else
    r  = R.pc_halves.(pop);
    ik = find(r.k_list == S.k_main);
    ri = strcmp(r.renderings, lower(group));
    x  = r.(field)(ik, ri, :, :);
    h  = r.half_size;
end
end

function [x, info] = pick_across(R, S, space, m, di, field, tcue, scue)
if strcmp(space, 'full')
    r = R.across.(m);
else
    r = R.pc_across.(m);
end
C = strcmp(r.combo_target_cue, tcue);
if ~isempty(scue), C = C & strcmp(r.combo_source_cue, scue); end
if strcmp(space, 'full')
    x = r.(field)(di, C, :, :);
    info = sprintf('full neuron space (V1 %d neurons, V2 subsampled to %d)', r.N_match, r.N_match);
else
    ik = find(r.k_list == S.k_main);
    x = r.(field)(di, C, ik, :);
    info = sprintf('top %d PCs per area (whole populations: V1 %d, V2 %d neurons)', S.k_main, r.N_area(1), r.N_area(2));
end
end

function series = make_series(S, edges, names, areas)
series = struct('name', {}, 'color', {}, 'edge', {}, 'alpha', {});
for f = 1:numel(edges)
    for a = 1:numel(areas)
        series(end + 1) = struct('name', sprintf('%s %s', areas{a}, names{f}), 'color', S.c.(areas{a}), ...
            'edge', edges{f}, 'alpha', []); %#ok<AGROW>
    end
end
end

function rows = draw_and_save(S, Ms, Es, groups, series, ttl, xl, width, name, out_dir)
% two stacked panels (KO top, FR bottom) with a common y axis; returns the plotted bars as table rows
fig = figure('Visible', 'off', 'Position', S.fig2, 'Color', 'w');
if ~isempty(width), fig.Position(3) = width; end
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
y_max = max(cellfun(@(M, E) max(M + E, [], 'all'), Ms, Es));
y_max = max(0.3, ceil(y_max * 1.08 * 10) / 10);
rows = {};
for im = 1:2
    ax = nexttile(tl);
    plot_grouped_bars(ax, Ms{im}, Es{im}, groups, series, S, y_max);
    title(ax, ttl{im});
    xlabel(ax, xl);
    if im == 1
        lg = legend(ax, [{series.name}, {'chance (2%)'}]);
        lg.Layout.Tile = 'east';
    end
    for g = 1:numel(groups)
        for i = 1:numel(series)
            rows(end + 1, :) = {name, S.monkeys{im}, groups{g}, series(i).name, Ms{im}(g, i), Es{im}(g, i)}; %#ok<AGROW>
        end
    end
end
save_figure(fig, S, out_dir, name);
end
