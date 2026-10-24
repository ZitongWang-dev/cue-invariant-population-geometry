%{
Filename: plot_paper_style_bars.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Fig 15 redrawn in the style of the published figure, now with the CORRECTED
    control as an extra pair of bars. Bar color = model (simple blue, complex
    orange, as published), bar edge = condition (plain = self, dashed = non-PT,
    thick solid = PT rotation only, dotted = control), groups = rendering pairs.
    Error bars are the SEM over all repeat estimates, as in the old scripts
    (they understate the uncertainty).

Output (results/figures/02_paper_style_bars/):
    fig15_gabor_cue_transfer_full   (a) both models at 48 cells (simple down-sampled to the
                                    complex-cell count, as published); (b) whole models
                                    (simple 192, complex 48)
    fig15_gabor_cue_transfer_k5     whole models, top 5 PCs per rendering
                                    (non-PT = source on the target's PCs)
    paper_style_bars.csv            every plotted bar (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '02_paper_style_bars';
fprintf('plot_paper_style_bars\n');

groups = {'EC-AC', 'EC-EX', 'AC-EX'};
F  = {'self', 'nopt', 'pt_rot', 'ctrl_m'};
E4 = {'self', 'nopt', 'pt', 'ctrl'};
rows = {};

%% full space: (a) 48 cells, (b) whole models
fig = figure('Visible', 'off', 'Position', [100 100 1250 900], 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
N4 = {'self-decoding', 'non-PT', 'PT (rotation only)', 'PT control (corrected)'};
designs = {'a', @(m) S.fig15a_N, 'both models at 48 cells (simple down-sampled to the complex-cell count)'; ...
           'b', @(m) R.cue.(m).N_total, 'whole models (simple 192, complex 48 cells)'};
for id = 1:2
    series = make_series(S, E4, N4);
    [M, E] = deal(zeros(numel(groups), numel(series)));
    for g = 1:numel(groups)
        for f = 1:numel(F)
            for j = 1:2
                m  = S.models{j};  r = R.cue.(m);
                iN = find(r.N_list == designs{id, 2}(m));
                P  = strcmp(r.pair_group, groups{g});
                [M(g, 2 * f - 2 + j), E(g, 2 * f - 2 + j)] = mean_sem(r.(F{f})(iN, P, 1:r.n_draws(iN), :));
            end
        end
    end
    ax = nexttile(tl);
    plot_grouped_bars(ax, M, E, groups, series, S, 1.1);
    xlabel(ax, 'Stimulus pair');
    title(ax, sprintf('(%s) Gabor-model cue-transfer decoding, full space: %s', designs{id, 1}, designs{id, 3}));
    if id == 1
        lg = legend(ax, [{series.name}, {'chance (2%)'}]);  lg.Layout.Tile = 'east';
    end
    rows = add_rows(rows, sprintf('fig15%s_full', designs{id, 1}), groups, series, M, E);
end
save_figure(fig, S, out_dir, 'fig15_gabor_cue_transfer_full');

%% top 5 PCs, whole models
N4 = {'self-decoding', 'non-PT (target PCs)', 'PT (rotation only)', 'PT control (corrected)'};
series = make_series(S, E4, N4);
[M, E] = deal(zeros(numel(groups), numel(series)));
for g = 1:numel(groups)
    for f = 1:numel(F)
        for j = 1:2
            r  = R.pc.(S.models{j});  ik = find(r.k_list == S.k_main);
            P  = strcmp(r.pair_group, groups{g});
            [M(g, 2 * f - 2 + j), E(g, 2 * f - 2 + j)] = mean_sem(r.(F{f})(ik, P, :));
        end
    end
end
fig = figure('Visible', 'off', 'Position', [100 100 1250 480], 'Color', 'w');
ax = axes(fig);
plot_grouped_bars(ax, M, E, groups, series, S, 1.1);
xlabel(ax, 'Stimulus pair');
title(ax, sprintf('Gabor-model cue-transfer decoding, top %d PCs per rendering (whole models: simple 192, complex 48 cells)', S.k_main));
legend(ax, [{series.name}, {'chance (2%)'}], 'Location', 'eastoutside');
rows = add_rows(rows, 'fig15_k5', groups, series, M, E);
save_figure(fig, S, out_dir, 'fig15_gabor_cue_transfer_k5');

T = cell2table(rows, 'VariableNames', {'figure', 'group', 'series', 'mean', 'sem'});
save_table(T, S, out_dir, 'paper_style_bars');

%% ===========================================================================
function series = make_series(S, edges, names)
% condition-major order, as the published figure: self simple, self complex, non-PT simple, ...
series = struct('name', {}, 'color', {}, 'edge', {}, 'alpha', {});
for f = 1:numel(edges)
    for j = 1:numel(S.models)
        m = S.models{j};
        series(end + 1) = struct('name', sprintf('%s %s', S.label.(m), names{f}), 'color', S.c.model.(m), ...
            'edge', edges{f}, 'alpha', []); %#ok<AGROW>
    end
end
end

function rows = add_rows(rows, name, groups, series, M, E)
for g = 1:numel(groups)
    for i = 1:numel(series)
        rows(end + 1, :) = {name, groups{g}, series(i).name, M(g, i), E(g, i)}; %#ok<AGROW>
    end
end
end
