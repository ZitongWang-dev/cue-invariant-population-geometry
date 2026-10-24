%{
Filename: plot_paper_style_curves.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Gabor-model cue-transfer decoding vs the number of model cells, full space,
    in the style of the published incremental figures (Figs 4c, 4d), with the
    CORRECTED control added: self-decoding (solid), PT rotation only (dashed),
    control (dotted) and non-PT (dash-dot), one panel per model, colors as
    the published Fig 15. Error bars are the SEM over all repeat estimates, as
    in the old scripts. The dotted vertical line marks the Fig 15a size (48).
    This is the Gabor curve of rebuttal point 9, Fig 1, now with the control.

Output (results/figures/03_paper_style_curves/):
    gabor_vs_cells              self, PT, control and non-PT vs number of cells
    paper_style_curves.csv      every plotted point (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '03_paper_style_curves';
fprintf('plot_paper_style_curves\n');

F  = {'self', 'pt_rot', 'ctrl_m', 'nopt'};
LN = {'self-decoding', 'PT (rotation only)', 'PT control (corrected)', 'non-PT'};
LS = {'-', '--', ':', '-.'};
LW = [1.5 1.5 1.8 1];

fig = figure('Visible', 'off', 'Position', [100 100 1400 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
for j = 1:2
    m = S.models{j};  r = R.cue.(m);  c = S.c.model.(m);
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(F));
    x = r.N_list(:)';
    for f = 1:numel(F)
        [mu, se] = deal(zeros(size(x)));
        for iN = 1:numel(x)
            [mu(iN), se(iN)] = mean_sem(r.(F{f})(iN, :, 1:r.n_draws(iN), :));
            rows(end + 1, :) = {S.label.(m), x(iN), LN{f}, mu(iN), se(iN)}; %#ok<SAGROW>
        end
        h(f) = errorbar(ax, x, mu, se, 'Color', c, 'LineStyle', LS{f}, 'LineWidth', LW(f), ...
            'Marker', '.', 'MarkerSize', 10, 'CapSize', 3);
    end
    xline(ax, S.fig15a_N, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    yline(ax, S.chance, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    ylim(ax, [0 1.05]);  xlim(ax, [0, r.N_total * 1.03]);  box(ax, 'on');
    xlabel(ax, 'Number of model cells');  ylabel(ax, 'Decoding accuracy');
    title(ax, sprintf('%s (%d cells)', S.label.(m), r.N_total));
    if j == 1, legend(ax, h, LN, 'Location', 'east'); end
end
title(tl, 'Gabor-model cue transfer vs number of cells, full space (mean +- SEM over repeat estimates; dotted line: Fig 15a size, 48)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'gabor_vs_cells');
save_table(cell2table(rows, 'VariableNames', {'model', 'N', 'series', 'mean', 'sem'}), S, out_dir, 'paper_style_curves');
