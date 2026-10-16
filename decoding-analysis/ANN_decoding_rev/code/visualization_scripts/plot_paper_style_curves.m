%{
Filename: plot_paper_style_curves.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Fig 4d (CNN decoding vs number of units) redrawn in the published style,
    full unit space, with the CORRECTED control added: per pooling layer,
    self-decoding (solid), PT rotation only (dashed) and control (dotted), AlexNet
    on top, VGG-19 below, layer colors as published. Error bars are the SEM over
    all repeat estimates, as in the old scripts. A second figure gives each
    layer its own panel with non-PT as well. The dotted vertical line marks the
    Fig 3d size (96 / 64 units).

Output (results/figures/03_paper_style_curves/):
    fig4d_cnn_vs_units            published layout: self, PT and control per layer
    fig4d_cnn_vs_units_per_layer  one panel per layer: self, PT, control and non-PT
    paper_style_curves.csv        every plotted point (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '03_paper_style_curves';
fprintf('plot_paper_style_curves\n');

F  = {'self', 'pt_rot', 'ctrl_m', 'nopt'};
LN = {'self-decoding', 'PT (rotation only)', 'PT control (corrected)', 'non-PT'};
LS = {'-', '--', ':', '-.'};
LW = [1.5 1.5 1.8 1];

%% curves of every layer: mean and SEM per unit count
C = struct();
rows = {};
for i = 1:numel(S.keys)
    r = R.cue.(S.keys{i});
    x = r.N_list(:)';
    for f = 1:numel(F)
        [mu, se] = deal(zeros(size(x)));
        for iN = 1:numel(x)
            [mu(iN), se(iN)] = mean_sem(r.(F{f})(iN, :, 1:r.n_draws(iN), :));
            rows(end + 1, :) = {r.population, x(iN), LN{f}, mu(iN), se(iN)}; %#ok<SAGROW>
        end
        C.(S.keys{i}).(F{f}) = struct('x', x, 'mu', mu, 'se', se);
    end
end
save_table(cell2table(rows, 'VariableNames', {'layer', 'N', 'series', 'mean', 'sem'}), S, out_dir, 'paper_style_curves');

%% published layout: one panel per network
fig = figure('Visible', 'off', 'Position', [100 100 1300 950], 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for in = 1:2
    net = S.nets{in};
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(0);  names = {};
    for i = find(strcmp(S.layers(:, 1), net))'
        key = S.keys{i};
        for f = 1:3
            c = C.(key).(F{f});
            h(end + 1) = errorbar(ax, c.x, c.mu, c.se, 'Color', S.c.layer.(key), 'LineStyle', LS{f}, 'LineWidth', LW(f), ...
                'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<SAGROW>
            names{end + 1} = sprintf('pool %d %s', S.layers{i, 2}, LN{f}); %#ok<SAGROW>
        end
    end
    xline(ax, S.fig3d_N.(net), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    yline(ax, S.chance, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    n_max = max(cellfun(@(k) R.cue.(k).N_total, S.keys(strcmp(S.layers(:, 1), net))));
    ylim(ax, [0 1.05]);  xlim(ax, [0, n_max * 1.02]);  box(ax, 'on');
    xlabel(ax, 'Number of units');  ylabel(ax, 'Decoding accuracy');
    title(ax, sprintf('%s cue-transfer decoding vs number of units (full space; dotted line: Fig 3d size, %d units)', ...
        S.net_label.(net), S.fig3d_N.(net)));
    legend(ax, h, names, 'Location', 'eastoutside', 'FontSize', 8);
end
save_figure(fig, S, out_dir, 'fig4d_cnn_vs_units');

%% one panel per layer, all four conditions
fig = figure('Visible', 'off', 'Position', [100 100 2000 800], 'Color', 'w');
tl  = tiledlayout(fig, 2, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
tiles = [1 2 3 6 7 8 9 10];                % AlexNet in row 1, VGG-19 in row 2
for i = 1:numel(S.keys)
    key = S.keys{i};  r = R.cue.(key);
    ax = nexttile(tl, tiles(i));  hold(ax, 'on');
    h = gobjects(1, numel(F));
    for f = 1:numel(F)
        c = C.(key).(F{f});
        h(f) = errorbar(ax, c.x, c.mu, c.se, 'Color', S.c.layer.(key), 'LineStyle', LS{f}, 'LineWidth', LW(f), ...
            'Marker', '.', 'MarkerSize', 10, 'CapSize', 2);
    end
    xline(ax, S.fig3d_N.(r.net), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    yline(ax, S.chance, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    ylim(ax, [0 1.05]);  xlim(ax, [0, r.N_total * 1.03]);  box(ax, 'on');
    title(ax, sprintf('%s (%d units)', r.population, r.N_total));
    xlabel(ax, 'Number of units');  ylabel(ax, 'Decoding accuracy');
    if i == 1, lg = legend(ax, h, LN); lg.Layout.Tile = 4; end
end
title(tl, 'CNN cue transfer vs number of units, per layer (mean +- SEM over repeat estimates; dotted line: Fig 3d size)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'fig4d_cnn_vs_units_per_layer');
