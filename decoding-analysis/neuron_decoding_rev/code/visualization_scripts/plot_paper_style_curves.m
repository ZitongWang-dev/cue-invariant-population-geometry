%{
Filename: plot_paper_style_curves.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    The incremental curves of the paper (Fig 4c, Fig 13) redrawn in the published
    style, full neuron space, with the CORRECTED control added: self-decoding,
    PT (rotation only), control and non-PT vs population size, V1 blue, V2
    orange, KO on top, FR below. Error bars are the SEM over all repeat
    estimates, as in the old scripts.

Output (results/figures/03_paper_style_curves/):
    fig4c_cue_transfer_vs_N    cue transfer vs number of neurons
    fig13_halves_vs_size       within-area halves vs neurons per half
    paper_style_curves.csv     every plotted point (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '03_paper_style_curves';
fprintf('plot_paper_style_curves\n');

rows = [plot_curves(R, S, 'cue', out_dir); plot_curves(R, S, 'halves', out_dir)];
T = cell2table(rows, 'VariableNames', {'figure', 'population', 'size', 'series', 'mean', 'sem'});
save_table(T, S, out_dir, 'paper_style_curves');

%% ===========================================================================
function rows = plot_curves(R, S, kind, out_dir)
F  = {'self', 'pt_rot', 'ctrl_m', 'nopt'};
LS = {'-', '--', ':', '-.'};
LW = [1.5 1.5 1.8 1];
if strcmp(kind, 'cue')
    LN = {'self-decoding', 'PT-decoding (rotation only)', 'PT control (corrected)', 'non-PT-decoding'};
    xl = 'Number of neurons';  name = 'fig4c_cue_transfer_vs_N';
else
    LN = {'self-decoding', 'PT-decoding (rotation only)', 'PT control (corrected)', 'non-PT (different neurons)'};
    xl = 'Neurons per half';   name = 'fig13_halves_vs_size';
end
fig = figure('Visible', 'off', 'Position', S.fig2, 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ymax = 0;
ax = gobjects(1, 2);
rows = {};
for im = 1:2
    m = S.monkeys{im};
    ax(im) = nexttile(tl);  hold(ax(im), 'on');
    h = gobjects(0);  names = {};
    for a = 1:2
        pop = sprintf('%s_V%d', m, a);
        r = R.(kind).(pop);
        if strcmp(kind, 'cue'), x = r.N_list(:)'; else, x = r.h_list(:)'; end
        for f = 1:numel(F)
            [mu, se] = deal(zeros(size(x)));
            for i = 1:numel(x)
                if strcmp(kind, 'cue')
                    v = r.(F{f})(i, :, 1:r.n_draws(i), :);
                else
                    v = r.(F{f})(i, :, :, :);
                end
                [mu(i), se(i)] = mean_sem(v);
                rows(end + 1, :) = {name, pop, x(i), LN{f}, mu(i), se(i)}; %#ok<AGROW>
            end
            h(end + 1) = errorbar(ax(im), x, mu, se, 'Color', S.c.(sprintf('V%d', a)), 'LineStyle', LS{f}, ...
                'LineWidth', LW(f), 'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<AGROW>
            names{end + 1} = sprintf('V%d %s', a, LN{f}); %#ok<AGROW>
            ymax = max(ymax, max(mu + se));
        end
    end
    yline(ax(im), S.chance, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    hold(ax(im), 'off');  box(ax(im), 'on');
    xlabel(ax(im), xl);  ylabel(ax(im), 'Decoding accuracy');
    if strcmp(kind, 'cue')
        title(ax(im), sprintf('%s cue-transfer decoding vs population size (full space)', m));
    else
        title(ax(im), sprintf('%s within-area halves vs half size (full space)', m));
    end
    if im == 1
        lg = legend(ax(im), h, names);  lg.Layout.Tile = 'east';
    end
end
for im = 1:2, ylim(ax(im), [0 ceil(ymax * 10) / 10]); xlim(ax(im), [0 max(xlim(ax(im)))]); end
save_figure(fig, S, out_dir, name);
end
