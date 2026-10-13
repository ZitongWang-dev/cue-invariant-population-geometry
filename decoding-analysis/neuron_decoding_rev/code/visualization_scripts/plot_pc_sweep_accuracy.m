%{
Filename: plot_pc_sweep_accuracy.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Decoding accuracy vs the number of retained PCs per side k (log axis,
    k = 1-10, 15-45, 49), for the three comparisons of the main figures:
    self-decoding, PT (rotation only), corrected control and non-PT (source
    trials projected onto the target basis; cue transfer only). Mean +- SEM over
    all repeat estimates; dotted line at k = 5.

Output (results/figures/06_pc_sweep_accuracy/):
    pc_cue_accuracy_vs_k        cue transfer, whole populations, one panel per population
    pc_halves_accuracy_vs_k     halves of floor(N/2) neurons, one panel per population
    pc_across_accuracy_vs_k     across areas, one panel per monkey and direction
    pc_sweep_accuracy.csv       every plotted point (mean, SEM)
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '06_pc_sweep_accuracy';
fprintf('plot_pc_sweep_accuracy\n');

rows = [acc_figure(R, S, 'cue', out_dir); acc_figure(R, S, 'halves', out_dir); acc_figure(R, S, 'across', out_dir)];
T = cell2table(rows, 'VariableNames', {'comparison', 'panel', 'k', 'series', 'mean', 'sem'});
save_table(T, S, out_dir, 'pc_sweep_accuracy');

%% ===========================================================================
function rows = acc_figure(R, S, kind, out_dir)
LS = {'-', '--', ':', '-.'};  LW = [1.5 1.5 1.8 1];
F  = {'self', 'pt_rot', 'ctrl_m', 'nopt'};
LN = {'self-decoding', 'PT (rotation only)', 'PT control (corrected)', 'non-PT (target basis)'};
fig = figure('Visible', 'off', 'Position', [100 100 1400 800], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
if strcmp(kind, 'across')
    panels = {};
    for im = 1:2, for di = 1:2, panels(end + 1, :) = {S.monkeys{im}, di}; end, end %#ok<AGROW>
else
    panels = [S.pops(:), num2cell(zeros(4, 1))];
end
rows = {};
for ip = 1:size(panels, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(0);  nm = {};
    switch kind
        case 'cue'
            r = R.pc_cue.(panels{ip, 1});  c = S.c.(panels{ip, 1}(4:5));
            getv = @(f, ik) r.(f)(ik, :, :);
            ttl = sprintf('%s cue transfer (%d neurons)', strrep(panels{ip, 1}, '_', ' '), r.N);
            pname = strrep(panels{ip, 1}, '_', ' ');
        case 'halves'
            r = R.pc_halves.(panels{ip, 1});  c = S.c.(panels{ip, 1}(4:5));
            getv = @(f, ik) r.(f)(ik, :, :, :);
            ttl = sprintf('%s halves (%d + %d neurons)', strrep(panels{ip, 1}, '_', ' '), r.half_size, r.half_size);
            pname = strrep(panels{ip, 1}, '_', ' ');
        case 'across'
            r = R.pc_across.(panels{ip, 1});  di = panels{ip, 2};  c = S.c.(sprintf('V%d', di));
            getv = @(f, ik) r.(f)(di, :, ik, :);
            ttl = sprintf('%s %s (V1 %d, V2 %d neurons)', panels{ip, 1}, r.dir_labels{di}, r.N_area(1), r.N_area(2));
            pname = sprintf('%s %s', panels{ip, 1}, r.dir_labels{di});
    end
    k = r.k_list(:)';
    for f = 1:numel(F)
        [mu, se] = deal(nan(size(k)));
        for ik = 1:numel(k)
            v = getv(F{f}, ik);
            if any(~isnan(v(:))), [mu(ik), se(ik)] = mean_sem(v); end
        end
        if all(isnan(mu)), continue; end              % no non-PT for halves / across areas
        h(end + 1) = errorbar(ax, k, mu, se, 'Color', c, 'LineStyle', LS{f}, 'LineWidth', LW(f), ...
            'Marker', '.', 'MarkerSize', 10, 'CapSize', 2); %#ok<AGROW>
        nm{end + 1} = LN{f}; %#ok<AGROW>
        for ik = 1:numel(k)
            rows(end + 1, :) = {kind, pname, k(ik), LN{f}, mu(ik), se(ik)}; %#ok<AGROW>
        end
    end
    yline(ax, S.chance, '--', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, S.k_main, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');  xlim(ax, [0.9 55]);  set(ax, 'XTick', [1 2 5 10 20 49]);
    ylim(ax, [0 max(0.3, ceil(max(ylim(ax)) * 10) / 10)]);
    box(ax, 'on');  xlabel(ax, 'Retained PCs per side, k');  ylabel(ax, 'Decoding accuracy');  title(ax, ttl);
    if ip == 1, legend(ax, h, nm, 'Location', 'northwest'); end
end
names = struct('cue', 'Cue transfer', 'halves', 'Within-area halves', 'across', 'Across areas');
title(tl, sprintf('%s vs number of PCs (mean +- SEM over repeat estimates; dotted line: k = %d)', names.(kind), S.k_main), 'FontSize', 11);
save_figure(fig, S, out_dir, sprintf('pc_%s_accuracy_vs_k', kind));
end
