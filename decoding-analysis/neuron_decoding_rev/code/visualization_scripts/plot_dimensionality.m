%{
Filename: plot_dimensionality.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Dimensionality of the recorded populations, from the participation ratio
    (PR) stored by the drivers: raw PR of the z-scored trial-averaged means (the
    rebuttal definition, inflated by trial noise) and the noise-corrected PR from
    cross-validated PCA (population_pr.m). Raw solid, noise-corrected dashed.

Output (results/figures/08_dimensionality/):
    pr_vs_N            raw and noise-corrected PR vs number of neurons (mean over the 3 renderings)
    pr_per_cue_vs_N    the same per rendering
    cumvar_vs_k        cumulative signal variance vs number of PCs, whole populations (raw and
                       noise-corrected spectra), with the k = 5 level marked
    metrics_vs_pr      PT/control and headroom of cue transfer vs PR (full space, along N)
    pr_vs_N.csv        PR per population, size and rendering (mean and SEM over neuron draws)
    cumvar_vs_k.csv    cumulative variance per population, rendering and k
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '08_dimensionality';
fprintf('plot_dimensionality\n');
pops = S.pops;

%% PR table: per rendering and averaged over renderings (per draw), SEM over draws
rows = {};
for ip = 1:numel(pops)
    r = R.cue.(pops{ip});
    for i = 1:numel(r.N_list)
        D = r.n_draws(i);
        for c = 0:3                                           % 0 = mean over renderings
            if c == 0
                xr = squeeze(mean(r.pr_raw(i, :, 1:D), 2));  xc = squeeze(mean(r.pr_cv(i, :, 1:D), 2));  cname = 'mean';
            else
                xr = squeeze(r.pr_raw(i, c, 1:D));           xc = squeeze(r.pr_cv(i, c, 1:D));           cname = upper(r.cue_names{c});
            end
            [mr, er] = mean_sem(xr);  [mc, ec] = mean_sem(xc);
            rows(end + 1, :) = {strrep(pops{ip}, '_', ' '), r.N_list(i), cname, mr, er, mc, ec, D}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'population', 'N', 'rendering', 'pr_raw', 'pr_raw_sem', 'pr_cv', 'pr_cv_sem', 'n_draws'});
save_table(T, S, out_dir, 'pr_vs_N');

%% PR vs N (mean over renderings)
fig = figure('Visible', 'off', 'Position', S.fig2, 'Color', 'w');
tl  = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for im = 1:2
    m = S.monkeys{im};
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(0);  nm = {};
    for a = 1:2
        pname = sprintf('%s V%d', m, a);  c = S.c.(sprintf('V%d', a));
        sel = strcmp(T.population, pname) & strcmp(T.rendering, 'mean');
        h(end + 1) = errorbar(ax, T.N(sel), T.pr_raw(sel), T.pr_raw_sem(sel), 'Color', c, 'LineStyle', '-', ...
            'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<SAGROW>
        h(end + 1) = errorbar(ax, T.N(sel), T.pr_cv(sel), T.pr_cv_sem(sel), 'Color', c, 'LineStyle', '--', ...
            'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<SAGROW>
        nm = [nm, {sprintf('V%d raw PR', a), sprintf('V%d noise-corrected PR', a)}]; %#ok<AGROW>
    end
    box(ax, 'on');  xlabel(ax, 'Number of neurons');  ylabel(ax, 'Participation ratio');
    title(ax, sprintf('%s: participation ratio of the trial-averaged signal manifold (mean over AC, EC, EX; SEM over neuron draws)', m));
    if im == 1, lg = legend(ax, h, nm);  lg.Layout.Tile = 'east'; end
end
save_figure(fig, S, out_dir, 'pr_vs_N');

%% PR per rendering
fig = figure('Visible', 'off', 'Position', [100 100 1400 800], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for ip = 1:numel(pops)
    r = R.cue.(pops{ip});  pname = strrep(pops{ip}, '_', ' ');
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(0);  nm = {};
    for c = 1:3
        sel = strcmp(T.population, pname) & strcmp(T.rendering, upper(r.cue_names{c}));
        h(end + 1) = plot(ax, T.N(sel), T.pr_raw(sel), 'Color', S.c.cues(c, :), 'LineStyle', '-', 'LineWidth', 1.5, ...
            'Marker', '.', 'MarkerSize', 10); %#ok<SAGROW>
        h(end + 1) = plot(ax, T.N(sel), T.pr_cv(sel), 'Color', S.c.cues(c, :), 'LineStyle', '--', 'LineWidth', 1.5, ...
            'Marker', '.', 'MarkerSize', 10); %#ok<SAGROW>
        nm = [nm, {sprintf('%s raw PR', upper(r.cue_names{c})), sprintf('%s noise-corrected PR', upper(r.cue_names{c}))}]; %#ok<AGROW>
    end
    box(ax, 'on');  xlabel(ax, 'Number of neurons');  ylabel(ax, 'Participation ratio');
    title(ax, pname);
    if ip == 1, legend(ax, h, nm, 'Location', 'northwest', 'NumColumns', 2); end
end
title(tl, 'Participation ratio per rendering vs number of neurons (mean over neuron draws)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'pr_per_cue_vs_N');

%% cumulative variance vs k (whole populations)
fig = figure('Visible', 'off', 'Position', [100 100 1400 800], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
for ip = 1:numel(pops)
    r = R.pc_cue.(pops{ip});
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(0);  nm = {};  at5 = zeros(2, 3);
    for c = 1:3
        cv = {r.cumvar{c}, r.cv_cumvar{c}};
        for t = 1:2
            ls = '-';  if t == 2, ls = '--'; end
            h(end + 1) = plot(ax, 1:numel(cv{t}), cv{t}, 'Color', S.c.cues(c, :), 'LineStyle', ls, 'LineWidth', 1.5); %#ok<SAGROW>
            at5(t, c) = cv{t}(S.k_main);
        end
        nm = [nm, {sprintf('%s raw', upper(r.cue_names{c})), sprintf('%s noise-corrected', upper(r.cue_names{c}))}]; %#ok<AGROW>
        nk = max(numel(cv{1}), numel(cv{2}));                 % the two spectra can differ in length
        cvp = cellfun(@(v) [v(:); nan(nk - numel(v), 1)], cv, 'UniformOutput', false);
        for k = 1:nk
            rows(end + 1, :) = {strrep(pops{ip}, '_', ' '), upper(r.cue_names{c}), k, cvp{1}(k), cvp{2}(k)}; %#ok<SAGROW>
        end
    end
    xline(ax, S.k_main, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    text(ax, S.k_main * 1.08, 0.08, sprintf('k = %d: raw %.0f-%.0f%%, noise-corrected %.0f-%.0f%%', S.k_main, ...
        100 * min(at5(1, :)), 100 * max(at5(1, :)), 100 * min(at5(2, :)), 100 * max(at5(2, :))), 'FontSize', 9);
    set(ax, 'XScale', 'log');  xlim(ax, [1 50]);  ylim(ax, [0 1.02]);  set(ax, 'XTick', [1 2 5 10 20 49]);
    box(ax, 'on');  xlabel(ax, 'Number of PCs, k');  ylabel(ax, 'Cumulative fraction of signal variance');
    title(ax, sprintf('%s (whole population, %d neurons)', strrep(pops{ip}, '_', ' '), r.N));
    if ip == 1
        lg = legend(ax, h, nm, 'Orientation', 'horizontal', 'NumColumns', 6);
        lg.Layout.Tile = 'south';
    end
end
title(tl, 'Variance of the trial-averaged class means captured by the top k PCs (raw and cross-validated spectra)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'cumvar_vs_k');
save_table(cell2table(rows, 'VariableNames', {'population', 'rendering', 'k', 'cumvar_raw', 'cumvar_noise_corrected'}), ...
    S, out_dir, 'cumvar_vs_k');

%% PT/control and headroom vs PR (cue transfer, full space, along N)
fig = figure('Visible', 'off', 'Position', [100 100 1300 900], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {'pt_over_ctrl', 'PT / control', 1; 'headroom', 'Headroom', 0};
xs   = {'pr_raw', 'pr_cv'};  xn = {'raw PR', 'noise-corrected PR'};
for s = 1:2
    for t = 1:2
        ax = nexttile(tl);  hold(ax, 'on');
        h = gobjects(1, numel(pops));
        for ip = 1:numel(pops)
            r = R.cue.(pops{ip});
            pooled = r.summary(strcmp(r.summary.pair, 'pooled'), :);
            ls = '-';  if startsWith(pops{ip}, 'FR'), ls = '--'; end
            c = S.c.(pops{ip}(4:5));
            h(ip) = plot(ax, pooled.(xs{t}), pooled.(spec{s, 1}), 'Color', c, 'LineStyle', ls, 'LineWidth', 1.5, ...
                'Marker', 'o', 'MarkerSize', 4, 'MarkerFaceColor', c);
        end
        yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        box(ax, 'on');  xlabel(ax, sprintf('Participation ratio (%s)', xn{t}));  ylabel(ax, spec{s, 2});
        if s == 1 && t == 1, legend(ax, h, strrep(pops, '_', ' '), 'Location', 'best'); end
    end
end
title(tl, 'Cue transfer (full space) vs dimensionality: each point is one population size N (pooled over pairs)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'metrics_vs_pr');
