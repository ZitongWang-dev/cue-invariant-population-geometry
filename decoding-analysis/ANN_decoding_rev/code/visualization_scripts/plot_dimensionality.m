%{
Filename: plot_dimensionality.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Dimensionality of the CNN layers, from the participation ratio (PR) stored
    by the drivers (population_pr.m): raw PR of the z-scored trial-averaged
    means (the rebuttal definition) and the noise-corrected PR from
    cross-validated PCA. Raw solid, noise-corrected dashed; layer colors as
    published.

Output (results/figures/08_dimensionality/):
    pr_vs_units        raw and noise-corrected PR vs number of units (mean over the 3 renderings;
                       SEM over unit draws)
    cumvar_vs_k        cumulative signal variance vs number of PCs, whole layers (mean over renderings)
    metrics_vs_pr      PT - control and headroom of cue transfer (full space) vs noise-corrected PR,
                       one point per unit count
    pr_vs_units.csv    PR per layer, unit count and rendering (mean and SEM over draws)
    cumvar_vs_k.csv    cumulative variance per layer, rendering and k
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '08_dimensionality';
fprintf('plot_dimensionality\n');

%% PR table: per rendering and averaged over renderings (per draw), SEM over draws
rows = {};
for i = 1:numel(S.keys)
    r = R.cue.(S.keys{i});
    for iN = 1:numel(r.N_list)
        D = r.n_draws(iN);
        for c = 0:3                                           % 0 = mean over renderings
            if c == 0
                xr = squeeze(mean(r.pr_raw(iN, :, 1:D), 2));  xc = squeeze(mean(r.pr_cv(iN, :, 1:D), 2));  cname = 'mean';
            else
                xr = squeeze(r.pr_raw(iN, c, 1:D));           xc = squeeze(r.pr_cv(iN, c, 1:D));           cname = upper(r.cue_names{c});
            end
            [mr, er] = mean_sem(xr);  [mc, ec] = mean_sem(xc);
            rows(end + 1, :) = {r.population, r.N_list(iN), cname, mr, er, mc, ec, D}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'layer', 'N', 'rendering', 'pr_raw', 'pr_raw_sem', 'pr_cv', 'pr_cv_sem', 'n_draws'});
save_table(T, S, out_dir, 'pr_vs_units');

%% PR vs number of units
fig = figure('Visible', 'off', 'Position', [100 100 1500 560], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for in = 1:2
    net = S.nets{in};
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(0);  nm = {};
    for i = find(strcmp(S.layers(:, 1), net))'
        key = S.keys{i};  c = S.c.layer.(key);  pname = R.cue.(key).population;
        sel = strcmp(T.layer, pname) & strcmp(T.rendering, 'mean');
        h(end + 1) = errorbar(ax, T.N(sel), T.pr_raw(sel), T.pr_raw_sem(sel), 'Color', c, 'LineStyle', '-', ...
            'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<SAGROW>
        h(end + 1) = errorbar(ax, T.N(sel), T.pr_cv(sel), T.pr_cv_sem(sel), 'Color', c, 'LineStyle', '--', ...
            'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<SAGROW>
        nm = [nm, {sprintf('pool %d raw', S.layers{i, 2}), sprintf('pool %d noise-corrected', S.layers{i, 2})}]; %#ok<AGROW>
    end
    n_max = max(cellfun(@(k) R.cue.(k).N_total, S.keys(strcmp(S.layers(:, 1), net))));
    xlim(ax, [0, n_max * 1.02]);  box(ax, 'on');  xlabel(ax, 'Number of units');  ylabel(ax, 'Participation ratio');
    title(ax, sprintf('%s', S.net_label.(net)));
    legend(ax, h, nm, 'Location', 'southeast', 'NumColumns', 2, 'FontSize', 8);
end
title(tl, 'Participation ratio of the trial-averaged signal manifold vs number of units (mean over AC, EC, EX; SEM over unit draws)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'pr_vs_units');

%% cumulative variance vs k (whole layers, mean over renderings)
fig = figure('Visible', 'off', 'Position', [100 100 1500 560], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
rows = {};
for in = 1:2
    net = S.nets{in};
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(0);  nm = {};  at5 = [];
    for i = find(strcmp(S.layers(:, 1), net))'
        key = S.keys{i};  r = R.pc_cue.(key);  c = S.c.layer.(key);
        nk  = max(cellfun(@numel, [r.cumvar, r.cv_cumvar]));
        pad = @(v) [v(:); repmat(v(end), nk - numel(v), 1)];   % a spectrum ends at its rank: the curve stays at 1
        raw = mean(cell2mat(cellfun(pad, r.cumvar, 'UniformOutput', false)), 2);
        cv  = mean(cell2mat(cellfun(pad, r.cv_cumvar, 'UniformOutput', false)), 2);
        h(end + 1) = plot(ax, 1:nk, raw, 'Color', c, 'LineStyle', '-', 'LineWidth', 1.5); %#ok<SAGROW>
        h(end + 1) = plot(ax, 1:nk, cv,  'Color', c, 'LineStyle', '--', 'LineWidth', 1.5); %#ok<SAGROW>
        nm = [nm, {sprintf('pool %d raw', S.layers{i, 2}), sprintf('pool %d noise-corrected', S.layers{i, 2})}]; %#ok<AGROW>
        at5 = [at5; raw(S.k_main), cv(S.k_main)]; %#ok<AGROW>
        for cc = 1:3
            for k = 1:nk
                rows(end + 1, :) = {r.population, upper(r.cue_names{cc}), k, val(r.cumvar{cc}, k), val(r.cv_cumvar{cc}, k)}; %#ok<SAGROW>
            end
        end
    end
    xline(ax, S.k_main, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    text(ax, S.k_main * 1.08, 0.06, sprintf('k = %d: raw %.0f-%.0f%%, noise-corrected %.0f-%.0f%%', S.k_main, ...
        100 * min(at5(:, 1)), 100 * max(at5(:, 1)), 100 * min(at5(:, 2)), 100 * max(at5(:, 2))), 'FontSize', 9);
    set(ax, 'XScale', 'log');  xlim(ax, [1 50]);  ylim(ax, [0 1.02]);  set(ax, 'XTick', [1 2 5 10 20 49]);
    box(ax, 'on');  xlabel(ax, 'Number of PCs, k');  ylabel(ax, 'Cumulative fraction of signal variance');
    title(ax, sprintf('%s (whole layers)', S.net_label.(net)));
    legend(ax, h, nm, 'Location', 'northwest', 'NumColumns', 2, 'FontSize', 8);
end
title(tl, 'Variance of the trial-averaged class means captured by the top k PCs (mean over AC, EC, EX; raw and cross-validated spectra)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'cumvar_vs_k');
save_table(cell2table(rows, 'VariableNames', {'layer', 'rendering', 'k', 'cumvar_raw', 'cumvar_noise_corrected'}), ...
    S, out_dir, 'cumvar_vs_k');

%% PT - control and headroom vs noise-corrected PR (full space, along the unit sweep)
fig = figure('Visible', 'off', 'Position', [100 100 1500 560], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {'pt_minus_ctrl', 'PT - control', 0; 'headroom', 'Headroom', 0};
for s = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(S.keys));
    for i = 1:numel(S.keys)
        r = R.cue.(S.keys{i});  net = S.layers{i, 1};  c = S.c.layer.(S.keys{i});
        pooled = r.summary(strcmp(r.summary.pair, 'pooled'), :);
        h(i) = plot(ax, pooled.pr_cv, pooled.(spec{s, 1}), 'Color', c, 'LineStyle', S.ls.net.(net), 'LineWidth', 1.5, ...
            'Marker', S.mk.net.(net), 'MarkerSize', 4, 'MarkerFaceColor', c);
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    box(ax, 'on');  xlabel(ax, 'Participation ratio (noise-corrected)');  ylabel(ax, spec{s, 2});
    if s == 1
        lg = legend(ax, h, cellfun(@(k) R.cue.(k).population, S.keys, 'UniformOutput', false), ...
            'Orientation', 'horizontal', 'NumColumns', 8);
        lg.Layout.Tile = 'south';
    end
end
title(tl, 'CNN cue transfer (full space) vs dimensionality: each point is one unit count (pooled over pairs; PR rises with the unit count)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'metrics_vs_pr');

%% ===========================================================================
function v = val(x, k)
% k-th entry of a spectrum, NaN beyond its end
if k <= numel(x), v = x(k); else, v = NaN; end
end
