%{
Filename: plot_dimensionality.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Dimensionality of the Gabor models, from the participation ratio (PR) stored
    by the drivers (population_pr.m): raw PR of the z-scored trial-averaged
    means (the rebuttal definition) and the noise-corrected PR from
    cross-validated PCA. Raw solid, noise-corrected dashed; colors as Fig 15.

Output (results/figures/08_dimensionality/):
    pr_vs_cells        raw and noise-corrected PR vs number of cells (mean over the 3 renderings;
                       SEM over cell draws)
    cumvar_vs_k        cumulative signal variance vs number of PCs, whole models (mean over renderings)
    metrics_vs_pr      PT - control and headroom (full space) vs noise-corrected PR, one point per size
    pr_vs_cells.csv    PR per model, size and rendering (mean and SEM over draws)
    cumvar_vs_k.csv    cumulative variance per model, rendering and k
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = gabor_figure_settings();
R = load_gabor_results(S);
out_dir = '08_dimensionality';
fprintf('plot_dimensionality\n');

%% PR table: per rendering and averaged over renderings (per draw), SEM over draws
rows = {};
for j = 1:2
    m = S.models{j};  r = R.cue.(m);
    for iN = 1:numel(r.N_list)
        D = r.n_draws(iN);
        for c = 0:3                                           % 0 = mean over renderings
            if c == 0
                xr = squeeze(mean(r.pr_raw(iN, :, 1:D), 2));  xc = squeeze(mean(r.pr_cv(iN, :, 1:D), 2));  cname = 'mean';
            else
                xr = squeeze(r.pr_raw(iN, c, 1:D));           xc = squeeze(r.pr_cv(iN, c, 1:D));           cname = upper(r.cue_names{c});
            end
            [mr, er] = mean_sem(xr);  [mc, ec] = mean_sem(xc);
            rows(end + 1, :) = {S.label.(m), r.N_list(iN), cname, mr, er, mc, ec, D}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'model', 'N', 'rendering', 'pr_raw', 'pr_raw_sem', 'pr_cv', 'pr_cv_sem', 'n_draws'});
save_table(T, S, out_dir, 'pr_vs_cells');

%% PR vs number of cells
fig = figure('Visible', 'off', 'Position', [100 100 900 560], 'Color', 'w');
ax = axes(fig);  hold(ax, 'on');
h = gobjects(0);  nm = {};
for j = 1:2
    m = S.models{j};  c = S.c.model.(m);
    sel = strcmp(T.model, S.label.(m)) & strcmp(T.rendering, 'mean');
    h(end + 1) = errorbar(ax, T.N(sel), T.pr_raw(sel), T.pr_raw_sem(sel), 'Color', c, 'LineStyle', '-', ...
        'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<SAGROW>
    h(end + 1) = errorbar(ax, T.N(sel), T.pr_cv(sel), T.pr_cv_sem(sel), 'Color', c, 'LineStyle', '--', ...
        'LineWidth', 1.5, 'Marker', '.', 'MarkerSize', 10, 'CapSize', 3); %#ok<SAGROW>
    nm = [nm, {[S.label.(m) ' raw'], [S.label.(m) ' noise-corrected']}]; %#ok<AGROW>
end
box(ax, 'on');  xlim(ax, [0 200]);  xlabel(ax, 'Number of model cells');  ylabel(ax, 'Participation ratio');
legend(ax, h, nm, 'Location', 'southeast');
title(ax, {'Participation ratio of the trial-averaged signal manifold vs number of cells', ...
    '(mean over AC, EC, EX; SEM over cell draws)'}, 'FontWeight', 'normal');
save_figure(fig, S, out_dir, 'pr_vs_cells');

%% cumulative variance vs k (whole models, mean over renderings)
fig = figure('Visible', 'off', 'Position', [100 100 900 560], 'Color', 'w');
ax = axes(fig);  hold(ax, 'on');
h = gobjects(0);  nm = {};  at5 = [];  rows = {};
for j = 1:2
    m = S.models{j};  r = R.pc.(m);  c = S.c.model.(m);
    nk  = max(cellfun(@numel, [r.cumvar, r.cv_cumvar]));
    pad = @(v) [v(:); repmat(v(end), nk - numel(v), 1)];   % a spectrum ends at its rank: the curve stays at 1
    raw = mean(cell2mat(cellfun(pad, r.cumvar, 'UniformOutput', false)), 2);
    cv  = mean(cell2mat(cellfun(pad, r.cv_cumvar, 'UniformOutput', false)), 2);
    h(end + 1) = plot(ax, 1:nk, raw, 'Color', c, 'LineStyle', '-', 'LineWidth', 1.5); %#ok<SAGROW>
    h(end + 1) = plot(ax, 1:nk, cv,  'Color', c, 'LineStyle', '--', 'LineWidth', 1.5); %#ok<SAGROW>
    nm = [nm, {[S.label.(m) ' raw'], [S.label.(m) ' noise-corrected']}]; %#ok<AGROW>
    at5 = [at5; raw(S.k_main), cv(S.k_main)]; %#ok<AGROW>
    for cc = 1:3
        for k = 1:nk
            rows(end + 1, :) = {S.label.(m), upper(r.cue_names{cc}), k, val(r.cumvar{cc}, k), val(r.cv_cumvar{cc}, k)}; %#ok<SAGROW>
        end
    end
end
xline(ax, S.k_main, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
text(ax, S.k_main * 1.08, 0.06, sprintf('k = %d: simple %.0f%%, complex %.0f%% (raw)', S.k_main, 100 * at5(1, 1), 100 * at5(2, 1)), 'FontSize', 9);
set(ax, 'XScale', 'log');  xlim(ax, [1 50]);  ylim(ax, [0 1.02]);  set(ax, 'XTick', [1 2 5 10 20 49]);
box(ax, 'on');  xlabel(ax, 'Number of PCs, k');  ylabel(ax, 'Cumulative fraction of signal variance');
legend(ax, h, nm, 'Location', 'northwest');
title(ax, 'Variance of the class means captured by the top k PCs, whole models (mean over AC, EC, EX)', 'FontWeight', 'normal');
save_figure(fig, S, out_dir, 'cumvar_vs_k');
save_table(cell2table(rows, 'VariableNames', {'model', 'rendering', 'k', 'cumvar_raw', 'cumvar_noise_corrected'}), ...
    S, out_dir, 'cumvar_vs_k');

%% PT - control and headroom vs noise-corrected PR (full space, along the cell sweep)
fig = figure('Visible', 'off', 'Position', [100 100 1300 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {'pt_minus_ctrl', 'PT - control'; 'headroom', 'Headroom'};
for s = 1:2
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, 2);
    for j = 1:2
        m = S.models{j};  r = R.cue.(m);  c = S.c.model.(m);
        pooled = r.summary(strcmp(r.summary.pair, 'pooled'), :);
        h(j) = plot(ax, pooled.pr_cv, pooled.(spec{s, 1}), '-o', 'Color', c, 'LineWidth', 1.5, 'MarkerSize', 4, 'MarkerFaceColor', c);
    end
    yline(ax, 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    box(ax, 'on');  xlabel(ax, 'Participation ratio (noise-corrected)');  ylabel(ax, spec{s, 2});
    if s == 1, legend(ax, h, {S.label.simple, S.label.complex}, 'Location', 'northwest'); end
end
title(tl, 'Gabor-model cue transfer (full space) vs dimensionality: each point is one number of cells (pooled over pairs)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'metrics_vs_pr');

%% ===========================================================================
function v = val(x, k)
% k-th entry of a spectrum, NaN beyond its end
if k <= numel(x), v = x(k); else, v = NaN; end
end
