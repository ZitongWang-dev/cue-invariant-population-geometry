%{
Filename: plot_pc_sweep_metrics.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Normalized transfer metrics (PT / control, headroom, PT - control; pooled
    means, transfer_metrics.m) vs the number of retained PCs per side k (log
    axis), for the three comparisons of the main figures. Shading is the
    percentile-bootstrap 95% CI (folds for cue transfer and across areas, which
    use whole populations once; half-splits for halves). Where the same neurons
    exist in full space, a "full" marker at the right of the axis gives the
    full-space value (filled KO, open FR).

Output (results/figures/07_pc_sweep_metrics/):
    pc_cue_metrics_vs_k             per population (+ full space, whole populations)
    pc_cue_metrics_vs_k_per_pair    PT/control and headroom per rendering pair
    pc_halves_metrics_vs_k          per population (+ full space, halves of floor(N/2))
    pc_across_metrics_vs_k          per monkey: each direction and both pooled
    pc_headroom_all_comparisons     cue transfer vs halves vs across areas, per monkey
    pc_sweep_metrics.csv            every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '07_pc_sweep_metrics';
fprintf('plot_pc_sweep_metrics\n');
seed = 20000;
pops = S.pops;
mstyle = struct('KO', '-', 'FR', '--');

%% ---- per-k metrics with bootstrap CIs ----
M = struct();
for ip = 1:numel(pops)
    r = R.pc_cue.(pops{ip});
    [M.cue.(pops{ip}), seed] = sweep_metrics(r.k_list, @(ik, P) deal(r.self(ik, P, :), r.pt_rot(ik, P, :), r.ctrl_m(ik, P, :)), ...
        1:numel(r.pair_labels), 3, S, seed);
    for g = {'EC-AC', 'EC-EX', 'AC-EX'}
        P = find(strcmp(r.pair_group, g{1}));
        [M.cue_pair.(pops{ip}).(strrep(g{1}, '-', '_')), seed] = sweep_metrics(r.k_list, ...
            @(ik, P) deal(r.self(ik, P, :), r.pt_rot(ik, P, :), r.ctrl_m(ik, P, :)), P, 3, S, seed);
    end
    r = R.pc_halves.(pops{ip});
    [M.halves.(pops{ip}), seed] = sweep_metrics(r.k_list, @(ik, P) deal(r.self(ik, P, :, :), r.pt_rot(ik, P, :, :), r.ctrl_m(ik, P, :, :)), ...
        1:numel(r.renderings), 3, S, seed);
end
dir_sets = {1, 2, [1 2]};
for im = 1:2
    m = S.monkeys{im};  r = R.pc_across.(m);
    for di = 1:3
        % [dirs x combos x k x folds] at one k: unit = folds (dim 4)
        [M.across.(m){di}, seed] = sweep_metrics(r.k_list, @(ik, P) deal(r.self(P, :, ik, :), r.pt_rot(P, :, ik, :), ...
            r.ctrl_m(P, :, ik, :)), dir_sets{di}, 4, S, seed);
    end
end

% full-space references with the same neurons: whole population (cue), floor(N/2) halves
for ip = 1:numel(pops)
    r = R.cue.(pops{ip});  iN = numel(r.N_list);  seed = seed + 1;
    FULL.cue.(pops{ip}) = bootstrap_metrics(r.self(iN, :, 1, :), r.pt_rot(iN, :, 1, :), r.ctrl_m(iN, :, 1, :), 4, S, seed);
    r = R.halves.(pops{ip});  iH = numel(r.h_list);  seed = seed + 1;
    FULL.halves.(pops{ip}) = bootstrap_metrics(r.self(iH, :, :, :), r.pt_rot(iH, :, :, :), r.ctrl_m(iH, :, :, :), 3, S, seed);
end

%% ---- metrics vs k, one figure per comparison ----
lines_cue = struct('name', {}, 'K', {}, 'color', {}, 'style', {}, 'full', {});
for ip = 1:numel(pops)
    lines_cue(end + 1) = struct('name', strrep(pops{ip}, '_', ' '), 'K', M.cue.(pops{ip}), ...
        'color', S.c.(pops{ip}(4:5)), 'style', mstyle.(pops{ip}(1:2)), 'full', FULL.cue.(pops{ip})); %#ok<SAGROW>
end
metrics_figure(S, lines_cue, 'Cue transfer vs number of PCs per cue (whole populations); "full" = full neuron space, same neurons', ...
    out_dir, 'pc_cue_metrics_vs_k');
lines_h = lines_cue;
for ip = 1:numel(pops)
    lines_h(ip).K = M.halves.(pops{ip});  lines_h(ip).full = FULL.halves.(pops{ip});
end
metrics_figure(S, lines_h, 'Within-area halves (floor(N/2) neurons each) vs number of PCs per half; "full" = full neuron space, same halves sizes', ...
    out_dir, 'pc_halves_metrics_vs_k');
lines_a = struct('name', {}, 'K', {}, 'color', {}, 'style', {}, 'full', {});
dname = {'V1-to-V2', 'V2-to-V1', 'both directions'};
dcol  = [S.c.V1; S.c.V2; 0.35 0.35 0.35];
for im = 1:2
    m = S.monkeys{im};
    for di = 1:3
        lines_a(end + 1) = struct('name', sprintf('%s %s', m, dname{di}), 'K', M.across.(m){di}, ...
            'color', dcol(di, :), 'style', mstyle.(m), 'full', []); %#ok<SAGROW>
    end
end
metrics_figure(S, lines_a, ['Across areas vs number of PCs per area (whole populations); per-direction headroom uses ' ...
    'the target area''s self-decoding and can exceed 1'], out_dir, 'pc_across_metrics_vs_k');

%% ---- cue transfer per rendering pair ----
fig = figure('Visible', 'off', 'Position', [100 100 1700 700], 'Color', 'w');
tl  = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
gnames = {'EC_AC', 'EC_EX', 'AC_EX'};
spec = {'ratio', 'PT / control', 1; 'headroom', 'Headroom', 0};
for s = 1:2
    for ip = 1:numel(pops)
        ax = nexttile(tl);  hold(ax, 'on');
        h = gobjects(1, 3);
        for g = 1:3
            h(g) = band_line(ax, M.cue_pair.(pops{ip}).(gnames{g}), spec{s, 1}, S.c.pairs(g, :), '-');
        end
        finish_k_axis(ax, spec{s, 3});
        ylabel(ax, spec{s, 2});
        if s == 1, title(ax, strrep(pops{ip}, '_', ' ')); end
        if s == 1 && ip == 1, legend(ax, h, strrep(gnames, '_', '-'), 'Location', 'best'); end
    end
end
title(tl, 'Cue transfer per rendering pair vs number of PCs (shading: bootstrap 95% CI over folds)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'pc_cue_metrics_vs_k_per_pair');

%% ---- all comparisons, per monkey ----
fig = figure('Visible', 'off', 'Position', [100 100 1400 800], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:2
    for im = 1:2
        m = S.monkeys{im};
        ax = nexttile(tl);  hold(ax, 'on');
        h = gobjects(0);  nm = {};
        for a = 1:2
            pop = sprintf('%s_V%d', m, a);  c = S.c.(sprintf('V%d', a));
            h(end + 1) = band_line(ax, M.cue.(pop), spec{s, 1}, c, '-');            nm{end + 1} = sprintf('cue transfer V%d', a); %#ok<SAGROW>
            h(end + 1) = band_line(ax, M.halves.(pop), spec{s, 1}, c, '--');        nm{end + 1} = sprintf('halves V%d', a); %#ok<SAGROW>
        end
        h(end + 1) = band_line(ax, M.across.(m){3}, spec{s, 1}, [0.35 0.35 0.35], ':');  nm{end + 1} = 'across areas (both directions)'; %#ok<SAGROW>
        finish_k_axis(ax, spec{s, 3});
        ylabel(ax, spec{s, 2});
        if s == 1, title(ax, m); legend(ax, h, nm, 'Location', 'best'); end
    end
end
title(tl, 'PC sweep, all comparisons: halves (same rendering, different neurons) vs cue transfer vs across areas', 'FontSize', 11);
save_figure(fig, S, out_dir, 'pc_headroom_all_comparisons');

%% ---- table ----
rows = {};
for ip = 1:numel(pops)
    rows = add_rows(rows, 'cue', strrep(pops{ip}, '_', ' '), M.cue.(pops{ip}));
    for g = 1:3
        rows = add_rows(rows, 'cue per pair', sprintf('%s %s', strrep(pops{ip}, '_', ' '), strrep(gnames{g}, '_', '-')), ...
            M.cue_pair.(pops{ip}).(gnames{g}));
    end
end
for ip = 1:numel(pops), rows = add_rows(rows, 'halves', strrep(pops{ip}, '_', ' '), M.halves.(pops{ip})); end
for im = 1:2
    m = S.monkeys{im};
    for di = 1:3, rows = add_rows(rows, 'across', sprintf('%s %s', m, dname{di}), M.across.(m){di}); end
end
T = cell2table(rows, 'VariableNames', {'comparison', 'population', 'k', 'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', ...
    'headroom', 'headroom_lo', 'headroom_hi', 'pt_minus_ctrl', 'diff_lo', 'diff_hi'});
save_table(T, S, out_dir, 'pc_sweep_metrics');

%% ===========================================================================
function [K, seed] = sweep_metrics(k_list, getter, P, unit_dim, S, seed)
% metrics and bootstrap CIs at every k; getter(ik, P) returns the self, PT and control arrays
n = numel(k_list);
K = struct('k', k_list(:)', 'self', zeros(1, n), 'pt', zeros(1, n), 'ctrl', zeros(1, n), ...
    'ratio', zeros(1, n), 'ratio_ci', zeros(n, 2), 'headroom', zeros(1, n), 'headroom_ci', zeros(n, 2), ...
    'diff', zeros(1, n), 'diff_ci', zeros(n, 2));
for ik = 1:n
    [Sa, Pa, Ca] = getter(ik, P);
    seed = seed + 1;
    B = bootstrap_metrics(Sa, Pa, Ca, unit_dim, S, seed);
    K.self(ik) = B.self;  K.pt(ik) = B.pt;  K.ctrl(ik) = B.ctrl;
    for f = {'ratio', 'headroom', 'diff'}
        K.(f{1})(ik) = B.(f{1});
        K.([f{1} '_ci'])(ik, :) = B.([f{1} '_ci']);
    end
end
end

function rows = add_rows(rows, kind, name, K)
for ik = 1:numel(K.k)
    rows(end + 1, :) = {kind, name, K.k(ik), K.self(ik), K.pt(ik), K.ctrl(ik), K.ratio(ik), K.ratio_ci(ik, 1), ...
        K.ratio_ci(ik, 2), K.headroom(ik), K.headroom_ci(ik, 1), K.headroom_ci(ik, 2), K.diff(ik), ...
        K.diff_ci(ik, 1), K.diff_ci(ik, 2)}; %#ok<AGROW>
end
end

function metrics_figure(S, L, ttl, out_dir, name)
% three panels (PT/control, headroom, PT - control) vs k, one line per entry of L
spec = {'ratio', 'PT / control', 1; 'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0; 'diff', 'PT - control (accuracy)', 0};
fig = figure('Visible', 'off', 'Position', [100 100 1650 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
k_full = 75;    % x position of the full-space markers
has_full = any(arrayfun(@(x) ~isempty(x.full), L));
for s = 1:3
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(L));
    for i = 1:numel(L)
        h(i) = band_line(ax, L(i).K, spec{s, 1}, L(i).color, L(i).style);
        if ~isempty(L(i).full)
            v = L(i).full.(spec{s, 1});  ci = L(i).full.([spec{s, 1} '_ci']);
            xo = k_full * (1 + 0.04 * (i - (numel(L) + 1) / 2));
            mf = 'none';  if strcmp(L(i).style, '-'), mf = L(i).color; end
            errorbar(ax, xo, v, v - ci(1), ci(2) - v, 'o', 'Color', L(i).color, 'MarkerFaceColor', mf, ...
                'MarkerSize', 6, 'LineWidth', 1.2, 'HandleVisibility', 'off');
        end
    end
    finish_k_axis(ax, spec{s, 3});
    if has_full
        xline(ax, 60, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        xlim(ax, [0.9 95]);  set(ax, 'XTick', [1 2 5 10 20 49 k_full], 'XTickLabel', {'1', '2', '5', '10', '20', '49', 'full'});
    end
    ylabel(ax, spec{s, 2});
    if s == 1
        lg = legend(ax, h, {L.name}, 'Orientation', 'horizontal', 'NumColumns', min(6, numel(L)));
        lg.Layout.Tile = 'south';
    end
end
title(tl, ttl, 'FontSize', 11);
if has_full
    subtitle(tl, 'lines: pooled means; shading: bootstrap 95% CI; markers at "full": full neuron space (filled KO, open FR)', 'FontSize', 9);
else
    subtitle(tl, 'lines: pooled means; shading: bootstrap 95% CI over folds', 'FontSize', 9);
end
save_figure(fig, S, out_dir, name);
end

function h = band_line(ax, K, f, c, ls)
ci = K.([f '_ci']);
fill(ax, [K.k, fliplr(K.k)], [ci(:, 1)', fliplr(ci(:, 2)')], c, 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h = plot(ax, K.k, K.(f), 'Color', c, 'LineStyle', ls, 'LineWidth', 1.6, 'Marker', '.', 'MarkerSize', 10);
end

function finish_k_axis(ax, ref)
yline(ax, ref, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
xline(ax, 5, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
set(ax, 'XScale', 'log');  xlim(ax, [0.9 55]);  set(ax, 'XTick', [1 2 5 10 20 49]);
box(ax, 'on');  xlabel(ax, 'Retained PCs per side, k');
end
