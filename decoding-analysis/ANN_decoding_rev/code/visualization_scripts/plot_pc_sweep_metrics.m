%{
Filename: plot_pc_sweep_metrics.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Normalized transfer metrics of the CNN layers (PT / control, PT - control,
    headroom; pooled means, transfer_metrics.m) vs the number of retained PCs
    per side k (log axis), whole layers. Fixing k fixes the rotation's freedom
    (k(k-1)/2 parameters) whatever the layer size, so the curves compare layers
    at equal freedom. Shading is the percentile-bootstrap 95% CI over folds. The
    markers at "full" give the full-space value of the same whole layer.
    AlexNet solid lines / circles, VGG-19 dashed / squares; layer colors as
    published.

Output (results/figures/07_pc_sweep_metrics/):
    pc_cue_metrics_vs_k         cue transfer, all eight layers
    pc_hierarchy_metrics_vs_k   AlexNet pool2 <-> pool5, per direction and both pooled
    pc_sweep_metrics.csv        every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '07_pc_sweep_metrics';
fprintf('plot_pc_sweep_metrics\n');
seed = 30000;

%% per-k metrics with bootstrap CIs, and the full-space reference of the same units
L = struct('name', {}, 'K', {}, 'color', {}, 'style', {}, 'marker', {}, 'full', {});
for i = 1:numel(S.keys)
    key = S.keys{i};  r = R.pc_cue.(key);  rf = R.cue.(key);  net = S.layers{i, 1};
    [K, seed] = sweep_metrics(r.k_list, @(ik) deal(r.self(ik, :, :), r.pt_rot(ik, :, :), r.ctrl_m(ik, :, :)), 3, S, seed);
    iN = numel(rf.N_list);  seed = seed + 1;
    Bf = bootstrap_metrics(rf.self(iN, :, 1, :), rf.pt_rot(iN, :, 1, :), rf.ctrl_m(iN, :, 1, :), 4, S, seed);
    L(end + 1) = struct('name', r.population, 'K', K, 'color', S.c.layer.(key), 'style', S.ls.net.(net), ...
        'marker', S.mk.net.(net), 'full', Bf); %#ok<SAGROW>
end
metrics_figure(S, L, 'CNN cue transfer vs number of PCs per rendering, whole layers; "full" = full unit space, same units', ...
    out_dir, 'pc_cue_metrics_vs_k');
rows = add_rows({}, 'cue', L);

r = R.pc_hier;  rf = R.hier;
sets  = {1, 2, [1 2]};
names = {'pool2-to-pool5', 'pool5-to-pool2', 'both directions'};
cols  = [S.c.dir; 0.35 0.35 0.35];
Lh = struct('name', {}, 'K', {}, 'color', {}, 'style', {}, 'marker', {}, 'full', {});
for s = 1:3
    di = sets{s};
    [K, seed] = sweep_metrics(r.k_list, @(ik) deal(r.self(di, :, ik, :), r.pt_rot(di, :, ik, :), r.ctrl_m(di, :, ik, :)), 4, S, seed);
    seed = seed + 1;
    Bf = bootstrap_metrics(rf.self(di, :, :, :), rf.pt_rot(di, :, :, :), rf.ctrl_m(di, :, :, :), 4, S, seed);
    Lh(end + 1) = struct('name', names{s}, 'K', K, 'color', cols(s, :), 'style', '-', 'marker', 'o', 'full', Bf); %#ok<SAGROW>
end
metrics_figure(S, Lh, sprintf(['AlexNet pool2 %s pool5 vs number of PCs per layer (whole layers); per-direction headroom ' ...
    'uses the training layer''s self-decoding; "full" = full unit space'], char(8596)), out_dir, 'pc_hierarchy_metrics_vs_k');
rows = add_rows(rows, 'hierarchy', Lh);

T = cell2table(rows, 'VariableNames', {'comparison', 'layer', 'k', 'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'headroom', 'headroom_lo', 'headroom_hi'});
save_table(T, S, out_dir, 'pc_sweep_metrics');

%% ===========================================================================
function [K, seed] = sweep_metrics(k_list, getter, unit_dim, S, seed)
% metrics and bootstrap CIs at every k; getter(ik) returns the self, PT and control arrays
n = numel(k_list);
K = struct('k', k_list(:)', 'self', zeros(1, n), 'pt', zeros(1, n), 'ctrl', zeros(1, n), ...
    'ratio', zeros(1, n), 'ratio_ci', zeros(n, 2), 'diff', zeros(1, n), 'diff_ci', zeros(n, 2), ...
    'headroom', zeros(1, n), 'headroom_ci', zeros(n, 2));
for ik = 1:n
    [Sa, Pa, Ca] = getter(ik);
    seed = seed + 1;
    B = bootstrap_metrics(Sa, Pa, Ca, unit_dim, S, seed);
    K.self(ik) = B.self;  K.pt(ik) = B.pt;  K.ctrl(ik) = B.ctrl;
    for f = {'ratio', 'diff', 'headroom'}
        K.(f{1})(ik) = B.(f{1});
        K.([f{1} '_ci'])(ik, :) = B.([f{1} '_ci']);
    end
end
end

function rows = add_rows(rows, kind, L)
for i = 1:numel(L)
    K = L(i).K;  B = L(i).full;
    for ik = 1:numel(K.k)
        rows(end + 1, :) = {kind, L(i).name, K.k(ik), K.self(ik), K.pt(ik), K.ctrl(ik), K.ratio(ik), K.ratio_ci(ik, 1), ...
            K.ratio_ci(ik, 2), K.diff(ik), K.diff_ci(ik, 1), K.diff_ci(ik, 2), K.headroom(ik), K.headroom_ci(ik, 1), ...
            K.headroom_ci(ik, 2)}; %#ok<AGROW>
    end
    rows(end + 1, :) = {kind, L(i).name, NaN, B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
        B.diff, B.diff_ci(1), B.diff_ci(2), B.headroom, B.headroom_ci(1), B.headroom_ci(2)}; %#ok<AGROW>   % k = NaN: full space
end
end

function metrics_figure(S, L, ttl, out_dir, name)
spec = {'ratio', 'PT / control', 1; 'diff', 'PT - control (accuracy)', 0; 'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0};
fig = figure('Visible', 'off', 'Position', [100 100 1700 560], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
k_full = 75;
for s = 1:3
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(L));
    for i = 1:numel(L)
        K = L(i).K;  ci = K.([spec{s, 1} '_ci']);
        fill(ax, [K.k, fliplr(K.k)], [ci(:, 1)', fliplr(ci(:, 2)')], L(i).color, 'FaceAlpha', 0.10, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
        h(i) = plot(ax, K.k, K.(spec{s, 1}), 'Color', L(i).color, 'LineStyle', L(i).style, 'LineWidth', 1.5, ...
            'Marker', L(i).marker, 'MarkerSize', 3, 'MarkerFaceColor', L(i).color);
        v = L(i).full.(spec{s, 1});  cf = L(i).full.([spec{s, 1} '_ci']);
        xo = k_full * (1 + 0.035 * (i - (numel(L) + 1) / 2));
        errorbar(ax, xo, v, v - cf(1), cf(2) - v, L(i).marker, 'Color', L(i).color, 'MarkerFaceColor', L(i).color, ...
            'MarkerSize', 5, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
    yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    xline(ax, S.k_main, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    xline(ax, 60, ':', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    set(ax, 'XScale', 'log');  xlim(ax, [0.9 95]);
    set(ax, 'XTick', [1 2 5 10 20 49 k_full], 'XTickLabel', {'1', '2', '5', '10', '20', '49', 'full'});
    box(ax, 'on');  xlabel(ax, 'Retained PCs per side, k');  ylabel(ax, spec{s, 2});
    if s == 1
        lg = legend(ax, h, {L.name}, 'Orientation', 'horizontal', 'NumColumns', min(8, numel(L)));
        lg.Layout.Tile = 'south';
    end
end
title(tl, ttl, 'FontSize', 11);
subtitle(tl, 'lines: pooled means; shading: bootstrap 95% CI over folds; markers at "full": full unit space', 'FontSize', 9);
save_figure(fig, S, out_dir, name);
end
