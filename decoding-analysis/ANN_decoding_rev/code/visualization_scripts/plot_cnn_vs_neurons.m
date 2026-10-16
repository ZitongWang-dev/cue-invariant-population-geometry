%{
Filename: plot_cnn_vs_neurons.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    CNN layers and recorded neurons on the same axes (the model figure of the
    rerun plan), all from the corrected-control reruns:
      - at k = 5 PCs per side, where the rotation has the same 10 free
        parameters for every population and layer;
      - along the PC sweep (k = 1-49);
      - in full space along the number of units / neurons;
      - against the noise-corrected participation ratio.
    Neurons: V1 blue, V2 orange (KO solid / dark, FR dashed / light). CNN layers:
    AlexNet in greens, VGG-19 in purples, darker = deeper.
    The two kinds of populations sit in different noise regimes: CNN trials
    are Poisson draws around the activations and self-decoding is near 1, while
    neural self-decoding is 0.2-0.5. Headroom divides by self - control, so the
    self-decoding panels are part of every comparison.

    Neuron results are read from neuron_decoding_rev/results/decoding_outputs
    through that folder's own helpers (figure_settings, load_all_results).

Output (results/figures/09_cnn_vs_neurons/):
    k5_cnn_vs_neurons           bars at k = 5, whole populations / layers (bootstrap 95% CI over folds)
    pc_sweep_cnn_vs_neurons     self, PT - control and headroom vs k
    vs_size_cnn_vs_neurons      self, PT - control and headroom vs number of units, full space
    vs_pr_cnn_vs_neurons        PT - control and headroom vs noise-corrected PR, full space
    cnn_vs_neurons_k5.csv       the k = 5 bars
    cnn_vs_neurons_curves.csv   every plotted curve point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S  = ann_figure_settings();
R  = load_ann_results(S);
Sn = figure_settings();             % neuron figures (neuron_decoding_rev helpers)
Rn = load_all_results(Sn);
out_dir = '09_cnn_vs_neurons';
fprintf('plot_cnn_vs_neurons\n');

%% the populations: name, results, color, line style, marker
P = struct('name', {}, 'kind', {}, 'pc', {}, 'cue', {}, 'color', {}, 'ls', {}, 'mk', {}, 'alpha', {});
for p = {'KO_V1', 'KO_V2', 'FR_V1', 'FR_V2'}
    pop = p{1};  isKO = startsWith(pop, 'KO');
    c = Sn.c.(pop(4:5));
    if isKO, ls = '-'; a = 1; else, ls = '--'; a = 0.45; end
    P(end + 1) = struct('name', strrep(pop, '_', ' '), 'kind', 'neurons', 'pc', Rn.pc_cue.(pop), 'cue', Rn.cue.(pop), ...
        'color', c, 'ls', ls, 'mk', 'o', 'alpha', a); %#ok<SAGROW>
end
shades = struct('alexnet', [0.70 0.86 0.55; 0.35 0.65 0.25; 0.08 0.38 0.08], ...
                'vgg',     [0.86 0.78 0.94; 0.72 0.58 0.86; 0.57 0.38 0.76; 0.42 0.20 0.62; 0.27 0.05 0.46]);
for in = 1:2
    net = S.nets{in};
    idx = find(strcmp(S.layers(:, 1), net))';
    for j = 1:numel(idx)
        key = S.keys{idx(j)};
        P(end + 1) = struct('name', R.cue.(key).population, 'kind', 'CNN', 'pc', R.pc_cue.(key), 'cue', R.cue.(key), ...
            'color', shades.(net)(j, :), 'ls', S.ls.net.(net), 'mk', S.mk.net.(net), 'alpha', 1); %#ok<SAGROW>
    end
end
nP = numel(P);

%% 1. bars at k = 5 (whole populations / layers)
seed = 40000;
B = cell(1, nP);
rows = {};
for i = 1:nP
    r = P(i).pc;  ik = find(r.k_list == S.k_main);
    seed = seed + 1;
    B{i} = bootstrap_metrics(r.self(ik, :, :), r.pt_rot(ik, :, :), r.ctrl_m(ik, :, :), 3, S, seed);
    b = B{i};
    rows(end + 1, :) = {P(i).kind, P(i).name, r.N, b.self, b.pt, b.ctrl, b.ratio, b.ratio_ci(1), b.ratio_ci(2), ...
        b.diff, b.diff_ci(1), b.diff_ci(2), b.headroom, b.headroom_ci(1), b.headroom_ci(2)}; %#ok<SAGROW>
end
save_table(cell2table(rows, 'VariableNames', {'kind', 'population', 'n_units', 'self', 'pt', 'ctrl', 'pt_over_ctrl', ...
    'ratio_lo', 'ratio_hi', 'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'headroom', 'headroom_lo', 'headroom_hi'}), ...
    S, out_dir, 'cnn_vs_neurons_k5');

spec = {'self', 'Self-decoding', NaN; 'diff', 'PT - control (accuracy)', 0; 'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0; ...
        'ratio', 'PT / control', 1};
fig = figure('Visible', 'off', 'Position', [100 100 1900 560], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
x = [1:4, 5.6:7.6, 9.2:13.2];                 % gaps between neurons, AlexNet and VGG
for s = 1:size(spec, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    for i = 1:nP
        v = B{i}.(spec{s, 1});
        bh = bar(ax, x(i), v, 0.8, 'FaceColor', P(i).color, 'EdgeColor', 'k', 'FaceAlpha', P(i).alpha);
        if strcmp(P(i).ls, '--'), bh.LineStyle = '--'; end
        if ~strcmp(spec{s, 1}, 'self')
            ci = B{i}.([spec{s, 1} '_ci']);
            errorbar(ax, x(i), v, v - ci(1), ci(2) - v, 'k', 'LineStyle', 'none', 'LineWidth', 1, 'CapSize', 3);
        end
    end
    if ~isnan(spec{s, 3}), yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3]); end
    set(ax, 'XTick', x, 'XTickLabel', {P.name}, 'XTickLabelRotation', 40, 'FontSize', 9);
    xlim(ax, [0.3, max(x) + 0.7]);  ylabel(ax, spec{s, 2});  box(ax, 'on');
end
title(tl, sprintf(['Neurons vs CNN layers at k = %d PCs per side (same rotation freedom), whole populations / layers; ' ...
    'error bars: bootstrap 95%% CI over folds'], S.k_main), 'FontSize', 11);
save_figure(fig, S, out_dir, 'k5_cnn_vs_neurons');

%% curves from the summary tables (pooled over the six ordered pairs)
crow = {};
curve_spec = {'self', 'Self-decoding'; 'pt_minus_ctrl', 'PT - control (accuracy)'; 'headroom', 'Headroom'};

% 2. vs k
[fig, ax, h] = curve_figure(S, P, 1, 3, @(p) pooled(p.pc.summary, 'set'), 'k', curve_spec, true);
for a = 1:numel(ax), xlabel(ax(a), 'Retained PCs per side, k'); xlim(ax(a), [0.9 55]); xline(ax(a), S.k_main, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); end
finish(fig, ax, h, P, sprintf(['Neurons vs CNN layers along the PC sweep (whole populations / layers); ' ...
    'dotted line: k = %d'], S.k_main));
save_figure(fig, S, out_dir, 'pc_sweep_cnn_vs_neurons');
crow = add_curves(crow, P, 'k', @(p) pooled(p.pc.summary, 'set'), 'k', curve_spec);

% 3. vs number of units, full space
[fig, ax, h] = curve_figure(S, P, 1, 3, @(p) pooled(p.cue.summary, 'pair'), 'N', curve_spec, true);
for a = 1:numel(ax), xlabel(ax(a), 'Number of neurons / units'); xlim(ax(a), [9 560]); end
finish(fig, ax, h, P, 'Neurons vs CNN layers in full space, vs population size (log axis)');
save_figure(fig, S, out_dir, 'vs_size_cnn_vs_neurons');
crow = add_curves(crow, P, 'N', @(p) pooled(p.cue.summary, 'pair'), 'N', curve_spec);

% 4. vs noise-corrected PR, full space
[fig, ax, h] = curve_figure(S, P, 1, 2, @(p) pooled(p.cue.summary, 'pair'), 'pr_cv', curve_spec(2:3, :), false);
for a = 1:numel(ax), xlabel(ax(a), 'Participation ratio (noise-corrected)'); end
finish(fig, ax, h, P, 'Neurons vs CNN layers in full space, vs dimensionality (each point is one population size)');
save_figure(fig, S, out_dir, 'vs_pr_cnn_vs_neurons');
crow = add_curves(crow, P, 'pr_cv', @(p) pooled(p.cue.summary, 'pair'), 'pr_cv', curve_spec);

save_table(cell2table(crow, 'VariableNames', {'kind', 'population', 'x_name', 'x', 'self', 'pt_minus_ctrl', 'headroom'}), ...
    S, out_dir, 'cnn_vs_neurons_curves');

%% ===========================================================================
function T = pooled(T, col)
% rows of a summary table pooled over the rendering pairs
T = T(strcmp(T.(col), 'pooled'), :);
end

function [fig, ax, h] = curve_figure(S, P, nr, nc, getter, xcol, spec, logx)
fig = figure('Visible', 'off', 'Position', [100 100 1800 560], 'Color', 'w');
tl  = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(1, size(spec, 1));
h  = gobjects(1, numel(P));
for s = 1:size(spec, 1)
    ax(s) = nexttile(tl);  hold(ax(s), 'on');
    for i = 1:numel(P)
        T = getter(P(i));
        h(i) = plot(ax(s), T.(xcol), T.(spec{s, 1}), 'Color', P(i).color, 'LineStyle', P(i).ls, 'LineWidth', 1.6, ...
            'Marker', P(i).mk, 'MarkerSize', 3, 'MarkerFaceColor', P(i).color);
    end
    if ~strcmp(spec{s, 1}, 'self'), yline(ax(s), 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off'); end
    if logx, set(ax(s), 'XScale', 'log'); end
    ylabel(ax(s), spec{s, 2});  box(ax(s), 'on');
end
if logx
    for s = 1:numel(ax), set(ax(s), 'XTick', [1 2 5 10 20 50 100 200 500]); end
end
end

function finish(fig, ax, h, P, ttl)
lg = legend(ax(1), h, {P.name}, 'Orientation', 'horizontal', 'NumColumns', 6);
lg.Layout.Tile = 'south';
title(fig.Children(end), ttl, 'FontSize', 11);
end

function crow = add_curves(crow, P, x_name, getter, xcol, spec) %#ok<INUSD>
for i = 1:numel(P)
    T = getter(P(i));
    for j = 1:height(T)
        crow(end + 1, :) = {P(i).kind, P(i).name, x_name, T.(xcol)(j), T.self(j), T.pt_minus_ctrl(j), T.headroom(j)}; %#ok<AGROW>
    end
end
end
