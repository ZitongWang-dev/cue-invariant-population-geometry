%{
Filename: plot_models_vs_neurons.m
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Gabor models, CNN layers and recorded neurons on the same axes, all from the
    corrected-control reruns (neuron_decoding_rev, ANN_decoding_rev,
    Gabor_decoding_rev):
      - rebuttal point 9 redone: neurons vs the Gabor filter bank vs population
        size and vs the noise-corrected participation ratio (the published point-9
        figures scored the Gabor bank against the old control);
      - the model comparison at k = 5 PCs per side, where the rotation has the
        same 10 free parameters for every population, model and layer;
      - the same along the PC sweep (k = 1-49).
    Neurons: V1 blue, V2 orange (KO solid, FR dashed). Gabor: simple black,
    complex gray (diamonds). CNN layers: AlexNet greens, VGG-19 purples, darker
    = deeper. Model cells and CNN units are nearly noiseless (self-decoding about
    1) while neural self-decoding is 0.1-0.5, so the self-decoding panel is part
    of every comparison.

Output (results/figures/09_models_vs_neurons/):
    point9_vs_cells              self, PT and control, PT - control, headroom vs number of
                                 neurons / cells, full space (neurons and Gabor)
    point9_vs_pr                 PT, PT - control and headroom vs noise-corrected PR (neurons and Gabor)
    k5_models_vs_neurons         bars at k = 5: 4 neural populations, 2 Gabor models, 8 CNN layers
                                 (bootstrap 95% CI over folds)
    pc_sweep_models_vs_neurons   self, PT - control and headroom vs k (all three families)
    models_vs_neurons_k5.csv     the k = 5 bars
    models_vs_neurons_curves.csv every plotted curve point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S  = gabor_figure_settings();
R  = load_gabor_results(S);
Sn = figure_settings();             % neuron figures (neuron_decoding_rev helpers)
Rn = load_all_results(Sn);
Sa = ann_figure_settings();         % CNN figures (ANN_decoding_rev helpers)
Ra = load_ann_results(Sa);
out_dir = '09_models_vs_neurons';
fprintf('plot_models_vs_neurons\n');

%% the populations: name, family, results, color, line style, marker, bar transparency
P = struct('name', {}, 'kind', {}, 'pc', {}, 'cue', {}, 'color', {}, 'ls', {}, 'mk', {}, 'alpha', {});
for p = {'KO_V1', 'KO_V2', 'FR_V1', 'FR_V2'}
    pop = p{1};  isKO = startsWith(pop, 'KO');
    if isKO, ls = '-'; a = 1; else, ls = '--'; a = 0.45; end
    P(end + 1) = struct('name', strrep(pop, '_', ' '), 'kind', 'neurons', 'pc', Rn.pc_cue.(pop), 'cue', Rn.cue.(pop), ...
        'color', Sn.c.(pop(4:5)), 'ls', ls, 'mk', 'o', 'alpha', a); %#ok<SAGROW>
end
gab_col = struct('simple', [0.12 0.12 0.12], 'complex', [0.55 0.55 0.55]);
for j = 1:2
    m = S.models{j};
    P(end + 1) = struct('name', S.label.(m), 'kind', 'Gabor', 'pc', R.pc.(m), 'cue', R.cue.(m), ...
        'color', gab_col.(m), 'ls', '-', 'mk', 'd', 'alpha', 1); %#ok<SAGROW>
end
shades = struct('alexnet', [0.70 0.86 0.55; 0.35 0.65 0.25; 0.08 0.38 0.08], ...
                'vgg',     [0.86 0.78 0.94; 0.72 0.58 0.86; 0.57 0.38 0.76; 0.42 0.20 0.62; 0.27 0.05 0.46]);
for in = 1:2
    net = Sa.nets{in};
    idx = find(strcmp(Sa.layers(:, 1), net))';
    for j = 1:numel(idx)
        key = Sa.keys{idx(j)};
        P(end + 1) = struct('name', Ra.cue.(key).population, 'kind', 'CNN', 'pc', Ra.pc_cue.(key), 'cue', Ra.cue.(key), ...
            'color', shades.(net)(j, :), 'ls', Sa.ls.net.(net), 'mk', Sa.mk.net.(net), 'alpha', 1); %#ok<SAGROW>
    end
end
isNG = ismember({P.kind}, {'neurons', 'Gabor'});       % rebuttal point 9: neurons and the Gabor bank
crow = {};

%% 1. point 9: vs number of neurons / cells, full space
Q = P(isNG);
fig = figure('Visible', 'off', 'Position', [100 100 1900 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = gobjects(1, 4);  h = gobjects(1, numel(Q));
for s = 1:4, ax(s) = nexttile(tl);  hold(ax(s), 'on'); end
for i = 1:numel(Q)
    T = pooled(Q(i).cue.summary, 'pair');
    st = {'Color', Q(i).color, 'LineStyle', Q(i).ls, 'LineWidth', 1.6, 'Marker', Q(i).mk, 'MarkerSize', 3, 'MarkerFaceColor', Q(i).color};
    h(i) = plot(ax(1), T.N, T.self, st{:});
    plot(ax(2), T.N, T.pt_rot, st{:});
    plot(ax(2), T.N, T.ctrl, 'Color', Q(i).color, 'LineStyle', ':', 'LineWidth', 1.4, 'HandleVisibility', 'off');
    plot(ax(3), T.N, T.pt_minus_ctrl, st{:});
    plot(ax(4), T.N, T.headroom, st{:});
    crow = add_curve_rows(crow, Q(i), 'N', T, T.N);
end
yl = {'Self-decoding', 'PT (lines) and control (dotted)', 'PT - control (accuracy)', 'Headroom  (PT - ctrl) / (self - ctrl)'};
for s = 1:4
    if s > 2, yline(ax(s), 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off'); end
    xlim(ax(s), [0 200]);  box(ax(s), 'on');  xlabel(ax(s), 'Number of neurons / model cells');  ylabel(ax(s), yl{s});
end
lg = legend(ax(1), h, {Q.name}, 'Orientation', 'horizontal', 'NumColumns', 6);  lg.Layout.Tile = 'south';
title(tl, 'Rebuttal point 9 with the corrected control: neurons vs Gabor models vs population size, full space (pooled over pairs)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'point9_vs_cells');

%% 2. point 9: vs noise-corrected participation ratio, full space
fig = figure('Visible', 'off', 'Position', [100 100 1700 520], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {'pt_rot', 'PT (rotation only)'; 'pt_minus_ctrl', 'PT - control (accuracy)'; 'headroom', 'Headroom'};
ax = gobjects(1, 3);
for s = 1:3
    ax(s) = nexttile(tl);  hold(ax(s), 'on');
    for i = 1:numel(Q)
        T = pooled(Q(i).cue.summary, 'pair');
        h(i) = plot(ax(s), T.pr_cv, T.(spec{s, 1}), 'Color', Q(i).color, 'LineStyle', Q(i).ls, 'LineWidth', 1.6, ...
            'Marker', Q(i).mk, 'MarkerSize', 3, 'MarkerFaceColor', Q(i).color);
    end
    if s > 1, yline(ax(s), 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off'); end
    box(ax(s), 'on');  xlabel(ax(s), 'Participation ratio (noise-corrected)');  ylabel(ax(s), spec{s, 2});
end
lg = legend(ax(1), h, {Q.name}, 'Orientation', 'horizontal', 'NumColumns', 6);  lg.Layout.Tile = 'south';
title(tl, 'Rebuttal point 9 with the corrected control: neurons vs Gabor models vs dimensionality (each point is one population size)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'point9_vs_pr');
for i = 1:numel(Q)
    T = pooled(Q(i).cue.summary, 'pair');
    crow = add_curve_rows(crow, Q(i), 'pr_cv', T, T.pr_cv);
end

%% 3. bars at k = 5 (whole populations / models / layers)
seed = 45000;
nP = numel(P);
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
    S, out_dir, 'models_vs_neurons_k5');

spec = {'self', 'Self-decoding', NaN; 'diff', 'PT - control (accuracy)', 0; 'headroom', 'Headroom  (PT - ctrl) / (self - ctrl)', 0; ...
        'ratio', 'PT / control', 1};
x = [1:4, 5.6:6.6, 8.2:10.2, 11.8:15.8];        % gaps between neurons, Gabor, AlexNet and VGG
fig = figure('Visible', 'off', 'Position', [100 100 1950 580], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
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
    set(ax, 'XTick', x, 'XTickLabel', {P.name}, 'XTickLabelRotation', 45, 'FontSize', 8);
    xlim(ax, [0.3, max(x) + 0.7]);  ylabel(ax, spec{s, 2});  box(ax, 'on');
end
title(tl, sprintf(['Neurons vs Gabor models vs CNN layers at k = %d PCs per side (same rotation freedom), whole populations; ' ...
    'error bars: bootstrap 95%% CI over folds'], S.k_main), 'FontSize', 11);
save_figure(fig, S, out_dir, 'k5_models_vs_neurons');

%% 4. along the PC sweep (all three families)
fig = figure('Visible', 'off', 'Position', [100 100 1800 600], 'Color', 'w');
tl  = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
spec = {'self', 'Self-decoding'; 'pt_minus_ctrl', 'PT - control (accuracy)'; 'headroom', 'Headroom'};
ax = gobjects(1, 3);  h = gobjects(1, nP);
for s = 1:3
    ax(s) = nexttile(tl);  hold(ax(s), 'on');
    for i = 1:nP
        T = pooled(P(i).pc.summary, 'set');
        lw = 1.2;  if ~strcmp(P(i).kind, 'CNN'), lw = 2; end
        h(i) = plot(ax(s), T.k, T.(spec{s, 1}), 'Color', P(i).color, 'LineStyle', P(i).ls, 'LineWidth', lw, ...
            'Marker', P(i).mk, 'MarkerSize', 3, 'MarkerFaceColor', P(i).color);
    end
    if s > 1, yline(ax(s), 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off'); end
    xline(ax(s), S.k_main, ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    set(ax(s), 'XScale', 'log', 'XTick', [1 2 5 10 20 49]);  xlim(ax(s), [0.9 55]);
    box(ax(s), 'on');  xlabel(ax(s), 'Retained PCs per side, k');  ylabel(ax(s), spec{s, 2});
end
lg = legend(ax(1), h, {P.name}, 'Orientation', 'horizontal', 'NumColumns', 7);  lg.Layout.Tile = 'south';
title(tl, 'Neurons vs Gabor models vs CNN layers along the PC sweep (whole populations; pooled over pairs; dotted line: k = 5)', 'FontSize', 11);
save_figure(fig, S, out_dir, 'pc_sweep_models_vs_neurons');
for i = 1:nP
    T = pooled(P(i).pc.summary, 'set');
    crow = add_curve_rows(crow, P(i), 'k', T, T.k);
end

save_table(cell2table(crow, 'VariableNames', {'kind', 'population', 'x_name', 'x', 'self', 'pt', 'ctrl', 'pt_minus_ctrl', 'headroom'}), ...
    S, out_dir, 'models_vs_neurons_curves');

%% ===========================================================================
function T = pooled(T, col)
% rows of a summary table pooled over the rendering pairs
T = T(strcmp(T.(col), 'pooled'), :);
end

function crow = add_curve_rows(crow, p, x_name, T, x)
for j = 1:height(T)
    crow(end + 1, :) = {p.kind, p.name, x_name, x(j), T.self(j), T.pt_rot(j), T.ctrl(j), T.pt_minus_ctrl(j), T.headroom(j)}; %#ok<AGROW>
end
end
