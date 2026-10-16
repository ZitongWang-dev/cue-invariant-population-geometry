%{
Filename: plot_normalized_metrics_vs_size.m
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Normalized transfer metrics of the CNN layers (pooled means,
    transfer_metrics.m) as a function of the number of units, full unit space:
    the test of the paper's "plateau at 250-300 units". Rows: PT / control,
    PT - control, headroom, and the headroom denominator self - control;
    columns: AlexNet, VGG-19; one line per layer (colors as published). Shading
    is the percentile-bootstrap 95% CI over unit draws (folds for the whole
    layer, the last point of each line). The dotted line marks the Fig 3d size.

Output (results/figures/05_normalized_metrics_vs_size/):
    cue_metrics_vs_units               the figure described above
    normalized_metrics_vs_size.csv     every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = ann_figure_settings();
R = load_ann_results(S);
out_dir = '05_normalized_metrics_vs_size';
fprintf('plot_normalized_metrics_vs_size\n');

spec = {'ratio', 'PT / control', 1; 'diff', 'PT - control', 0; 'headroom', 'Headroom', 0; 'denom', 'self - control', 0};
seed = 6000;
rows = {};
V = struct();
for i = 1:numel(S.keys)
    key = S.keys{i};  r = R.cue.(key);
    x = r.N_list(:)';  nN = numel(x);
    V.(key) = struct('x', x, 'ratio', nan(1, nN), 'ratio_ci', nan(nN, 2), 'diff', nan(1, nN), 'diff_ci', nan(nN, 2), ...
        'headroom', nan(1, nN), 'headroom_ci', nan(nN, 2), 'denom', nan(1, nN), 'denom_ci', nan(nN, 2));
    for iN = 1:nN
        D = r.n_draws(iN);
        ud = 3 + (D == 1);                               % draws, or folds when one draw
        seed = seed + 1;
        B = bootstrap_metrics(r.self(iN, :, 1:D, :), r.pt_rot(iN, :, 1:D, :), r.ctrl_m(iN, :, 1:D, :), ud, S, seed);
        for f = {'ratio', 'diff', 'headroom'}
            V.(key).(f{1})(iN) = B.(f{1});
            V.(key).([f{1} '_ci'])(iN, :) = B.([f{1} '_ci']);
        end
        V.(key).denom(iN) = B.self - B.ctrl;
        rows(end + 1, :) = {r.population, x(iN), B.self, B.pt, B.ctrl, B.ratio, B.ratio_ci(1), B.ratio_ci(2), ...
            B.diff, B.diff_ci(1), B.diff_ci(2), B.headroom, B.headroom_ci(1), B.headroom_ci(2), B.self - B.ctrl, B.n_units}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'layer', 'N', 'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'headroom', 'headroom_lo', 'headroom_hi', 'self_minus_ctrl', 'n_units_boot'});
save_table(T, S, out_dir, 'normalized_metrics_vs_size');

%% figure
fig = figure('Visible', 'off', 'Position', [100 100 1300 1150], 'Color', 'w');
tl  = tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for in = 1:2
    net  = S.nets{in};
    idx  = find(strcmp(S.layers(:, 1), net))';
    n_max = max(arrayfun(@(i) R.cue.(S.keys{i}).N_total, idx));
    for s = 1:size(spec, 1)
        ax = nexttile(tl, (s - 1) * 2 + in);  hold(ax, 'on');
        h = gobjects(1, numel(idx));
        for j = 1:numel(idx)
            key = S.keys{idx(j)};  c = S.c.layer.(key);  v = V.(key);
            ci = v.([spec{s, 1} '_ci']);
            if ~all(isnan(ci(:)))
                fill(ax, [v.x, fliplr(v.x)], [ci(:, 1)', fliplr(ci(:, 2)')], c, 'FaceAlpha', 0.15, ...
                    'EdgeColor', 'none', 'HandleVisibility', 'off');
            end
            h(j) = plot(ax, v.x, v.(spec{s, 1}), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3, 'LineWidth', 1.5);
        end
        yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        xline(ax, S.fig3d_N.(net), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
        xlim(ax, [0, n_max * 1.02]);  box(ax, 'on');  ylabel(ax, spec{s, 2});
        if s == 1
            title(ax, S.net_label.(net));
            legend(ax, h, arrayfun(@(i) sprintf('pool %d', S.layers{i, 2}), idx, 'UniformOutput', false), 'Location', 'best');
        end
        if s == size(spec, 1), xlabel(ax, 'Number of units'); end
    end
end
title(tl, {'CNN cue transfer vs number of units, full unit space (shading: bootstrap 95% CI over unit draws; dotted line: Fig 3d size)', ...
    'headroom = (PT - control) / (self - control); its denominator (bottom row) shrinks as the control approaches self-decoding'}, 'FontSize', 11);
save_figure(fig, S, out_dir, 'cue_metrics_vs_units');
