%{
Filename: plot_normalized_metrics_vs_size.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Normalized transfer metrics (PT / control, headroom, PT - control; pooled
    means, transfer_metrics.m) as a function of population size, full neuron
    space. Rows: the three metrics; columns: KO, FR; V1 blue, V2 orange. Shading
    is the percentile-bootstrap 95% CI over neuron draws (folds when the whole
    population is used once) or over half-splits.

Output (results/figures/05_normalized_metrics_vs_size/):
    cue_metrics_vs_N                    cue transfer vs number of neurons
    halves_metrics_vs_size              within-area halves vs neurons per half
    normalized_metrics_vs_size.csv      every plotted point
%}

clc; clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
R = load_all_results(S);
out_dir = '05_normalized_metrics_vs_size';
fprintf('plot_normalized_metrics_vs_size\n');

rows = [vs_size_figure(R, S, 'cue', out_dir); vs_size_figure(R, S, 'halves', out_dir)];
T = cell2table(rows, 'VariableNames', {'comparison', 'population', 'size', ...
    'self', 'pt', 'ctrl', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', 'headroom', 'headroom_lo', 'headroom_hi', ...
    'pt_minus_ctrl', 'diff_lo', 'diff_hi', 'n_units'});
save_table(T, S, out_dir, 'normalized_metrics_vs_size');

%% ===========================================================================
function rows = vs_size_figure(R, S, kind, out_dir)
spec = {'ratio', 'PT / control', 1; 'headroom', 'Headroom', 0; 'diff', 'PT - control', 0};
fig = figure('Visible', 'off', 'Position', [100 100 1250 900], 'Color', 'w');
tl  = tiledlayout(fig, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
seed = 5000;
rows = {};
for im = 1:2
    m = S.monkeys{im};
    [x, V] = deal(cell(1, 2));
    for a = 1:2
        pop = sprintf('%s_V%d', m, a);
        r = R.(kind).(pop);
        if strcmp(kind, 'cue'), x{a} = r.N_list(:)'; else, x{a} = r.h_list(:)'; end
        V{a} = struct('ratio', [], 'ratio_ci', [], 'headroom', [], 'headroom_ci', [], 'diff', [], 'diff_ci', []);
        for i = 1:numel(x{a})
            if strcmp(kind, 'cue'), D = r.n_draws(i); else, D = size(r.self, 3); end
            ud = 3 + (D == 1);                               % draws, or folds when one draw
            seed = seed + 1;
            B = bootstrap_metrics(r.self(i, :, 1:D, :), r.pt_rot(i, :, 1:D, :), r.ctrl_m(i, :, 1:D, :), ud, S, seed);
            for f = {'ratio', 'headroom', 'diff'}
                V{a}.(f{1})(i) = B.(f{1});
                V{a}.([f{1} '_ci'])(i, :) = B.([f{1} '_ci']);
            end
            rows(end + 1, :) = {kind, strrep(pop, '_', ' '), x{a}(i), B.self, B.pt, B.ctrl, ...
                B.ratio, B.ratio_ci(1), B.ratio_ci(2), B.headroom, B.headroom_ci(1), B.headroom_ci(2), ...
                B.diff, B.diff_ci(1), B.diff_ci(2), B.n_units}; %#ok<AGROW>
        end
    end
    for s = 1:3
        ax = nexttile(tl, (s - 1) * 2 + im);  hold(ax, 'on');
        h = gobjects(1, 2);
        for a = 1:2
            c = S.c.(sprintf('V%d', a));
            ci = V{a}.([spec{s, 1} '_ci']);
            fill(ax, [x{a}, fliplr(x{a})], [ci(:, 1)', fliplr(ci(:, 2)')], c, 'FaceAlpha', 0.15, ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
            h(a) = plot(ax, x{a}, V{a}.(spec{s, 1}), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 3, 'LineWidth', 1.5);
        end
        yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
        box(ax, 'on');  ylabel(ax, spec{s, 2});
        if s == 1, title(ax, m); legend(ax, h, {'V1', 'V2'}, 'Location', 'best'); end
        if s == 3
            if strcmp(kind, 'cue'), xlabel(ax, 'Number of neurons'); else, xlabel(ax, 'Neurons per half'); end
        end
    end
end
if strcmp(kind, 'cue')
    title(tl, 'Cue transfer vs population size, full neuron space (shading: bootstrap 95% CI over neuron draws)', 'FontSize', 11);
    save_figure(fig, S, out_dir, 'cue_metrics_vs_N');
else
    title(tl, 'Within-area halves vs half size, full neuron space (shading: bootstrap 95% CI over half-splits)', 'FontSize', 11);
    save_figure(fig, S, out_dir, 'halves_metrics_vs_size');
end
end
