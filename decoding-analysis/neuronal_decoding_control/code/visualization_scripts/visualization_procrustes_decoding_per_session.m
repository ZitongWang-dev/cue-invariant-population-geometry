%{
Filename: visualization_procrustes_decoding_per_session.m
Author: Zitong Wang
Date: 2026-04-30

Description:
    Visualizes per-session Procrustes decoding results produced by
    procrustes_decoding_per_session.m.

    Six figures are produced:
      Figure 1 - perturb_mode = 'none'   (real noise correlations preserved)
      Figure 2 - perturb_mode = 'affine' (within-session correlations broken)
      Figure 3 - delta (affine - none) per session per decoding type
      Figure 4 - self-decoding vs PT-decoding scatter, none + affine
                 overlaid; affine markers carry a black edge; lines
                 connect the same session's none and affine points.
      Figure 5 - normalized delta = 2*(affine-none)/(affine+none) per
                 session, with Wilcoxon signed-rank tests on the
                 across-session median, FDR-corrected across all
                 (decoding type x stim pair x area) tests.
      Figure 6 - pooled across stim pairs and across both monkeys:
                 per (monkey, area, decoding type) the N_sessions x 3
                 stim-pair normalized deltas are pooled into a single
                 sample, Wilcoxon-tested against H0: median = 0, and
                 FDR-corrected across all 12 tests. Independent of the
                 top-level monkey setting (loads both FR and KO).

    Layout (Figures 1-3):
      - 3 subplots side-by-side: EC-AC, EC-EX, AC-EX
      - x-axis: 5 decoding types {self, non-PT, PT(rot), PT ctrl, non-PT ctrl}
      - V1 sessions jittered LEFT of each x-tick (blue, graded shades)
      - V2 sessions jittered RIGHT of each x-tick (orange, graded shades)
      - Lines connect points from the same session across decoding types,
        but never across stimulus pairs (subplot boundary enforces this)

    Layout (Figure 4):
      - 3 subplots side-by-side: EC-AC, EC-EX, AC-EX
      - x-axis: self-decoding accuracy
      - y-axis: PT-decoding accuracy (rotation only)
      - Each session contributes 2 points (none + affine) connected by a
        session-colored line. Affine markers carry a black edge.
      - Diagonal y=x dashed reference line.

    Layout (Figure 5):
      - 3 subplots side-by-side: EC-AC, EC-EX, AC-EX
      - x-axis: 3 decoding types {self, non-PT, PT(rot)} (controls dropped)
      - y-axis: normalized delta
      - Significance markers above each (decoding type, area) cluster.

    Layout (Figure 6):
      - 2 subplots side-by-side: FR, KO
      - x-axis: 3 decoding types (self, non-PT, PT (rot))
      - V1 dodged LEFT of each x-tick, V2 dodged RIGHT
      - Each cluster: N_sessions x 3 scatter points
          - Marker shape: stim pair  (o = EC-AC, s = EC-EX, ^ = AC-EX)
          - Color shade : session    (light = lowest index, dark = highest)
      - Cluster overlays: thick dark median bar, thin black mean tick
      - Asterisks above each cluster from FDR-corrected q.
%}

clc; clear;
close all
%% Configuration
monkey        = 'KO';
file_location = ['..\..\results\decoding_outputs\','procrustes_decoding_per_session_results'];
fig_save_dir  = ['..\..\results\figures\','procrustes_decoding_per_session_', monkey];
% if ~exist(fig_save_dir, 'dir'), mkdir(fig_save_dir); end

%% Decoding type columns to display
% Original 8-column accuracy order (from pro_decoding):
%   1=genAcc, 2=non-PT, 3=full PT, 4=rand PT ctrl, 5=scale-only,
%   6=rotation-only, 7=translation-only, 8=non-PT ctrl
dec_idx          = [1, 2, 6, 4, 8];
decoding_labels  = {'self', 'non-PT', 'PT (rot)', 'PT ctrl', 'non-PT ctrl'};

%% Load per-session data for both perturb modes
[none_v1,   none_v1_meta,   none_v2,   none_v2_meta]   = ...
    load_per_session_data(monkey, 'none',   file_location, dec_idx);
[affine_v1, affine_v1_meta, affine_v2, affine_v2_meta] = ...
    load_per_session_data(monkey, 'affine', file_location, dec_idx);

assert(none_v1_meta.n_sessions == affine_v1_meta.n_sessions, ...
    'V1 session counts differ between none and affine.');
assert(none_v2_meta.n_sessions == affine_v2_meta.n_sessions, ...
    'V2 session counts differ between none and affine.');

%% Pack data structs (used by all plotting functions below)
data_none   = pack_data(none_v1,   none_v1_meta,   none_v2,   none_v2_meta);
data_affine = pack_data(affine_v1, affine_v1_meta, affine_v2, affine_v2_meta);

%% Figure 1: perturb_mode = none
plot_per_session(data_none, monkey, 'none (correlations preserved)', false, decoding_labels);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_none.png']));

%% Figure 2: perturb_mode = affine
plot_per_session(data_affine, monkey, 'affine (within-session correlations broken)', false, decoding_labels);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_affine.png']));

%% Figure 3: delta = affine - none
delta_v1   = cellfun(@(a, b) a - b, affine_v1, none_v1, 'UniformOutput', false);
delta_v2   = cellfun(@(a, b) a - b, affine_v2, none_v2, 'UniformOutput', false);
data_delta = pack_data(delta_v1, none_v1_meta, delta_v2, none_v2_meta);
plot_per_session(data_delta, monkey, '\Delta (affine - none)', true, decoding_labels);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_delta.png']));

%% Figure 4: self vs PT-decoding scatter (none + affine on the same axes)
plot_self_vs_pt(data_none, data_affine, monkey);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_self_vs_pt.png']));

%% Figure 5: normalized delta per session, with significance annotations
plot_normalized_delta(data_none, data_affine, monkey);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_norm_delta.png']));

%% Figure 6: pooled across stim pairs AND both monkeys (FR + KO)
% Independent of the top-level `monkey` setting: this section always
% loads FR and KO internally and runs one Wilcoxon test per
% (monkey, area, decoding type) on the pooled N_sessions x 3 deltas.
plot_pooled_pairs(file_location, dec_idx);
% saveas(gcf, fullfile('..\..\results\figures\', ...
%     'procrustes_decoding_per_session_pooled_pairs.png'));

%% =================== LOCAL FUNCTIONS ===================
function data = pack_data(v1_data, v1_meta, v2_data, v2_meta)
    data.v1_data = v1_data;
    data.v1_meta = v1_meta;
    data.v2_data = v2_data;
    data.v2_meta = v2_meta;
end

function [v1_data, v1_meta, v2_data, v2_meta] = ...
    load_per_session_data(monkey, perturb_mode, file_location, dec_idx)
    pair_names = {'acec', 'ecac', 'ecex', 'exec', 'acex', 'exac'};

    v1_pair_data = struct();
    for k = 1:numel(pair_names)
        pn    = pair_names{k};
        fname = fullfile(pwd, file_location, perturb_mode, monkey, 'V1', [pn '_results.mat']);
        S     = load(fname);
        v1_pair_data.(pn) = S.([pn '_results']);
    end
    n_v1     = numel(v1_pair_data.acec);
    v1_files = {v1_pair_data.acec.session_file};

    v1_data    = cell(1, 3);
    v1_data{1} = combine_per_session(v1_pair_data.acec, v1_pair_data.ecac, dec_idx);  % EC-AC
    v1_data{2} = combine_per_session(v1_pair_data.ecex, v1_pair_data.exec, dec_idx);  % EC-EX
    v1_data{3} = combine_per_session(v1_pair_data.acex, v1_pair_data.exac, dec_idx);  % AC-EX
    v1_meta    = struct('n_sessions', n_v1, 'session_files', {v1_files});

    v2_pair_data = struct();
    for k = 1:numel(pair_names)
        pn    = pair_names{k};
        fname = fullfile(pwd, file_location, perturb_mode, monkey, 'V2', [pn '_results.mat']);
        S     = load(fname);
        v2_pair_data.(pn) = S.([pn '_results']);
    end
    n_v2     = numel(v2_pair_data.acec);
    v2_files = {v2_pair_data.acec.session_file};

    v2_data    = cell(1, 3);
    v2_data{1} = combine_per_session(v2_pair_data.acec, v2_pair_data.ecac, dec_idx);
    v2_data{2} = combine_per_session(v2_pair_data.ecex, v2_pair_data.exec, dec_idx);
    v2_data{3} = combine_per_session(v2_pair_data.acex, v2_pair_data.exac, dec_idx);
    v2_meta    = struct('n_sessions', n_v2, 'session_files', {v2_files});
end

function combined = combine_per_session(pair1_struct, pair2_struct, dec_idx)
    n        = numel(pair1_struct);
    combined = zeros(n, numel(dec_idx));
    for s = 1:n
        stacked       = [pair1_struct(s).accuracy; pair2_struct(s).accuracy];
        col_mean      = mean(stacked, 1);
        combined(s,:) = col_mean(dec_idx);
    end
end

function plot_per_session(data, monkey, title_str, is_delta, decoding_labels)
    pair_labels = {'EC-AC', 'EC-EX', 'AC-EX'};
    n_dec       = numel(decoding_labels);
    x_base      = 1:n_dec;

    v1_base = [0      0.4470 0.7410];
    v2_base = [0.8500 0.3250 0.0980];

    n_v1 = data.v1_meta.n_sessions;
    n_v2 = data.v2_meta.n_sessions;

    v1_colors = make_shades(v1_base, n_v1);
    v2_colors = make_shades(v2_base, n_v2);

    if n_v1 > 1, v1_jitter = linspace(-0.30, -0.05, n_v1); else, v1_jitter = -0.175; end
    if n_v2 > 1, v2_jitter = linspace( 0.05,  0.30, n_v2); else, v2_jitter =  0.175; end

    figure('Position', [100 100 1300 480]);
    ax_handles = gobjects(1, 3);

    for p = 1:3
        ax_handles(p) = subplot(1, 3, p);
        hold on;

        v1_pair = data.v1_data{p};
        v2_pair = data.v2_data{p};

        for s = 1:size(v1_pair, 1)
            plot(x_base + v1_jitter(s), v1_pair(s, :), '-o', ...
                 'Color', v1_colors(s, :), 'MarkerFaceColor', v1_colors(s, :), ...
                 'MarkerEdgeColor', v1_colors(s, :), 'MarkerSize', 5, ...
                 'LineWidth', 1.0, 'HandleVisibility','off');
        end
        for s = 1:size(v2_pair, 1)
            plot(x_base + v2_jitter(s), v2_pair(s, :), '-o', ...
                 'Color', v2_colors(s, :), 'MarkerFaceColor', v2_colors(s, :), ...
                 'MarkerEdgeColor', v2_colors(s, :), 'MarkerSize', 5, ...
                 'LineWidth', 1.0, 'HandleVisibility','off');
        end

        if is_delta
            yline(0,    '--k', 'LineWidth', 1, 'HandleVisibility','off');
        else
            yline(0.02, '--k', 'LineWidth', 1, 'HandleVisibility','off');
        end

        for x = x_base
            xline(x, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');
        end

        xlim([0.4, n_dec + 0.6]);
        xticks(x_base);
        xticklabels(decoding_labels);
        xtickangle(20);
        title(pair_labels{p});

        if p == 1
            if is_delta, ylabel('\Delta decoding accuracy');
            else,        ylabel('Decoding accuracy'); end
        end

        grid on; box on;

        if p == 3
            hl_v1 = plot(NaN, NaN, '-o', 'Color', v1_base, 'MarkerFaceColor', v1_base, ...
                'MarkerSize', 5, 'LineWidth', 1.0);
            hl_v2 = plot(NaN, NaN, '-o', 'Color', v2_base, 'MarkerFaceColor', v2_base, ...
                'MarkerSize', 5, 'LineWidth', 1.0);
            hl_chance = plot(NaN, NaN, '--k', 'LineWidth', 1);
            if is_delta
                legend([hl_v1, hl_v2, hl_chance], ...
                    {sprintf('V1 (n=%d sessions)', n_v1), ...
                     sprintf('V2 (n=%d sessions)', n_v2), 'no effect'}, ...
                    'Location', 'best');
            else
                legend([hl_v1, hl_v2, hl_chance], ...
                    {sprintf('V1 (n=%d sessions)', n_v1), ...
                     sprintf('V2 (n=%d sessions)', n_v2), 'chance (2%)'}, ...
                    'Location', 'best');
            end
        end
        hold off;
    end

    linkaxes(ax_handles, 'y');
    sgtitle(sprintf('%s | per-session Procrustes decoding | %s', monkey, title_str));
end

function plot_self_vs_pt(data_none, data_affine, monkey)
    self_col = 1;
    pt_col   = 3;

    pair_labels = {'EC-AC', 'EC-EX', 'AC-EX'};

    v1_base = [0      0.4470 0.7410];
    v2_base = [0.8500 0.3250 0.0980];

    n_v1 = data_none.v1_meta.n_sessions;
    n_v2 = data_none.v2_meta.n_sessions;

    v1_colors = make_shades(v1_base, n_v1);
    v2_colors = make_shades(v2_base, n_v2);

    figure('Position', [100 100 1300 480]);
    ax_handles = gobjects(1, 3);

    global_max = 0;
    for p = 1:3
        ax_handles(p) = subplot(1, 3, p);
        hold on;

        for s = 1:n_v1
            x_n = data_none.v1_data{p}(s,   self_col);
            y_n = data_none.v1_data{p}(s,   pt_col);
            x_a = data_affine.v1_data{p}(s, self_col);
            y_a = data_affine.v1_data{p}(s, pt_col);
            plot([x_n, x_a], [y_n, y_a], '-', 'Color', v1_colors(s,:), ...
                 'LineWidth', 0.8, 'HandleVisibility','off');
            plot(x_n, y_n, 'o', 'MarkerFaceColor', v1_colors(s,:), ...
                 'MarkerEdgeColor', v1_colors(s,:), 'MarkerSize', 7, 'HandleVisibility','off');
            plot(x_a, y_a, 'o', 'MarkerFaceColor', v1_colors(s,:), ...
                 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility','off');
            global_max = max([global_max, x_n, y_n, x_a, y_a]);
        end
        for s = 1:n_v2
            x_n = data_none.v2_data{p}(s,   self_col);
            y_n = data_none.v2_data{p}(s,   pt_col);
            x_a = data_affine.v2_data{p}(s, self_col);
            y_a = data_affine.v2_data{p}(s, pt_col);
            plot([x_n, x_a], [y_n, y_a], '-', 'Color', v2_colors(s,:), ...
                 'LineWidth', 0.8, 'HandleVisibility','off');
            plot(x_n, y_n, 'o', 'MarkerFaceColor', v2_colors(s,:), ...
                 'MarkerEdgeColor', v2_colors(s,:), 'MarkerSize', 7, 'HandleVisibility','off');
            plot(x_a, y_a, 'o', 'MarkerFaceColor', v2_colors(s,:), ...
                 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility','off');
            global_max = max([global_max, x_n, y_n, x_a, y_a]);
        end
    end

    axis_max = ceil(global_max * 20) / 20;
    if axis_max == 0, axis_max = 0.1; end

    for p = 1:3
        subplot(1, 3, p); hold on;
        plot([0, axis_max], [0, axis_max], '--', 'Color', [0.7 0.7 0.7], ...
             'LineWidth', 1, 'HandleVisibility','off');
        xlim([0, axis_max]); ylim([0, axis_max]); axis square;
        xlabel('Self-decoding accuracy');
        if p == 1, ylabel('PT-decoding accuracy (rotation only)'); end
        title(pair_labels{p}); grid on; box on;

        if p == 3
            hl_v1 = plot(NaN, NaN, 'o', 'MarkerFaceColor', v1_base, ...
                'MarkerEdgeColor', v1_base, 'MarkerSize', 7);
            hl_v2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', v2_base, ...
                'MarkerEdgeColor', v2_base, 'MarkerSize', 7);
            hl_none = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], ...
                'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 7);
            hl_aff  = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.5 0.5 0.5], ...
                'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'LineWidth', 1.2);
            hl_diag = plot(NaN, NaN, '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
            legend([hl_v1, hl_v2, hl_none, hl_aff, hl_diag], ...
                {sprintf('V1 (n=%d)', n_v1), sprintf('V2 (n=%d)', n_v2), ...
                 'none (preserved)', 'affine (broken)', 'y = x'}, ...
                'Location', 'best');
        end
        hold off;
    end

    linkaxes(ax_handles, 'xy');
    sgtitle(sprintf('%s | self vs PT decoding | none + affine per session', monkey));
end

function plot_normalized_delta(data_none, data_affine, monkey)
% Figure 5: per-session normalized delta = 2*(affine - none) / (affine + none),
% with one-sample Wilcoxon signed-rank tests across sessions per (area, pair,
% decoding type), FDR-corrected across all 18 tests.

    dec_col_keep = [1, 2, 3];                       % indices in the 5-col stored data
    dec_labels   = {'self', 'non-PT', 'PT (rot)'};
    pair_labels  = {'EC-AC', 'EC-EX', 'AC-EX'};
    n_dec        = numel(dec_labels);
    n_pairs      = 3;

    n_v1 = data_none.v1_meta.n_sessions;
    n_v2 = data_none.v2_meta.n_sessions;

    %% Compute per-session normalized delta per stim pair
    v1_norm = cell(1, n_pairs);
    v2_norm = cell(1, n_pairs);
    for p = 1:n_pairs
        a_v1 = data_affine.v1_data{p}(:, dec_col_keep);
        b_v1 = data_none.v1_data{p}(:,   dec_col_keep);
        v1_norm{p} = 2 .* (a_v1 - b_v1) ./ (a_v1 + b_v1);

        a_v2 = data_affine.v2_data{p}(:, dec_col_keep);
        b_v2 = data_none.v2_data{p}(:,   dec_col_keep);
        v2_norm{p} = 2 .* (a_v2 - b_v2) ./ (a_v2 + b_v2);
    end

    %% One-sample Wilcoxon signed-rank tests
    p_v1 = nan(n_pairs, n_dec);
    p_v2 = nan(n_pairs, n_dec);
    for p = 1:n_pairs
        for d = 1:n_dec
            x_v1 = v1_norm{p}(:, d);
            x_v2 = v2_norm{p}(:, d);
            x_v1 = x_v1(~isnan(x_v1));
            x_v2 = x_v2(~isnan(x_v2));
            if numel(x_v1) >= 1, p_v1(p, d) = signrank(x_v1); end
            if numel(x_v2) >= 1, p_v2(p, d) = signrank(x_v2); end
        end
    end

    %% FDR correction across all 18 tests
    all_p = [p_v1(:); p_v2(:)];
    all_q = bh_fdr(all_p);
    q_v1  = reshape(all_q(1:numel(p_v1)),    size(p_v1));
    q_v2  = reshape(all_q(numel(p_v1)+1:end), size(p_v2));

    %% Print test results to console
    fprintf('\n=== Wilcoxon signed-rank tests on normalized delta (H0: median = 0) ===\n');
    fprintf('FDR-corrected (Benjamini-Hochberg) across %d tests\n\n', numel(all_p));
    fprintf('%-7s | %-9s | V1 (n=%d)   raw_p     q       sig | V2 (n=%d)   raw_p     q       sig\n', ...
        'Pair', 'Decoding', n_v1, n_v2);
    fprintf('%s\n', repmat('-', 1, 105));
    for p = 1:n_pairs
        for d = 1:n_dec
            fprintf('%-7s | %-9s |             %.4f   %.4f   %-3s |             %.4f   %.4f   %s\n', ...
                pair_labels{p}, dec_labels{d}, ...
                p_v1(p,d), q_v1(p,d), sig_marker(q_v1(p,d)), ...
                p_v2(p,d), q_v2(p,d), sig_marker(q_v2(p,d)));
        end
    end
    fprintf('\n');

    %% Plotting parameters
    v1_base = [0      0.4470 0.7410];
    v2_base = [0.8500 0.3250 0.0980];

    v1_colors = make_shades(v1_base, n_v1);
    v2_colors = make_shades(v2_base, n_v2);

    if n_v1 > 1, v1_jitter = linspace(-0.30, -0.05, n_v1); else, v1_jitter = -0.175; end
    if n_v2 > 1, v2_jitter = linspace( 0.05,  0.30, n_v2); else, v2_jitter =  0.175; end

    x_base = 1:n_dec;

    %% Symmetric y-limit
    all_vals = [];
    for p = 1:n_pairs
        all_vals = [all_vals; v1_norm{p}(:); v2_norm{p}(:)]; %#ok<AGROW>
    end
    all_vals = all_vals(~isnan(all_vals) & isfinite(all_vals));
    y_max    = max(0.5, max(abs(all_vals)) * 1.25);

    %% Plot
    figure('Position', [100 100 1150 480]);
    ax_handles = gobjects(1, n_pairs);

    for p = 1:n_pairs
        ax_handles(p) = subplot(1, n_pairs, p);
        hold on;

        for s = 1:n_v1
            plot(x_base + v1_jitter(s), v1_norm{p}(s, :), '-o', ...
                 'Color', v1_colors(s, :), 'MarkerFaceColor', v1_colors(s, :), ...
                 'MarkerEdgeColor', v1_colors(s, :), 'MarkerSize', 5, ...
                 'LineWidth', 1.0, 'HandleVisibility','off');
        end
        for s = 1:n_v2
            plot(x_base + v2_jitter(s), v2_norm{p}(s, :), '-o', ...
                 'Color', v2_colors(s, :), 'MarkerFaceColor', v2_colors(s, :), ...
                 'MarkerEdgeColor', v2_colors(s, :), 'MarkerSize', 5, ...
                 'LineWidth', 1.0, 'HandleVisibility','off');
        end

        yline(0, '--k', 'LineWidth', 1, 'HandleVisibility','off');
        for x = x_base
            xline(x, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');
        end

        xlim([0.4, n_dec + 0.6]);
        ylim([-y_max, y_max]);
        xticks(x_base);
        xticklabels(dec_labels);
        xtickangle(20);
        title(pair_labels{p});
        if p == 1
            ylabel('Normalized \Delta: 2(affine - none)/(affine + none)');
        end
        grid on; box on;

        sig_y = y_max * 0.92;
        for d = 1:n_dec
            text(d - 0.175, sig_y, sig_marker(q_v1(p, d)), ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', v1_base);
            text(d + 0.175, sig_y, sig_marker(q_v2(p, d)), ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', v2_base);
        end

        if p == n_pairs
            hl_v1 = plot(NaN, NaN, '-o', 'Color', v1_base, ...
                'MarkerFaceColor', v1_base, 'MarkerSize', 5, 'LineWidth', 1.0);
            hl_v2 = plot(NaN, NaN, '-o', 'Color', v2_base, ...
                'MarkerFaceColor', v2_base, 'MarkerSize', 5, 'LineWidth', 1.0);
            hl_zero = plot(NaN, NaN, '--k', 'LineWidth', 1);
            legend([hl_v1, hl_v2, hl_zero], ...
                {sprintf('V1 (n=%d sessions)', n_v1), ...
                 sprintf('V2 (n=%d sessions)', n_v2), ...
                 'no effect'}, ...
                'Location', 'best');
        end
        hold off;
    end

    linkaxes(ax_handles, 'y');
    sgtitle(sprintf(['%s | normalized \\Delta per session  |  ' ...
                     'Wilcoxon signed-rank (FDR-corrected):  * q<0.05,  ** q<0.01,  *** q<0.001'], ...
                    monkey));
end

function plot_pooled_pairs(file_location, dec_idx)
% Figure 6: pooled across stim pairs AND across both monkeys.
% Independent of the top-level `monkey` setting: loads FR and KO data
% internally. For each (monkey, area, decoding type), the N_sessions x 3
% stim-pair normalized deltas are pooled into a single sample, and a
% one-sample Wilcoxon signed-rank test is run against H0: median = 0.
% FDR-corrected across all 12 tests.

    monkeys      = {'FR', 'KO'};
    areas        = {'V1', 'V2'};
    dec_col_keep = [1, 2, 3];                       % cols of the 5-col loaded data
    dec_labels   = {'self', 'non-PT', 'PT (rot)'};
    pair_labels  = {'EC-AC', 'EC-EX', 'AC-EX'};
    pair_shapes  = {'o', 's', '^'};

    n_monkeys = numel(monkeys);
    n_areas   = numel(areas);
    n_dec     = numel(dec_labels);
    n_pairs   = numel(pair_labels);

    %% Load and build per-(monkey, area) normalized-delta tensor
    % norm_delta.(monkey).(area) = [n_sessions x n_pairs x n_dec]
    norm_delta = struct();
    n_sessions = struct();

    for m = 1:n_monkeys
        mk = monkeys{m};
        [none_v1, none_v1_meta, none_v2, none_v2_meta] = ...
            load_per_session_data(mk, 'none',   file_location, dec_idx);
        [aff_v1,  aff_v1_meta,  aff_v2,  aff_v2_meta ] = ...
            load_per_session_data(mk, 'affine', file_location, dec_idx);

        assert(none_v1_meta.n_sessions == aff_v1_meta.n_sessions, ...
            '%s V1 session counts differ between modes.', mk);
        assert(none_v2_meta.n_sessions == aff_v2_meta.n_sessions, ...
            '%s V2 session counts differ between modes.', mk);

        for a = 1:n_areas
            area = areas{a};
            if strcmp(area, 'V1')
                none_data = none_v1; aff_data = aff_v1;
                n_s = none_v1_meta.n_sessions;
            else
                none_data = none_v2; aff_data = aff_v2;
                n_s = none_v2_meta.n_sessions;
            end

            tensor = nan(n_s, n_pairs, n_dec);
            for p = 1:n_pairs
                a_mat = aff_data{p}(:,  dec_col_keep);   % [n_s x 3]
                b_mat = none_data{p}(:, dec_col_keep);
                tensor(:, p, :) = 2 .* (a_mat - b_mat) ./ (a_mat + b_mat);
            end

            norm_delta.(mk).(area) = tensor;
            n_sessions.(mk).(area) = n_s;
        end
    end

    %% Pooled Wilcoxon signed-rank tests
    pvals = nan(n_monkeys, n_areas, n_dec);
    for m = 1:n_monkeys
        mk = monkeys{m};
        for a = 1:n_areas
            area = areas{a};
            tensor = norm_delta.(mk).(area);
            for d = 1:n_dec
                vec = tensor(:, :, d);
                vec = vec(:);
                vec = vec(~isnan(vec) & isfinite(vec));
                if numel(vec) >= 1
                    pvals(m, a, d) = signrank(vec);
                end
            end
        end
    end

    %% FDR correction across all 12 tests
    qvals = bh_fdr(pvals(:));
    qvals = reshape(qvals, n_monkeys, n_areas, n_dec);

    %% Print summary
    fprintf('\n=== Pooled-across-stim-pairs Wilcoxon signed-rank (H0: median = 0) ===\n');
    fprintf('FDR-corrected (Benjamini-Hochberg) across %d tests\n\n', numel(pvals));
    fprintf('Monkey | Area | Decoding   |  N  |  raw p   |    q     | sig\n');
    fprintf('%s\n', repmat('-', 1, 62));
    for m = 1:n_monkeys
        mk = monkeys{m};
        for a = 1:n_areas
            area = areas{a};
            n_pooled = n_sessions.(mk).(area) * n_pairs;
            for d = 1:n_dec
                fprintf('%-6s | %-4s | %-9s  | %3d | %.4f   |  %.4f  | %s\n', ...
                    mk, area, dec_labels{d}, n_pooled, ...
                    pvals(m,a,d), qvals(m,a,d), sig_marker(qvals(m,a,d)));
            end
        end
    end
    fprintf('\n');

    %% Y-axis range (symmetric, shared across both subplots)
    all_vals = [];
    for m = 1:n_monkeys
        for a = 1:n_areas
            t = norm_delta.(monkeys{m}).(areas{a});
            all_vals = [all_vals; t(:)]; %#ok<AGROW>
        end
    end
    all_vals = all_vals(~isnan(all_vals) & isfinite(all_vals));
    y_max    = max(0.5, max(abs(all_vals)) * 1.25);

    %% Plot
    v1_base = [0      0.4470 0.7410];
    v2_base = [0.8500 0.3250 0.0980];

    rng(0);   % reproducible jitter

    figure('Position', [100 100 1300 520]);
    ax_handles = gobjects(1, n_monkeys);

    for m = 1:n_monkeys
        mk = monkeys{m};
        ax_handles(m) = subplot(1, n_monkeys, m);
        hold on;

        for a = 1:n_areas
            area = areas{a};
            if strcmp(area, 'V1'), base_color = v1_base; dodge = -0.175;
            else,                  base_color = v2_base; dodge =  0.175;
            end

            n_s    = n_sessions.(mk).(area);
            shades = make_shades(base_color, n_s);
            tensor = norm_delta.(mk).(area);

            for d = 1:n_dec
                x_center = d + dodge;

                % Scatter N_s x N_pairs points (jittered)
                for s = 1:n_s
                    for p = 1:n_pairs
                        x_jit = x_center + (rand - 0.5) * 0.14;
                        y_val = tensor(s, p, d);
                        plot(x_jit, y_val, pair_shapes{p}, ...
                             'MarkerFaceColor', shades(s, :), ...
                             'MarkerEdgeColor', shades(s, :), ...
                             'MarkerSize',      6, ...
                             'HandleVisibility','off');
                    end
                end

                % Cluster center-tendency overlays
                cluster_vals = tensor(:, :, d);
                cluster_vals = cluster_vals(:);
                cluster_vals = cluster_vals(~isnan(cluster_vals) & isfinite(cluster_vals));
                if ~isempty(cluster_vals)
                    med_y  = median(cluster_vals);
                    mean_y = mean(cluster_vals);

                    % Median: thick dark bar
                    plot([x_center - 0.13, x_center + 0.13], [med_y, med_y], '-', ...
                         'Color', [0.15 0.15 0.15], 'LineWidth', 3.0, ...
                         'HandleVisibility','off');
                    % Mean: thin black tick
                    plot([x_center - 0.09, x_center + 0.09], [mean_y, mean_y], '-', ...
                         'Color', [0 0 0], 'LineWidth', 1.2, ...
                         'HandleVisibility','off');
                end

                % Significance asterisk (area-colored)
                sig_y = y_max * 0.92;
                text(x_center, sig_y, sig_marker(qvals(m, a, d)), ...
                     'HorizontalAlignment', 'center', 'FontSize', 12, ...
                     'FontWeight', 'bold', 'Color', base_color);
            end
        end

        yline(0, '--k', 'LineWidth', 1, 'HandleVisibility','off');
        for x = 1:n_dec
            xline(x, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');
        end

        xlim([0.4, n_dec + 0.6]);
        ylim([-y_max, y_max]);
        xticks(1:n_dec);
        xticklabels(dec_labels);
        xtickangle(20);

        title(sprintf('%s   (V1: n=%d sessions, V2: n=%d sessions; %d pairs each)', ...
            mk, n_sessions.(mk).V1, n_sessions.(mk).V2, n_pairs));

        if m == 1
            ylabel('Normalized \Delta: 2(affine - none)/(affine + none)');
        end
        grid on; box on;

        % Legend in the rightmost subplot
        if m == n_monkeys
            hl_v1 = plot(NaN, NaN, 'o', 'MarkerFaceColor', v1_base, ...
                         'MarkerEdgeColor', v1_base, 'MarkerSize', 6);
            hl_v2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', v2_base, ...
                         'MarkerEdgeColor', v2_base, 'MarkerSize', 6);
            hl_p1 = plot(NaN, NaN, pair_shapes{1}, 'MarkerFaceColor', [0.5 0.5 0.5], ...
                         'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 6);
            hl_p2 = plot(NaN, NaN, pair_shapes{2}, 'MarkerFaceColor', [0.5 0.5 0.5], ...
                         'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 6);
            hl_p3 = plot(NaN, NaN, pair_shapes{3}, 'MarkerFaceColor', [0.5 0.5 0.5], ...
                         'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerSize', 6);
            hl_med = plot(NaN, NaN, '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 3.0);
            hl_mn  = plot(NaN, NaN, '-', 'Color', [0 0 0],           'LineWidth', 1.2);
            legend([hl_v1, hl_v2, hl_p1, hl_p2, hl_p3, hl_med, hl_mn], ...
                {'V1', 'V2', pair_labels{1}, pair_labels{2}, pair_labels{3}, ...
                 'median', 'mean'}, ...
                'Location', 'eastoutside', 'NumColumns', 1);
        end
        hold off;
    end

    linkaxes(ax_handles, 'y');
    sgtitle(['Pooled across stim pairs  |  ' ...
             'Wilcoxon signed-rank (FDR-corrected):  * q<0.05,  ** q<0.01,  *** q<0.001']);
end

function s = sig_marker(q)
    if isnan(q),       s = '';
    elseif q < 0.001,  s = '***';
    elseif q < 0.01,   s = '**';
    elseif q < 0.05,   s = '*';
    else,              s = 'n.s.';
    end
end

function q = bh_fdr(p)
% Benjamini-Hochberg FDR correction.
% Input  : vector of raw p-values (NaNs allowed and preserved).
% Output : vector of q-values, same shape as input.
    sz   = size(p);
    pvec = p(:);
    n    = numel(pvec);

    valid = ~isnan(pvec);
    q     = nan(n, 1);

    pv         = pvec(valid);
    m          = numel(pv);
    [ps, sidx] = sort(pv);
    qs         = ps .* m ./ (1:m)';

    % Enforce monotonicity (right to left)
    for k = (m-1):-1:1
        qs(k) = min(qs(k), qs(k+1));
    end
    qs = min(qs, 1);

    q_valid       = zeros(m, 1);
    q_valid(sidx) = qs;
    q(valid)      = q_valid;

    q = reshape(q, sz);
end

function colors = make_shades(base_color, n)
% Generate n shades of base_color by varying saturation only.
% All shades share the same hue (interpolated white -> base_color).
% Session 1 (lightest) -> 40% saturation; session n (darkest) -> base_color.
    if n <= 1
        colors = base_color;
        return;
    end
    white  = [1 1 1];
    t      = linspace(0.4, 1.0, n)';
    colors = (1 - t) .* white + t .* base_color;
end