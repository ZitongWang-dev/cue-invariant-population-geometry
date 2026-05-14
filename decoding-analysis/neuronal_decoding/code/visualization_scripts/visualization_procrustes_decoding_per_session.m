%{
Filename: visualization_procrustes_decoding_per_session.m
Author: Zitong Wang
Date: 2026-04-30

Description:
    Visualizes per-session Procrustes decoding results produced by
    procrustes_decoding_per_session.m.

    Five figures are produced per monkey:
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
%}

clc; clear;

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
%
% Drops the two control columns (PT ctrl, non-PT ctrl) which sit at chance
% and would produce unstable ratios.

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

        % V1 sessions
        for s = 1:n_v1
            plot(x_base + v1_jitter(s), v1_norm{p}(s, :), '-o', ...
                 'Color', v1_colors(s, :), 'MarkerFaceColor', v1_colors(s, :), ...
                 'MarkerEdgeColor', v1_colors(s, :), 'MarkerSize', 5, ...
                 'LineWidth', 1.0, 'HandleVisibility','off');
        end
        % V2 sessions
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

        % Significance markers above each (decoding type, area) cluster
        sig_y = y_max * 0.92;
        for d = 1:n_dec
            text(d - 0.175, sig_y, sig_marker(q_v1(p, d)), ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', v1_base);
            text(d + 0.175, sig_y, sig_marker(q_v2(p, d)), ...
                'HorizontalAlignment', 'center', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', v2_base);
        end

        % Legend in last subplot
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