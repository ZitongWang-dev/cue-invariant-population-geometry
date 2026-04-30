%{
Filename: visualization_procrustes_decoding_per_session.m
Author: Zitong Wang
Date: 2026-04-30

Description:
    Visualizes per-session Procrustes decoding results produced by
    procrustes_decoding_per_session.m. Each session contributes a single
    point per decoding type per stimulus pair (mean across trial samples,
    averaged over the two complementary directions).

    Three figures are produced per monkey:
      Figure 1 - perturb_mode = 'none'   (real noise correlations preserved)
      Figure 2 - perturb_mode = 'affine' (within-session correlations broken)
      Figure 3 - delta (affine - none) per session per decoding type

    Layout (per figure):
      - 3 subplots side-by-side: EC-AC, EC-EX, AC-EX
      - x-axis: 5 decoding types {self, non-PT, PT(rot), PT ctrl, non-PT ctrl}
      - V1 sessions jittered LEFT of each x-tick (blue, graded shades)
      - V2 sessions jittered RIGHT of each x-tick (orange, graded shades)
      - Lines connect points from the same session across decoding types,
        but never across stimulus pairs (subplot boundary enforces this)
      - 2 percent chance line on raw figures, zero line on delta figure
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
% Display order matches the existing bar plot:
%   self -> non-PT -> PT(rot) -> PT ctrl -> non-PT ctrl
dec_idx          = [1, 2, 6, 4, 8];
decoding_labels  = {'self', 'non-PT', 'PT (rot)', 'PT ctrl', 'non-PT ctrl'};

%% Load per-session data for both perturb modes
[none_v1,   none_v1_meta,   none_v2,   none_v2_meta]   = ...
    load_per_session_data(monkey, 'none',   file_location, dec_idx);
[affine_v1, affine_v1_meta, affine_v2, affine_v2_meta] = ...
    load_per_session_data(monkey, 'affine', file_location, dec_idx);

% Sanity: matching session counts across modes
assert(none_v1_meta.n_sessions == affine_v1_meta.n_sessions, ...
    'V1 session counts differ between none and affine.');
assert(none_v2_meta.n_sessions == affine_v2_meta.n_sessions, ...
    'V2 session counts differ between none and affine.');

%% Figure 1: perturb_mode = none
data_none = pack_data(none_v1, none_v1_meta, none_v2, none_v2_meta);
plot_per_session(data_none, monkey, 'none (correlations preserved)', false, decoding_labels);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_none.png']));

%% Figure 2: perturb_mode = affine
data_affine = pack_data(affine_v1, affine_v1_meta, affine_v2, affine_v2_meta);
plot_per_session(data_affine, monkey, 'affine (within-session correlations broken)', false, decoding_labels);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_affine.png']));

%% Figure 3: delta = affine - none
delta_v1 = cellfun(@(a, b) a - b, affine_v1, none_v1, 'UniformOutput', false);
delta_v2 = cellfun(@(a, b) a - b, affine_v2, none_v2, 'UniformOutput', false);
data_delta = pack_data(delta_v1, none_v1_meta, delta_v2, none_v2_meta);
plot_per_session(data_delta, monkey, '\Delta (affine - none)', true, decoding_labels);
% saveas(gcf, fullfile(fig_save_dir, ...
%     ['procrustes_decoding_per_session_', monkey, '_delta.png']));

%% =================== LOCAL FUNCTIONS ===================
function data = pack_data(v1_data, v1_meta, v2_data, v2_meta)
    data.v1_data = v1_data;
    data.v1_meta = v1_meta;
    data.v2_data = v2_data;
    data.v2_meta = v2_meta;
end

function [v1_data, v1_meta, v2_data, v2_meta] = ...
    load_per_session_data(monkey, perturb_mode, file_location, dec_idx)
% Loads all 6 per-pair files for V1 and V2, combines complementary
% directions per session, and returns a {3} cell array per area where
% each cell is [n_sessions x numel(dec_idx)].
%
% Output cell order: {EC-AC, EC-EX, AC-EX}
%   EC-AC <- combine acec + ecac
%   EC-EX <- combine ecex + exec
%   AC-EX <- combine acex + exac

    pair_names = {'acec', 'ecac', 'ecex', 'exec', 'acex', 'exac'};

    % --- V1 ---
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

    % --- V2 ---
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
% Per session: stack the two complementary directions vertically
% [trial_sample_repeat x 8] each -> [2*trial_sample_repeat x 8],
% then take the column mean and select the dec_idx columns.
    n        = numel(pair1_struct);
    combined = zeros(n, numel(dec_idx));
    for s = 1:n
        stacked      = [pair1_struct(s).accuracy; pair2_struct(s).accuracy];
        col_mean     = mean(stacked, 1);
        combined(s,:) = col_mean(dec_idx);
    end
end

function plot_per_session(data, monkey, title_str, is_delta, decoding_labels)
    pair_labels = {'EC-AC', 'EC-EX', 'AC-EX'};
    n_dec       = numel(decoding_labels);
    x_base      = 1:n_dec;

    % Base colors (matching the existing visualization)
    v1_base = [0      0.4470 0.7410];
    v2_base = [0.8500 0.3250 0.0980];

    n_v1 = data.v1_meta.n_sessions;
    n_v2 = data.v2_meta.n_sessions;

    v1_colors = make_shades(v1_base, n_v1);
    v2_colors = make_shades(v2_base, n_v2);

    % Per-session x-jitter (fixed, not random) so lines look near-vertical
    if n_v1 > 1
        v1_jitter = linspace(-0.30, -0.05, n_v1);
    else
        v1_jitter = -0.175;
    end
    if n_v2 > 1
        v2_jitter = linspace( 0.05,  0.30, n_v2);
    else
        v2_jitter = 0.175;
    end

    figure('Position', [100 100 1300 480]);
    ax_handles = gobjects(1, 3);

    for p = 1:3
        ax_handles(p) = subplot(1, 3, p);
        hold on;

        v1_pair = data.v1_data{p};   % [n_v1 x n_dec]
        v2_pair = data.v2_data{p};   % [n_v2 x n_dec]

        % V1 sessions: line + markers per session
        for s = 1:size(v1_pair, 1)
            plot(x_base + v1_jitter(s), v1_pair(s, :), '-o', ...
                 'Color',           v1_colors(s, :), ...
                 'MarkerFaceColor', v1_colors(s, :), ...
                 'MarkerEdgeColor', v1_colors(s, :), ...
                 'MarkerSize',      5, ...
                 'LineWidth',       1.0, ...
                 'HandleVisibility','off');
        end
        % V2 sessions
        for s = 1:size(v2_pair, 1)
            plot(x_base + v2_jitter(s), v2_pair(s, :), '-o', ...
                 'Color',           v2_colors(s, :), ...
                 'MarkerFaceColor', v2_colors(s, :), ...
                 'MarkerEdgeColor', v2_colors(s, :), ...
                 'MarkerSize',      5, ...
                 'LineWidth',       1.0, ...
                 'HandleVisibility','off');
        end

        % Reference line
        if is_delta
            yline(0,    '--k', 'LineWidth', 1, 'HandleVisibility','off');
        else
            yline(0.02, '--k', 'LineWidth', 1, 'HandleVisibility','off');
        end

        % Vertical separators between V1 (left) and V2 (right) at each x-tick
        % (kept light so they don't dominate)
        for x = x_base
            xline(x, ':', 'Color', [0.85 0.85 0.85], 'HandleVisibility','off');
        end

        xlim([0.4, n_dec + 0.6]);
        xticks(x_base);
        xticklabels(decoding_labels);
        xtickangle(20);
        title(pair_labels{p});

        if p == 1
            if is_delta
                ylabel('\Delta decoding accuracy');
            else
                ylabel('Decoding accuracy');
            end
        end

        grid on; box on;

        % Legend lives in the rightmost subplot, using the V1/V2 base colors
        if p == 3
            hl_v1 = plot(NaN, NaN, '-o', ...
                'Color', v1_base, 'MarkerFaceColor', v1_base, ...
                'MarkerSize', 5, 'LineWidth', 1.0);
            hl_v2 = plot(NaN, NaN, '-o', ...
                'Color', v2_base, 'MarkerFaceColor', v2_base, ...
                'MarkerSize', 5, 'LineWidth', 1.0);
            hl_chance = plot(NaN, NaN, '--k', 'LineWidth', 1);
            if is_delta
                legend([hl_v1, hl_v2, hl_chance], ...
                    {sprintf('V1 (n=%d sessions)', n_v1), ...
                     sprintf('V2 (n=%d sessions)', n_v2), ...
                     'no effect'}, ...
                    'Location', 'best');
            else
                legend([hl_v1, hl_v2, hl_chance], ...
                    {sprintf('V1 (n=%d sessions)', n_v1), ...
                     sprintf('V2 (n=%d sessions)', n_v2), ...
                     'chance (2%)'}, ...
                    'Location', 'best');
            end
        end

        hold off;
    end

    % Share y-axis across the three pair-panels
    linkaxes(ax_handles, 'y');

    sgtitle(sprintf('%s | per-session Procrustes decoding | %s', monkey, title_str));
end

function colors = make_shades(base_color, n)
% Generate n shades of base_color by varying saturation only.
% All shades have the same hue as base_color (white -> base_color line in RGB).
% Session 1 (lightest) -> 40% saturation; session n (darkest) -> base_color exactly.
    if n <= 1
        colors = base_color;
        return;
    end
    white  = [1 1 1];
    t      = linspace(0.4, 1.0, n)';        % saturation factor: 0=white, 1=base
    colors = (1 - t) .* white + t .* base_color;
end