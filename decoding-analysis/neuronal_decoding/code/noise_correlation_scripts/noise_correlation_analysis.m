%{
Filename: noise_correlation_analysis.m
Author: Zitong Wang
Date: 2025-03-25

Description:
    Characterizes noise correlation structure in pseudopopulation data from
    V1 and V2. For each neuron pair and each stimulus, computes pairwise
    noise correlation and geometric mean firing rate. All three renderings
    (AC, EC, EX) are processed in a single run.

    Stage 1: Examines the relationship between geometric mean firing rate
    and noise correlation across stimuli for each neuron pair. Produces
    example scatter plots (40 pairs per rendering) and a combined
    distribution of per-pair geo-mean-NC correlations across renderings.

    Stage 2: For each neuron pair, computes the mean and variance of noise
    correlation across all stimuli. Plots mean vs. variance with marginal
    distributions.

    Shuffle analysis: Optionally shuffles trial order independently per
    neuron per stimulus to destroy noise correlations while preserving
    marginal firing rate distributions. Produces comparison figures and
    significance tests.

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat

Outputs:
    - MAT-files: results/noise_correlation_outputs/<monkey>/<area>/
    - Figures:   results/figures/noise_correlation_<monkey>_<area>/
%}

%% Initialization
clc; clear;
close all

%% Configuration
monkey     = 'FR';       % 'FR' or 'KO'
vp         = 'V1';       % 'V1' or 'V2'
timewindow = [330 630];  % spike-count window (ms)

% Shuffle configuration
run_shuffle = true;
n_shuffles  = 10;

renderings = {'ac', 'ec', 'ex'};
rendering_colors = struct('ac', [0.85 0.33 0.10], ...   % orange
                          'ec', [0.00 0.45 0.74], ...   % blue
                          'ex', [0.47 0.67 0.19]);      % green
N_renderings = numel(renderings);

%% Load data
data_file = fullfile('..','..','neuronal_data', monkey, vp, ...
    sprintf('%s_%s_allstim.mat', monkey, vp));
tmp = load(data_file, 'three_stim_array');
spike_data = tmp.three_stim_array;
clear tmp;

%% Prepare output directories
save_path_data = fullfile('..','..','results','noise_correlation_outputs', ...
    monkey, vp);
save_path_figs = fullfile('..','..','results','figures', ...
    sprintf('noise_correlation_%s_%s', monkey, vp));
if ~exist(save_path_data, 'dir'), mkdir(save_path_data); end
if ~exist(save_path_figs, 'dir'), mkdir(save_path_figs); end

fig_prefix = sprintf('nc_%s_%s', monkey, vp);

%% Constants
name2idx = struct('ac',1,'ec',2,'ex',3);
N_stimuli = 50;
N_trials_per_stim = 10;
rate_floor = 0.1;

%% Get neuron count and pair indices
[trial_data_tmp, ~] = multiclass_svmloader_PT(spike_data{1}, timewindow);
N_neurons = size(trial_data_tmp, 2);
N_pairs = nchoosek(N_neurons, 2);
pair_idx = nchoosek(1:N_neurons, 2);
clear trial_data_tmp;

fprintf('Monkey %s | %s | %d neurons | %d pairs\n', ...
    monkey, vp, N_neurons, N_pairs);

%% Select 40 example pairs (shared across renderings)
n_examples = 40;
rng(42);
example_pair_indices = randsample(N_pairs, n_examples);

%% Preallocate storage
all_results = struct();

%% ===== MAIN LOOP: Process each rendering =====
for r = 1:N_renderings
    rend = renderings{r};
    rendering_idx = name2idx.(rend);
    fprintf('\n========== Processing %s ==========\n', upper(rend));

    % Load and format trial data
    [trial_data, labels] = multiclass_svmloader_PT(spike_data{rendering_idx}, timewindow);

    % Reorganize: [N_stimuli x N_trials_per_stim x N_neurons]
    data_by_stim = zeros(N_stimuli, N_trials_per_stim, N_neurons);
    for s = 1:N_stimuli
        data_by_stim(s, :, :) = trial_data(labels == s, :);
    end

    % Mean firing rate per stimulus per neuron: [N_stimuli x N_neurons]
    mean_rates = squeeze(mean(data_by_stim, 2));
    mean_rates_floored = max(mean_rates, rate_floor);

    % ---- Original data: compute NC and geo mean ----
    [nc_per_stim, geomean_per_stim] = compute_nc_geomean( ...
        data_by_stim, mean_rates_floored, pair_idx, N_pairs, N_stimuli);

    % Per-pair correlation between geo mean and NC
    geomean_nc_corr = compute_geomean_nc_corr(geomean_per_stim, nc_per_stim, N_pairs);

    % Mean and variance of NC across stimuli
    mean_nc = mean(nc_per_stim, 2);
    var_nc  = var(nc_per_stim, 0, 2);

    % Store original results
    all_results.(rend).nc_per_stim      = nc_per_stim;
    all_results.(rend).geomean_per_stim = geomean_per_stim;
    all_results.(rend).geomean_nc_corr  = geomean_nc_corr;
    all_results.(rend).mean_nc          = mean_nc;
    all_results.(rend).var_nc           = var_nc;

    % ---- Shuffle analysis ----
    if run_shuffle
        fprintf('Running %d shuffles for %s...\n', n_shuffles, upper(rend));

        shuf_geomean_nc_corr = zeros(N_pairs, n_shuffles);
        shuf_mean_nc         = zeros(N_pairs, n_shuffles);
        shuf_var_nc          = zeros(N_pairs, n_shuffles);

        % Summary statistics per shuffle
        shuf_pop_mean_nc            = zeros(n_shuffles, 1);
        shuf_pop_median_geomean_nc  = zeros(n_shuffles, 1);

        for sh = 1:n_shuffles
            % Shuffle: independently permute 10 trials per neuron per stimulus
            data_by_stim_shuf = zeros(size(data_by_stim));
            for n = 1:N_neurons
                for s = 1:N_stimuli
                    data_by_stim_shuf(s, :, n) = data_by_stim(s, randperm(N_trials_per_stim), n);
                end
            end

            % Compute NC (geo means unchanged, reuse mean_rates_floored)
            [nc_shuf, ~] = compute_nc_geomean( ...
                data_by_stim_shuf, mean_rates_floored, pair_idx, N_pairs, N_stimuli);

            % Per-pair geo-mean-NC correlation
            shuf_geomean_nc_corr(:, sh) = compute_geomean_nc_corr( ...
                geomean_per_stim, nc_shuf, N_pairs);

            % Mean and variance of NC across stimuli
            shuf_mean_nc(:, sh) = mean(nc_shuf, 2);
            shuf_var_nc(:, sh)  = var(nc_shuf, 0, 2);

            % Population-level summary statistics
            shuf_pop_mean_nc(sh)           = mean(mean(nc_shuf, 2));
            shuf_pop_median_geomean_nc(sh) = median(shuf_geomean_nc_corr(:, sh));

            if mod(sh, 50) == 0
                fprintf('  Shuffle %d / %d done\n', sh, n_shuffles);
            end
        end

        % Average across shuffles (per-pair)
        all_results.(rend).shuf_geomean_nc_corr_avg = mean(shuf_geomean_nc_corr, 2);
        all_results.(rend).shuf_mean_nc_avg         = mean(shuf_mean_nc, 2);
        all_results.(rend).shuf_var_nc_avg          = mean(shuf_var_nc, 2);

        % Full per-shuffle data (for pooled visualization)
        all_results.(rend).shuf_mean_nc_all = shuf_mean_nc; % [N_pairs x n_shuffles]
        all_results.(rend).shuf_var_nc_all  = shuf_var_nc;  % [N_pairs x n_shuffles]

        % Summary stat distributions
        all_results.(rend).shuf_pop_mean_nc           = shuf_pop_mean_nc;
        all_results.(rend).shuf_pop_median_geomean_nc = shuf_pop_median_geomean_nc;

        % Real summary stats for comparison
        all_results.(rend).real_pop_mean_nc           = mean(mean_nc);
        all_results.(rend).real_pop_median_geomean_nc = median(geomean_nc_corr);

        fprintf('Shuffle analysis complete for %s.\n', upper(rend));
    end
end

%% ==================== ORIGINAL FIGURES ====================

%% Determine shared x-axis range for example scatter plots
all_geomeans = [];
for r = 1:N_renderings
    rend = renderings{r};
    all_geomeans = [all_geomeans; ...
        all_results.(rend).geomean_per_stim(example_pair_indices, :)];
end
shared_x_max = ceil(max(all_geomeans(:)));
shared_xlim_examples = [0 shared_x_max];
shared_ylim_examples = [-1 1];

%% FIGURE 1-3: Stage 1 Example Scatter Plots (one per rendering)
n_cols = 8;
n_rows = ceil(n_examples / n_cols);

for r = 1:N_renderings
    rend = renderings{r};
    clr = rendering_colors.(rend);

    fig = figure('Position', [50 50 2000 n_rows * 200]);
    for k = 1:n_examples
        p = example_pair_indices(k);
        ni = pair_idx(p, 1);
        nj = pair_idx(p, 2);

        subplot(n_rows, n_cols, k);
        scatter(all_results.(rend).geomean_per_stim(p, :), ...
                all_results.(rend).nc_per_stim(p, :), ...
                15, clr, 'filled', 'MarkerFaceAlpha', 0.6);
        hold on;
        coeffs = polyfit(all_results.(rend).geomean_per_stim(p, :), ...
                         all_results.(rend).nc_per_stim(p, :), 1);
        x_line = linspace(shared_xlim_examples(1), shared_xlim_examples(2), 100);
        plot(x_line, polyval(coeffs, x_line), 'r-', 'LineWidth', 1);
        yline(0, 'k:', 'LineWidth', 0.5);
        hold off;

        xlim(shared_xlim_examples);
        ylim(shared_ylim_examples);

        title(sprintf('(%d,%d) r=%.2f', ni, nj, ...
            all_results.(rend).geomean_nc_corr(p)), 'FontSize', 7);
        if k > (n_rows - 1) * n_cols
            xlabel('Geo mean', 'FontSize', 7);
        end
        if mod(k - 1, n_cols) == 0
            ylabel('NC', 'FontSize', 7);
        end
        set(gca, 'FontSize', 6);
    end
    sgtitle(sprintf('%s %s %s — Geo Mean vs Noise Correlation (40 example pairs)', ...
        monkey, vp, upper(rend)), 'FontSize', 12);
    % saveas(fig, fullfile(save_path_figs, ...
    %     sprintf('%s_%s_stage1_example_scatter.png', fig_prefix, rend)));
end

%% FIGURE 4: Stage 1 Distribution of Geo Mean-NC Correlations (3 renderings)
fig4 = figure('Position', [100 100 900 600]);
hold on;
legend_entries = cell(1, N_renderings);
for r = 1:N_renderings
    rend = renderings{r};
    clr = rendering_colors.(rend);
    histogram(all_results.(rend).geomean_nc_corr, 50, ...
        'FaceColor', clr, 'EdgeColor', 'w', 'FaceAlpha', 0.5);
    legend_entries{r} = sprintf('%s (med=%.3f)', upper(rend), ...
        median(all_results.(rend).geomean_nc_corr));
end
xline(0, 'k--', 'LineWidth', 1);
hold off;
xlabel('Correlation (geo mean vs NC)');
ylabel('Number of neuron pairs');
title(sprintf('%s %s — Distribution of Geo Mean-NC Correlations', monkey, vp));
legend(legend_entries, 'Location', 'best');
% saveas(fig4, fullfile(save_path_figs, ...
%     sprintf('%s_stage1_geomean_nc_corr_distribution.png', fig_prefix)));

%% FIGURE 5a-c: Stage 2 Scatterhist (one per rendering)
all_mean_nc_cat = [];
all_var_nc_cat  = [];
for r = 1:N_renderings
    rend = renderings{r};
    all_mean_nc_cat = [all_mean_nc_cat; all_results.(rend).mean_nc];
    all_var_nc_cat  = [all_var_nc_cat;  all_results.(rend).var_nc];
end
shared_xlim_s2 = [floor(min(all_mean_nc_cat)*20)/20, ceil(max(all_mean_nc_cat)*20)/20];
shared_ylim_s2 = [floor(min(all_var_nc_cat)*100)/100, ceil(max(all_var_nc_cat)*100)/100];

for r = 1:N_renderings
    rend = renderings{r};
    clr = rendering_colors.(rend);
    mn = all_results.(rend).mean_nc;
    vn = all_results.(rend).var_nc;

    fig5 = figure('Position', [50 + (r-1)*500, 100, 800, 800]);
    h = scatterhist(mn, vn, ...
        'Direction', 'out', 'Color', clr, 'Marker', '.', 'MarkerSize', 8);

    xlim(h(1), shared_xlim_s2);
    ylim(h(1), shared_ylim_s2);

    % Regression line
    hold(h(1), 'on');
    valid = ~isnan(mn) & ~isnan(vn);
    coeffs = polyfit(mn(valid), vn(valid), 1);
    x_line = linspace(shared_xlim_s2(1), shared_xlim_s2(2), 100);
    plot(h(1), x_line, polyval(coeffs, x_line), 'r-', 'LineWidth', 1.5);
    R2 = corrcoef(mn(valid), vn(valid));
    text(h(1), 0.05, 0.95, sprintf('r = %.3f', R2(1,2)), ...
        'Units', 'normalized', 'FontSize', 11, 'VerticalAlignment', 'top');
    hold(h(1), 'off');

    % Marginal annotations: mean (red) and median (black)
    hold(h(2), 'on');
    y_top = ylim(h(2));
    plot(h(2), [mean(mn) mean(mn)], y_top, 'r-', 'LineWidth', 1.5);
    plot(h(2), [median(mn) median(mn)], y_top, 'k--', 'LineWidth', 1.5);
    text(h(2), 0.02, 0.95, sprintf('mean: %.3f', mean(mn)), ...
        'Units', 'normalized', 'FontSize', 10, 'Color', 'r', 'VerticalAlignment', 'top');
    text(h(2), 0.02, 0.75, sprintf('median: %.3f', median(mn)), ...
        'Units', 'normalized', 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'top');
    hold(h(2), 'off');

    hold(h(3), 'on');
    y_right = ylim(h(3));
    plot(h(3), [mean(vn) mean(vn)], y_right, 'r-', 'LineWidth', 1.5);
    plot(h(3), [median(vn) median(vn)], y_right, 'k--', 'LineWidth', 1.5);
    text(h(3), 0.02, 0.95, sprintf('mean: %.4f', mean(vn)), ...
        'Units', 'normalized', 'FontSize', 10, 'Color', 'r', 'VerticalAlignment', 'top');
    text(h(3), 0.02, 0.75, sprintf('median: %.4f', median(vn)), ...
        'Units', 'normalized', 'FontSize', 10, 'Color', 'k', 'VerticalAlignment', 'top');
    hold(h(3), 'off');

    xlabel(h(1), 'Mean noise correlation');
    ylabel(h(1), 'Variance of noise correlation');
    title(h(1), sprintf('%s %s %s — Mean vs Variance of NC (%d pairs)', ...
        monkey, vp, upper(rend), N_pairs));

    % saveas(fig5, fullfile(save_path_figs, ...
    %     sprintf('%s_%s_stage2_mean_vs_variance_nc.png', fig_prefix, rend)));
end

%% ==================== SHUFFLE FIGURES ====================
if run_shuffle

%% FIGURE 6: Stage 1 Geo Mean-NC Correlation: Original vs Shuffle (3 subplots)
fig6 = figure('Position', [100 100 2500 500]);
for r = 1:N_renderings
    rend = renderings{r};
    clr = rendering_colors.(rend);

    orig_vals = all_results.(rend).geomean_nc_corr;
    shuf_vals = all_results.(rend).shuf_geomean_nc_corr_avg;

    subplot(1, 3, r);
    hold on;
    % Shuffle-averaged distribution (black)
    histogram(shuf_vals, 50, ...
        'FaceColor', [0 0 0], 'EdgeColor', 'w', 'FaceAlpha', 1);
    % Original distribution (colored)
    histogram(orig_vals, 50, ...
        'FaceColor', clr, 'EdgeColor', 'w', 'FaceAlpha', 0.5);

    xline(0, 'k:', 'LineWidth', 0.5);
    hold off;

    xlabel('Correlation (geo mean vs NC)');
    ylabel('Number of neuron pairs');
    title(sprintf('%s', upper(rend)));

    % Text annotations
    text(0.02, 0.95, sprintf('Original mean: %.3f', mean(orig_vals)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', clr, ...
        'VerticalAlignment', 'top', 'FontWeight', 'bold');
    text(0.02, 0.87, sprintf('Original median: %.3f', median(orig_vals)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', clr, ...
        'VerticalAlignment', 'top');
    text(0.02, 0.79, sprintf('Shuffle mean: %.3f', mean(shuf_vals)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', [0.3 0.3 0.3], ...
        'VerticalAlignment', 'top', 'FontWeight', 'bold');
    text(0.02, 0.71, sprintf('Shuffle median: %.3f', median(shuf_vals)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', [0.3 0.3 0.3], ...
        'VerticalAlignment', 'top');

    legend('Shuffle (avg per pair)', 'Original', 'Location', 'northeast');
end
sgtitle(sprintf('%s %s — Geo Mean-NC Correlation: Original vs Shuffle (%d shuffles)', ...
    monkey, vp, n_shuffles), 'FontSize', 12);
% saveas(fig6, fullfile(save_path_figs, ...
%     sprintf('%s_stage1_geomean_nc_corr_original_vs_shuffle.png', fig_prefix)));

%% FIGURE 7a-c: Stage 2 Scatterhist Original vs Pooled Shuffle (one per rendering)
for r = 1:N_renderings
    rend = renderings{r};
    clr = rendering_colors.(rend);

    % Original data
    mn_orig = all_results.(rend).mean_nc;
    vn_orig = all_results.(rend).var_nc;

    % Pool all shuffle data
    mn_shuf_pooled = reshape(all_results.(rend).shuf_mean_nc_all, [], 1);
    vn_shuf_pooled = reshape(all_results.(rend).shuf_var_nc_all, [], 1);

    % Concatenate: shuffle first so original plots on top
    all_mn = [mn_shuf_pooled; mn_orig];
    all_vn = [vn_shuf_pooled; vn_orig];
    group_labels = [repmat({'Shuffle'}, length(mn_shuf_pooled), 1); ...
                    repmat({'Original'}, length(mn_orig), 1)];

    fig7 = figure('Position', [50 + (r-1)*500, 100, 800, 800]);
    h = scatterhist(all_mn, all_vn, 'Group', group_labels, ...
        'Color', [0.7 0.7 0.7; clr], ...
        'Marker', '..', ...
        'MarkerSize', [4 8], ...
        'Kernel', 'on', ...
        'Direction', 'out');

    xlim(h(1), shared_xlim_s2);
    ylim(h(1), shared_ylim_s2);

    % Annotate top marginal (mean NC distribution)
    text(h(2), 0.8, 0.95, sprintf('Orig mean: %.3f', mean(mn_orig)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', clr, ...
        'VerticalAlignment', 'top', 'FontWeight', 'bold');
    text(h(2), 0.8, 0.80, sprintf('Orig median: %.3f', median(mn_orig)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', clr, ...
        'VerticalAlignment', 'top');
    text(h(2), 0.8, 0.65, sprintf('Shuf mean: %.3f', mean(mn_shuf_pooled)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
        'VerticalAlignment', 'top', 'FontWeight', 'bold');
    text(h(2), 0.8, 0.50, sprintf('Shuf median: %.3f', median(mn_shuf_pooled)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
        'VerticalAlignment', 'top');

    % Annotate right marginal (variance of NC distribution)
    text(h(3), 0.02, 0.95, sprintf('Orig mean: %.4f', mean(vn_orig)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', clr, ...
        'VerticalAlignment', 'top', 'FontWeight', 'bold');
    text(h(3), 0.02, 0.9, sprintf('Orig median: %.4f', median(vn_orig)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', clr, ...
        'VerticalAlignment', 'top');
    text(h(3), 0.02, 0.85, sprintf('Shuf mean: %.4f', mean(vn_shuf_pooled)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
        'VerticalAlignment', 'top', 'FontWeight', 'bold');
    text(h(3), 0.02, 0.80, sprintf('Shuf median: %.4f', median(vn_shuf_pooled)), ...
        'Units', 'normalized', 'FontSize', 9, 'Color', [0.4 0.4 0.4], ...
        'VerticalAlignment', 'top');

    xlabel(h(1), 'Mean noise correlation');
    ylabel(h(1), 'Variance of noise correlation');
    title(h(1), sprintf('%s %s %s — Original vs Shuffle (%d pairs, %d shuffles pooled)', ...
        monkey, vp, upper(rend), N_pairs, n_shuffles));
    legend(h(1), 'Location', 'best');

    % saveas(fig7, fullfile(save_path_figs, ...
    %     sprintf('%s_%s_stage2_mean_vs_variance_nc_vs_shuffle.png', fig_prefix, rend)));
end

%% FIGURE 8: Significance Test — Real vs Shuffle Distribution (summary stats)
fig8 = figure('Position', [100 100 2000 800]);

for r = 1:N_renderings
    rend = renderings{r};
    clr = rendering_colors.(rend);

    % --- Top row: Population mean of noise correlation ---
    subplot(2, 3, r);
    hold on;
    histogram(all_results.(rend).shuf_pop_mean_nc, 30, ...
        'FaceColor', [0 0 0], 'EdgeColor', 'w');
    xline(all_results.(rend).real_pop_mean_nc, '-', 'Color', clr, 'LineWidth', 2);
    hold off;

    real_val = all_results.(rend).real_pop_mean_nc;
    shuf_vals = all_results.(rend).shuf_pop_mean_nc;
    p_val = mean(shuf_vals >= real_val);
    xlabel('Population mean NC');
    ylabel('Shuffle count');
    title(sprintf('%s — Pop. mean NC\nreal=%.4f, p=%.3f', upper(rend), real_val, p_val));

    % --- Bottom row: Median of geo-mean-NC correlation ---
    subplot(2, 3, r + 3);
    hold on;
    histogram(all_results.(rend).shuf_pop_median_geomean_nc, 30, ...
        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'w');
    xline(all_results.(rend).real_pop_median_geomean_nc, '-', 'Color', clr, 'LineWidth', 2);
    hold off;

    real_val = all_results.(rend).real_pop_median_geomean_nc;
    shuf_vals = all_results.(rend).shuf_pop_median_geomean_nc;
    p_val = mean(abs(shuf_vals) >= abs(real_val));
    xlabel('Median geo-mean-NC correlation');
    ylabel('Shuffle count');
    title(sprintf('%s — Med. geo-NC corr\nreal=%.4f, p=%.3f', upper(rend), real_val, p_val));
end

sgtitle(sprintf('%s %s — Real vs Shuffle Distribution (%d shuffles)', ...
    monkey, vp, n_shuffles), 'FontSize', 12);
% saveas(fig8, fullfile(save_path_figs, ...
%     sprintf('%s_significance_real_vs_shuffle.png', fig_prefix)));

end % end shuffle figures

%% Save all results
results = struct();
results.monkey     = monkey;
results.vp         = vp;
results.timewindow = timewindow;
results.pair_idx   = pair_idx;
results.example_pair_indices = example_pair_indices;
results.n_shuffles = n_shuffles;
results.renderings = all_results;

% save(fullfile(save_path_data, 'noise_correlation_results.mat'), 'results', '-v7.3');
% fprintf('\nAll results saved to %s\n', save_path_data);

%% ==================== HELPER FUNCTIONS ====================

function [nc_per_stim, geomean_per_stim] = compute_nc_geomean( ...
    data_by_stim, mean_rates_floored, pair_idx, N_pairs, N_stimuli)
% Compute pairwise noise correlation and geometric mean for all pairs and stimuli

    nc_per_stim      = zeros(N_pairs, N_stimuli);
    geomean_per_stim = zeros(N_pairs, N_stimuli);

    for p = 1:N_pairs
        ni = pair_idx(p, 1);
        nj = pair_idx(p, 2);
        for s = 1:N_stimuli
            ri = data_by_stim(s, :, ni)';
            rj = data_by_stim(s, :, nj)';
            R = corrcoef(ri, rj);
            nc_per_stim(p, s) = R(1, 2);
            geomean_per_stim(p, s) = sqrt(mean_rates_floored(s, ni) * mean_rates_floored(s, nj));
        end
    end

    % Handle NaN
    nan_count = sum(isnan(nc_per_stim(:)));
    if nan_count > 0
        fprintf('  Warning: %d NaN noise correlations. Setting to 0.\n', nan_count);
        nc_per_stim(isnan(nc_per_stim)) = 0;
    end
end

function geomean_nc_corr = compute_geomean_nc_corr(geomean_per_stim, nc_per_stim, N_pairs)
% Compute per-pair correlation between geometric mean and noise correlation
    geomean_nc_corr = zeros(N_pairs, 1);
    for p = 1:N_pairs
        R = corrcoef(geomean_per_stim(p, :), nc_per_stim(p, :));
        geomean_nc_corr(p) = R(1, 2);
    end
    geomean_nc_corr(isnan(geomean_nc_corr)) = 0;
end