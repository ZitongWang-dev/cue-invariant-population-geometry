%{
Filename: noise_correlation_analysis.m
Author: Zitong Wang
Date: 2025-03-25

Description:
    Characterizes noise correlation structure in pseudopopulation data from
    V1 and V2. For each neuron pair and each stimulus, computes pairwise
    noise correlation and geometric mean firing rate.

    Stage 1: Examines the relationship between geometric mean firing rate
    and noise correlation across stimuli for each neuron pair. Produces
    example scatter plots and a distribution of per-pair geo-mean–NC
    correlations.

    Stage 2: For each neuron pair, computes the mean and variance of noise
    correlation across all stimuli. Plots mean vs. variance for the full
    population of pairs.

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat
        Contains variable 'three_stim_array': cell array of [trial x neuron]
        matrices for each stimulus rendering.

Outputs:
    - Figures saved to:
      results/noise_correlation_results/<monkey>/<area>/<rendering>/
%}

%% Initialization
clc; clear;

%% Configuration
monkey    = 'KO';       % 'FR' or 'KO'
vp        = 'V1';       % 'V1' or 'V2'
rendering = 'ec';       % 'ac', 'ec', or 'ex'
timewindow = [330 630]; % spike-count window (ms)

%% Load data
data_file = fullfile('..','..','neuronal_data', monkey, vp, ...
    sprintf('%s_%s_allstim.mat', monkey, vp));
tmp = load(data_file, 'three_stim_array');
spike_data = tmp.three_stim_array;
clear tmp;

%% Prepare output directories
save_path_data = fullfile('..','..','results','noise_correlation_outputs', ...
    monkey, vp, rendering);
save_path_figs = fullfile('..','..','results','figures');
if ~exist(save_path_data, 'dir'), mkdir(save_path_data); end
if ~exist(save_path_figs, 'dir'), mkdir(save_path_figs); end

% Figure filename prefix for this configuration
fig_prefix = sprintf('nc_%s_%s_%s', monkey, vp, rendering);

%% Select rendering and load trial data
name2idx = struct('ac',1,'ec',2,'ex',3);
rendering_idx = name2idx.(rendering);

[trial_data, labels] = multiclass_svmloader_PT(spike_data{rendering_idx}, timewindow);
% trial_data: [500 x N_neurons], labels: [500 x 1] (1..50, 10 trials each)

[N_trials, N_neurons] = size(trial_data);
N_stimuli = 50;
N_trials_per_stim = 10;
N_pairs = nchoosek(N_neurons, 2);

fprintf('Monkey %s | %s | %s | %d neurons | %d pairs\n', ...
    monkey, vp, upper(rendering), N_neurons, N_pairs);

%% Reorganize data: [N_stimuli x N_trials_per_stim x N_neurons]
data_by_stim = zeros(N_stimuli, N_trials_per_stim, N_neurons);
for s = 1:N_stimuli
    trial_idx = (labels == s);
    data_by_stim(s, :, :) = trial_data(trial_idx, :);
end

%% Compute mean firing rate per stimulus per neuron: [N_stimuli x N_neurons]
mean_rates = squeeze(mean(data_by_stim, 2)); % [50 x N_neurons]

%% Floor small/zero firing rates for geometric mean computation
rate_floor = 0.1; % small floor to avoid zero geometric means
mean_rates_floored = max(mean_rates, rate_floor);

%% Generate all neuron pair indices
pair_idx = nchoosek(1:N_neurons, 2); % [N_pairs x 2]

%% Stage 1 & 2: Compute per-pair, per-stimulus noise correlation and geo mean
% Preallocate
nc_per_stim    = zeros(N_pairs, N_stimuli); % noise corr per pair per stimulus
geomean_per_stim = zeros(N_pairs, N_stimuli); % geo mean per pair per stimulus

fprintf('Computing pairwise noise correlations...\n');
for p = 1:N_pairs
    ni = pair_idx(p, 1);
    nj = pair_idx(p, 2);

    for s = 1:N_stimuli
        % Extract 10 trials for both neurons for this stimulus
        ri = data_by_stim(s, :, ni)'; % [10 x 1]
        rj = data_by_stim(s, :, nj)'; % [10 x 1]

        % Noise correlation (Pearson r of trial-to-trial fluctuations)
        R = corrcoef(ri, rj);
        nc_per_stim(p, s) = R(1, 2);

        % Geometric mean of mean firing rates
        geomean_per_stim(p, s) = sqrt(mean_rates_floored(s, ni) * mean_rates_floored(s, nj));
    end

    if mod(p, 1000) == 0
        fprintf('  %d / %d pairs done\n', p, N_pairs);
    end
end
fprintf('Done.\n');

%% Handle NaN noise correlations (constant responses across trials)
nan_count = sum(isnan(nc_per_stim(:)));
if nan_count > 0
    fprintf('Warning: %d NaN noise correlations (constant trial responses). Setting to 0.\n', nan_count);
    nc_per_stim(isnan(nc_per_stim)) = 0;
end

%% ===== STAGE 1 =====
% For each pair: correlation between geo mean and NC across 50 stimuli

geomean_nc_corr = zeros(N_pairs, 1);
for p = 1:N_pairs
    R = corrcoef(geomean_per_stim(p, :), nc_per_stim(p, :));
    geomean_nc_corr(p) = R(1, 2);
end
geomean_nc_corr(isnan(geomean_nc_corr)) = 0;

% --- Stage 1 Plot A: Example scatter plots for a few pairs ---
n_examples = 4;
rng(42);
example_pairs = randsample(N_pairs, n_examples);

fig1a = figure('Position', [100 100 1000 250]);
for k = 1:n_examples
    p = example_pairs(k);
    ni = pair_idx(p, 1);
    nj = pair_idx(p, 2);

    subplot(1, n_examples, k);
    scatter(geomean_per_stim(p, :), nc_per_stim(p, :), 20, 'filled', ...
        'MarkerFaceAlpha', 0.6);
    xlabel('Geometric mean (spk/s)');
    ylabel('Noise correlation');
    title(sprintf('Pair (%d, %d)\nr = %.2f', ni, nj, geomean_nc_corr(p)));
    hold on;
    % Add regression line
    coeffs = polyfit(geomean_per_stim(p, :), nc_per_stim(p, :), 1);
    x_range = xlim;
    x_line = linspace(x_range(1), x_range(2), 100);
    plot(x_line, polyval(coeffs, x_line), 'r-', 'LineWidth', 1.5);
    hold off;
end
sgtitle(sprintf('%s %s %s — Example pairs: Geo Mean vs Noise Correlation', ...
    monkey, vp, upper(rendering)));
saveas(fig1a, fullfile(save_path_figs, [fig_prefix '_stage1_example_scatter.png']));

% --- Stage 1 Plot B: Distribution of geo-mean–NC correlations ---
fig1b = figure('Position', [100 100 500 400]);
histogram(geomean_nc_corr, 50, 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'w');
xlabel('Correlation (geo mean vs NC)');
ylabel('Number of neuron pairs');
title(sprintf('%s %s %s — Distribution of Geo Mean–NC correlations\nMedian = %.3f, Mean = %.3f', ...
    monkey, vp, upper(rendering), median(geomean_nc_corr), mean(geomean_nc_corr)));
xline(0, 'k--', 'LineWidth', 1);
xline(median(geomean_nc_corr), 'r-', 'LineWidth', 1.5);
legend('', 'Zero', 'Median', 'Location', 'best');
saveas(fig1b, fullfile(save_path_figs, [fig_prefix '_stage1_geomean_nc_corr_distribution.png']));

%% ===== STAGE 2 =====
% For each pair: mean and variance of noise correlation across 50 stimuli

mean_nc  = mean(nc_per_stim, 2);   % [N_pairs x 1]
var_nc   = var(nc_per_stim, 0, 2);  % [N_pairs x 1]

% --- Stage 2 Plot: Mean NC vs Variance of NC ---
fig2 = figure('Position', [100 100 550 450]);
scatter(mean_nc, var_nc, 10, 'filled', 'MarkerFaceAlpha', 0.4);
xlabel('Mean noise correlation');
ylabel('Variance of noise correlation');
title(sprintf('%s %s %s — Mean vs Variance of NC across stimuli\n%d neuron pairs', ...
    monkey, vp, upper(rendering), N_pairs));
hold on;
% Add regression line
valid = ~isnan(mean_nc) & ~isnan(var_nc);
coeffs2 = polyfit(mean_nc(valid), var_nc(valid), 1);
x_range2 = xlim;
x_line2 = linspace(x_range2(1), x_range2(2), 100);
plot(x_line2, polyval(coeffs2, x_line2), 'r-', 'LineWidth', 1.5);
hold off;

% Add correlation annotation
R2 = corrcoef(mean_nc(valid), var_nc(valid));
text(0.05, 0.95, sprintf('r = %.3f', R2(1,2)), ...
    'Units', 'normalized', 'FontSize', 12, 'VerticalAlignment', 'top');
saveas(fig2, fullfile(save_path_figs, [fig_prefix '_stage2_mean_vs_variance_nc.png']));

%% Save results
results = struct();
results.monkey = monkey;
results.vp = vp;
results.rendering = rendering;
results.timewindow = timewindow;
results.pair_idx = pair_idx;
results.nc_per_stim = nc_per_stim;
results.geomean_per_stim = geomean_per_stim;
results.geomean_nc_corr = geomean_nc_corr;
results.mean_nc = mean_nc;
results.var_nc = var_nc;

save(fullfile(save_path_data, 'noise_correlation_results.mat'), 'results');
fprintf('Results saved to %s\n', save_path_data);