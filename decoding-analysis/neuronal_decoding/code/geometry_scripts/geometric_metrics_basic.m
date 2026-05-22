%{
Filename: geometric_metrics_basic.m
Author:   Zitong Wang
Date:     2026-05-21

Description:
    Computes simple geometric metrics on trial-averaged neural population
    responses to directly test Procrustes-invariant quantities (relative
    distances and angles), addressing reviewer 1's concern that the
    "shape" of the manifold can be assessed without training a decoder.

    For each rendering condition (ac, ec, ex), produces:
        - D   : 50x50 pairwise Euclidean distance matrix.
        - Cos : 50x50 centroid-anchored cosine similarity matrix.
                Cos(i,j) = cos(theta_C), where theta_C is the angle at the
                population centroid C in the triangle (stim_i, stim_j, C).
                Equivalently, Cos(i,j) is the cosine of the centered
                stimulus mean vectors (i-C) and (j-C). Diagonal is 1.

    Two estimation modes are computed in every run:
        - trial_averaged  : single 50xN mean per condition.
        - trial_resampled : bootstrap 10 trials with replacement per
                            stimulus, average the 50x50 metric across R
                            resamples. (Default visualization target.)

    Three normalization modes, switched via `normalize_mode`:
        - 'raw'    : no normalization.
        - 'zscore' : per-neuron z-score on the trial-level 500xN matrix
                     (matches procrustes_decoding_basic.m convention).
        - 'pt'     : Procrustes-style preprocessing on the 50xN stimulus
                     means -- subtract centroid, divide by Frobenius norm.

    Three distributions per (rendering pair, metric, statistic):
        - Observed : per-resample cross-condition correlation, R values,
                     drawn from the same bootstrap iterations that populate
                     the resampled heatmaps (no extra cost).
        - Null     : stimulus-label shuffle nested in bootstrap. At each
                     iteration, a fresh joint bootstrap resample is drawn
                     for all conditions, then stimulus IDs in one condition
                     of each pair are permuted before correlating. This
                     matches the noise structure of the observed and
                     ceiling distributions, so all three violins are
                     directly comparable. Centered near 0 under no
                     correspondence.
        - Ceiling  : two independent bootstrap resamples WITHIN a single
                     condition; correlate their distance/cosine matrices.
                     One ceiling distribution per condition (AC, EC, EX),
                     since SNR differs across renderings.

    Reports Pearson r and Spearman rho for all three distributions, plus
    point estimates from the trial-averaged matrices.

Inputs:
    neuronal_data/<monkey>/<vp>/<monkey>_<vp>_allstim.mat
        variable 'three_stim_array' : 1x3 cell {ac, ec, ex} of
                                       [trials x neurons] spike data.

Outputs:
    results/geometric_metrics_outputs/<monkey>/<vp>/geom_<normalize_mode>_results.mat
    results/figures/geometric_metrics_<monkey>_<vp>/*.fig + *.png
        (saveas calls are commented out during active development.)
%}

%% Initialization
clc; clear;
close all
%% Configuration
monkey              = 'KO';        % 'FR' or 'KO'
vp                  = 'V2';        % 'V1' or 'V2'
timewindow          = [330 630];   % spike-count window (ms)
normalize_mode      = 'raw';        % 'raw' | 'zscore' | 'pt'
trial_sample_repeat = 100;         % # bootstrap resamples for observed distribution
ceiling_repeat      = 100;         % # bootstrap-pair iterations per condition for ceiling
shuffle_repeat      = 1000;        % # stimulus-label shuffles for null
rng(1);

%% Paths
data_file = fullfile('..','..','neuronal_data', monkey, vp, ...
    sprintf('%s_%s_allstim.mat', monkey, vp));
result_root = fullfile('..','..','results','geometric_metrics_outputs', monkey, vp);
fig_root    = fullfile('..','..','results','figures', ...
    sprintf('geometric_metrics_%s_%s', monkey, vp));
if ~exist(result_root,'dir'); mkdir(result_root); end
if ~exist(fig_root,'dir');    mkdir(fig_root);    end

%% Load data
tmp = load(data_file, 'three_stim_array');
spike_data = tmp.three_stim_array;
nconds = numel(spike_data);

data_trial  = cell(1, nconds);
label_trial = cell(1, nconds);
for i = 1:nconds
    [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, timewindow);
end
labels = label_trial{1};

cond_names = {'ac','ec','ex'};
nstim = 50;

%% Trial-level normalization (only applied for z-score mode)
data_trial_proc = cell(1, nconds);
for i = 1:nconds
    if strcmp(normalize_mode, 'zscore')
        data_trial_proc{i} = zscore(data_trial{i});
    else
        data_trial_proc{i} = data_trial{i};
    end
end

%% Rendering pairs
pair_list  = {[1 2], [1 3], [2 3]};
pair_names = {'acec', 'acex', 'ecex'};
n_pairs    = numel(pair_list);

%% Mode A: trial-averaged matrices and point-estimate correlations
fprintf('Trial-averaged mode (normalize=%s)\n', normalize_mode);
D_avg   = cell(1, nconds);
Cos_avg = cell(1, nconds);
for i = 1:nconds
    M = trial_mean(data_trial_proc{i}, labels);
    if strcmp(normalize_mode, 'pt'), M = pt_normalize(M); end
    [D_avg{i}, Cos_avg{i}] = geom_metrics(M);
end

stats_avg = struct();
for p = 1:n_pairs
    i = pair_list{p}(1); j = pair_list{p}(2); pn = pair_names{p};
    [stats_avg.(pn).dist_r, stats_avg.(pn).dist_rho] = pair_corr(D_avg{i},   D_avg{j});
    [stats_avg.(pn).cos_r,  stats_avg.(pn).cos_rho]  = pair_corr(Cos_avg{i}, Cos_avg{j});
end

%% Mode B: trial-resampled (bootstrap) — averaged matrices + per-resample observed distribution
fprintf('Trial-resampled mode (normalize=%s, R=%d)\n', normalize_mode, trial_sample_repeat);

D_rs_sum   = cell(1, nconds);
Cos_rs_sum = cell(1, nconds);
for i = 1:nconds
    D_rs_sum{i}   = zeros(nstim, nstim);
    Cos_rs_sum{i} = zeros(nstim, nstim);
end

dist_r_obs   = zeros(trial_sample_repeat, n_pairs);
dist_rho_obs = zeros(trial_sample_repeat, n_pairs);
cos_r_obs    = zeros(trial_sample_repeat, n_pairs);
cos_rho_obs  = zeros(trial_sample_repeat, n_pairs);

for r = 1:trial_sample_repeat
    D_r   = cell(1, nconds);
    Cos_r = cell(1, nconds);
    for i = 1:nconds
        M = trial_mean_bootstrap(data_trial_proc{i}, labels);
        if strcmp(normalize_mode, 'pt'), M = pt_normalize(M); end
        [D_r{i}, Cos_r{i}] = geom_metrics(M);
        D_rs_sum{i}   = D_rs_sum{i}   + D_r{i};
        Cos_rs_sum{i} = Cos_rs_sum{i} + Cos_r{i};
    end
    % per-resample cross-condition correlations (observed distribution)
    for p = 1:n_pairs
        i = pair_list{p}(1); j = pair_list{p}(2);
        [dist_r_obs(r,p), dist_rho_obs(r,p)] = pair_corr(D_r{i},   D_r{j});
        [cos_r_obs(r,p),  cos_rho_obs(r,p)]  = pair_corr(Cos_r{i}, Cos_r{j});
    end
end
D_rs   = cellfun(@(x) x/trial_sample_repeat, D_rs_sum,   'UniformOutput', false);
Cos_rs = cellfun(@(x) x/trial_sample_repeat, Cos_rs_sum, 'UniformOutput', false);

% point-estimate correlations on the resample-averaged matrices (for the heatmap caption)
stats_rs = struct();
for p = 1:n_pairs
    i = pair_list{p}(1); j = pair_list{p}(2); pn = pair_names{p};
    [stats_rs.(pn).dist_r, stats_rs.(pn).dist_rho] = pair_corr(D_rs{i},   D_rs{j});
    [stats_rs.(pn).cos_r,  stats_rs.(pn).cos_rho]  = pair_corr(Cos_rs{i}, Cos_rs{j});
end

%% Ceiling: bootstrap-pair within each condition
fprintf('Ceiling (bootstrap-pair, R=%d per condition)\n', ceiling_repeat);
dist_r_ceil   = zeros(ceiling_repeat, nconds);
dist_rho_ceil = zeros(ceiling_repeat, nconds);
cos_r_ceil    = zeros(ceiling_repeat, nconds);
cos_rho_ceil  = zeros(ceiling_repeat, nconds);

for r = 1:ceiling_repeat
    for i = 1:nconds
        M1 = trial_mean_bootstrap(data_trial_proc{i}, labels);
        M2 = trial_mean_bootstrap(data_trial_proc{i}, labels);
        if strcmp(normalize_mode, 'pt')
            M1 = pt_normalize(M1); M2 = pt_normalize(M2);
        end
        [D1, C1] = geom_metrics(M1);
        [D2, C2] = geom_metrics(M2);
        [dist_r_ceil(r,i), dist_rho_ceil(r,i)] = pair_corr(D1, D2);
        [cos_r_ceil(r,i),  cos_rho_ceil(r,i)]  = pair_corr(C1, C2);
    end
end

%% Null: stimulus-label shuffle nested in bootstrap resamples
% Each iteration: fresh joint bootstrap resample of all conditions, then
% permute stimulus IDs in one condition of each pair before correlating.
% This matches the noise structure of the observed and ceiling
% distributions, so the three violins are directly comparable -- the only
% thing distinguishing them is whether i<->i correspondence is preserved,
% broken via permutation, or replaced by a second independent draw of the
% same condition.
fprintf('Null (shuffle nested in bootstrap, R=%d)\n', shuffle_repeat);
dist_r_null   = zeros(shuffle_repeat, n_pairs);
dist_rho_null = zeros(shuffle_repeat, n_pairs);
cos_r_null    = zeros(shuffle_repeat, n_pairs);
cos_rho_null  = zeros(shuffle_repeat, n_pairs);

for r = 1:shuffle_repeat
    % Joint bootstrap resample of all conditions
    D_r   = cell(1, nconds);
    Cos_r = cell(1, nconds);
    for i = 1:nconds
        M = trial_mean_bootstrap(data_trial_proc{i}, labels);
        if strcmp(normalize_mode, 'pt'), M = pt_normalize(M); end
        [D_r{i}, Cos_r{i}] = geom_metrics(M);
    end
    % Shuffle stimulus IDs in the second condition of each pair before correlating
    for p = 1:n_pairs
        i = pair_list{p}(1); j = pair_list{p}(2);
        perm = randperm(nstim);
        Dj_sh = D_r{j}(perm, perm);
        Cj_sh = Cos_r{j}(perm, perm);
        [dist_r_null(r,p), dist_rho_null(r,p)] = pair_corr(D_r{i},   Dj_sh);
        [cos_r_null(r,p),  cos_rho_null(r,p)]  = pair_corr(Cos_r{i}, Cj_sh);
    end
end

%% Distribution summary stats (mean and 95% percentile CI)
dist_stats = struct();
for p = 1:n_pairs
    pn = pair_names{p};
    dist_stats.(pn).obs.dist_r   = summarize(dist_r_obs(:,p));
    dist_stats.(pn).obs.dist_rho = summarize(dist_rho_obs(:,p));
    dist_stats.(pn).obs.cos_r    = summarize(cos_r_obs(:,p));
    dist_stats.(pn).obs.cos_rho  = summarize(cos_rho_obs(:,p));
    dist_stats.(pn).null.dist_r   = summarize(dist_r_null(:,p));
    dist_stats.(pn).null.dist_rho = summarize(dist_rho_null(:,p));
    dist_stats.(pn).null.cos_r    = summarize(cos_r_null(:,p));
    dist_stats.(pn).null.cos_rho  = summarize(cos_rho_null(:,p));
end
for i = 1:nconds
    cn = cond_names{i};
    dist_stats.ceil.(cn).dist_r   = summarize(dist_r_ceil(:,i));
    dist_stats.ceil.(cn).dist_rho = summarize(dist_rho_ceil(:,i));
    dist_stats.ceil.(cn).cos_r    = summarize(cos_r_ceil(:,i));
    dist_stats.ceil.(cn).cos_rho  = summarize(cos_rho_ceil(:,i));
end

%% Save results
save(fullfile(result_root, sprintf('geom_%s_results.mat', normalize_mode)), ...
    'D_avg','Cos_avg','D_rs','Cos_rs', ...
    'stats_avg','stats_rs','dist_stats', ...
    'dist_r_obs','dist_rho_obs','cos_r_obs','cos_rho_obs', ...
    'dist_r_null','dist_rho_null','cos_r_null','cos_rho_null', ...
    'dist_r_ceil','dist_rho_ceil','cos_r_ceil','cos_rho_ceil', ...
    'normalize_mode','trial_sample_repeat','ceiling_repeat','shuffle_repeat', ...
    'timewindow','cond_names','pair_names','pair_list');

%% Visualization — heatmaps and scatter (unchanged)
fprintf('Generating figures...\n');
plot_heatmaps(D_rs,    'Distance (resampled)',         cond_names, fig_root, sprintf('dist_heatmap_rs_%s',  normalize_mode), 'linear');
plot_heatmaps(Cos_rs,  'Centroid cosine (resampled)',  cond_names, fig_root, sprintf('cos_heatmap_rs_%s',   normalize_mode), 'cosine');
plot_heatmaps(D_avg,   'Distance (averaged)',          cond_names, fig_root, sprintf('dist_heatmap_avg_%s', normalize_mode), 'linear');
plot_heatmaps(Cos_avg, 'Centroid cosine (averaged)',   cond_names, fig_root, sprintf('cos_heatmap_avg_%s',  normalize_mode), 'cosine');

plot_scatter(D_rs,    'Distance (resampled)',          cond_names, pair_list, pair_names, stats_rs,  'dist', fig_root, sprintf('dist_scatter_rs_%s',  normalize_mode));
plot_scatter(Cos_rs,  'Centroid cosine (resampled)',   cond_names, pair_list, pair_names, stats_rs,  'cos',  fig_root, sprintf('cos_scatter_rs_%s',   normalize_mode));
plot_scatter(D_avg,   'Distance (averaged)',           cond_names, pair_list, pair_names, stats_avg, 'dist', fig_root, sprintf('dist_scatter_avg_%s', normalize_mode));
plot_scatter(Cos_avg, 'Centroid cosine (averaged)',    cond_names, pair_list, pair_names, stats_avg, 'cos',  fig_root, sprintf('cos_scatter_avg_%s',  normalize_mode));

%% Visualization — distribution figures (null / ceiling / observed)
plot_distributions(dist_r_obs, dist_r_null, dist_r_ceil, ...
    pair_list, pair_names, cond_names, stats_avg, 'dist', ...
    'Distance — Pearson r', fig_root, sprintf('dist_distributions_%s', normalize_mode));

plot_distributions(cos_r_obs, cos_r_null, cos_r_ceil, ...
    pair_list, pair_names, cond_names, stats_avg, 'cos', ...
    'Centroid cosine — Pearson r', fig_root, sprintf('cos_distributions_%s', normalize_mode));

%% Print summary
fprintf('\n=== %s %s  [normalize=%s, R_obs=%d, R_ceil=%d, R_null=%d] ===\n', ...
    monkey, vp, normalize_mode, trial_sample_repeat, ceiling_repeat, shuffle_repeat);

fprintf('\n-- Trial-averaged point estimates --\n');
fprintf('%-6s | %-26s | %-26s\n', 'Pair', 'Distance r/rho', 'Cosine r/rho');
for p = 1:n_pairs
    pn = pair_names{p};
    fprintf('%-6s | r=%+.3f  rho=%+.3f     | r=%+.3f  rho=%+.3f\n', ...
        pn, stats_avg.(pn).dist_r, stats_avg.(pn).dist_rho, ...
            stats_avg.(pn).cos_r,  stats_avg.(pn).cos_rho);
end

fprintf('\n-- Observed (trial-resampled) Pearson r: mean [2.5%%, 97.5%%] --\n');
for p = 1:n_pairs
    pn = pair_names{p};
    s = dist_stats.(pn).obs.dist_r; t = dist_stats.(pn).obs.cos_r;
    fprintf('%-6s | dist  %+.3f [%+.3f, %+.3f]   | cos  %+.3f [%+.3f, %+.3f]\n', ...
        pn, s.mean, s.ci(1), s.ci(2), t.mean, t.ci(1), t.ci(2));
end

fprintf('\n-- Null (shuffle) Pearson r: mean [2.5%%, 97.5%%] --\n');
for p = 1:n_pairs
    pn = pair_names{p};
    s = dist_stats.(pn).null.dist_r; t = dist_stats.(pn).null.cos_r;
    fprintf('%-6s | dist  %+.3f [%+.3f, %+.3f]   | cos  %+.3f [%+.3f, %+.3f]\n', ...
        pn, s.mean, s.ci(1), s.ci(2), t.mean, t.ci(1), t.ci(2));
end

fprintf('\n-- Ceiling (bootstrap-pair) Pearson r: mean [2.5%%, 97.5%%] --\n');
for i = 1:nconds
    cn = cond_names{i};
    s = dist_stats.ceil.(cn).dist_r; t = dist_stats.ceil.(cn).cos_r;
    fprintf('%-6s | dist  %+.3f [%+.3f, %+.3f]   | cos  %+.3f [%+.3f, %+.3f]\n', ...
        upper(cn), s.mean, s.ci(1), s.ci(2), t.mean, t.ci(1), t.ci(2));
end


%% ====================== Local functions ======================

function M = trial_mean(X, labels)
% Trial-averaged mean response per stimulus.
nstim = max(labels);
M = zeros(nstim, size(X,2));
for s = 1:nstim
    M(s,:) = mean(X(labels==s, :), 1);
end
end

function M = trial_mean_bootstrap(X, labels)
% Bootstrap-mean response per stimulus: sample n_trials with replacement.
nstim = max(labels);
M = zeros(nstim, size(X,2));
for s = 1:nstim
    idx = find(labels == s);
    boot_idx = idx(randi(numel(idx), numel(idx), 1));
    M(s,:) = mean(X(boot_idx, :), 1);
end
end

function M = pt_normalize(M)
% Procrustes-style normalization: center on stimulus centroid, divide by
% Frobenius norm of the centered matrix.
c  = mean(M, 1);
M  = M - c;
fn = norm(M, 'fro');
if fn > 0, M = M / fn; end
end

function [D, Cos] = geom_metrics(M)
% D   : pairwise Euclidean distances.
% Cos : centroid-anchored cosine similarity.
nstim = size(M, 1);
D = squareform(pdist(M, 'euclidean'));

c  = mean(M, 1);
Mc = M - c;
nrm = sqrt(sum(Mc.^2, 2));
Mc_unit = Mc ./ max(nrm, eps);
Cos = Mc_unit * Mc_unit';
Cos = max(min(Cos, 1), -1);
Cos(1:nstim+1:end) = 1;
end

function [r, rho] = pair_corr(M1, M2)
% Correlation across upper-triangular off-diagonal entries (1225 pairs).
n = size(M1, 1);
mask = triu(true(n), 1);
v1 = M1(mask); v2 = M2(mask);
r   = corr(v1, v2, 'type', 'Pearson');
rho = corr(v1, v2, 'type', 'Spearman');
end

function s = summarize(v)
% Return mean and 2.5/97.5 percentile interval of a distribution.
s.mean = mean(v);
s.ci   = prctile(v, [2.5, 97.5]);
end

function plot_heatmaps(M_cell, ttl, cond_names, fig_root, fname, scale_kind)
% scale_kind : 'linear' uses data-driven shared limits; 'cosine' uses [-1, 1].
fig = figure('Position', [100 100 1500 450], 'Color', 'w');
nconds = numel(M_cell);
if strcmp(scale_kind, 'cosine')
    clim = [-1, 1];
else
    all_vals = vertcat(M_cell{:});
    clim = [min(all_vals(:)), max(all_vals(:))];
    if diff(clim) == 0, clim = clim + [-eps, eps]; end
end
for i = 1:nconds
    subplot(1, nconds, i);
    imagesc(M_cell{i}, clim);
    axis square; colorbar;
    title(upper(cond_names{i}));
    xlabel('stimulus'); ylabel('stimulus');
end
sgtitle(ttl);
% saveas(fig, fullfile(fig_root, [fname '.fig']));
saveas(fig, fullfile(fig_root, [fname '.png']));
end

function plot_scatter(M_cell, ttl, cond_names, pair_list, pair_names, stats, field, fig_root, fname)
fig = figure('Position', [100 100 1500 450], 'Color', 'w');
nstim = size(M_cell{1}, 1);
mask = triu(true(nstim), 1);
for p = 1:numel(pair_list)
    i = pair_list{p}(1); j = pair_list{p}(2); pn = pair_names{p};
    v1 = M_cell{i}(mask); v2 = M_cell{j}(mask);
    subplot(1, numel(pair_list), p);
    scatter(v1, v2, 12, 'filled', 'MarkerFaceAlpha', 0.4);
    axis square; hold on;
    mn = min([v1; v2]); mx = max([v1; v2]);
    plot([mn mx], [mn mx], 'k--');
    xlabel(upper(cond_names{i})); ylabel(upper(cond_names{j}));
    if strcmp(field, 'dist')
        rr = stats.(pn).dist_r; rho = stats.(pn).dist_rho;
    else
        rr = stats.(pn).cos_r;  rho = stats.(pn).cos_rho;
    end
    title(sprintf('%s vs %s : r=%.3f, \\rho=%.3f', ...
        upper(cond_names{i}), upper(cond_names{j}), rr, rho));
end
sgtitle(ttl);
% saveas(fig, fullfile(fig_root, [fname '.fig']));
saveas(fig, fullfile(fig_root, [fname '.png']));
end

function plot_distributions(obs, null_dist, ceil_dist, pair_list, pair_names, ...
                            cond_names, stats_avg, field, ttl, fig_root, fname)
% obs       : [trial_sample_repeat x n_pairs]   per-resample cross-condition r
% null_dist : [shuffle_repeat x n_pairs]        shuffle-null r
% ceil_dist : [ceiling_repeat x nconds]         bootstrap-pair r within each condition
% stats_avg : struct of point estimates from trial-averaged matrices (marker overlay)
% field     : 'dist' or 'cos' (selects the right point-estimate field)

fig = figure('Position', [100 100 1500 500], 'Color', 'w');
n_pairs = numel(pair_list);

col_null = [0.70 0.70 0.70];
col_ceil = [0.55 0.78 0.65];
col_obs  = [0.30 0.50 0.85];

for p = 1:n_pairs
    i = pair_list{p}(1); j = pair_list{p}(2); pn = pair_names{p};
    subplot(1, n_pairs, p);
    hold on;

    width = 0.65;
    violin_at(null_dist(:,p),  1, width, col_null);
    violin_at(ceil_dist(:,i),  2, width, col_ceil);
    violin_at(ceil_dist(:,j),  3, width, col_ceil);
    violin_at(obs(:,p),        4, width, col_obs);

    % point-estimate marker on observed
    if strcmp(field, 'dist'), avg_val = stats_avg.(pn).dist_r;
    else,                     avg_val = stats_avg.(pn).cos_r;
    end
    plot(4, avg_val, 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

    yline(0, 'k:', 'LineWidth', 0.5);
    xlim([0.4, 4.6]);
    ylim([-0.3, 1.05]);
    xticks(1:4);
    xticklabels({'Null', [upper(cond_names{i}) ' self'], ...
                 [upper(cond_names{j}) ' self'], 'Obs'});
    ylabel('Pearson r');
    title(sprintf('%s vs %s', upper(cond_names{i}), upper(cond_names{j})));
    set(gca, 'Box', 'on');
end
sgtitle(ttl);
% saveas(fig, fullfile(fig_root, [fname '.fig']));
saveas(fig, fullfile(fig_root, [fname '.png']));
end

function violin_at(data, x_center, width, color)
% Manual violin: kernel density estimate of `data`, mirrored, centered at
% x_center, max half-width = width/2. Adds a median line.
data = data(:);
if numel(unique(data)) < 3 || std(data) == 0
    m = mean(data);
    plot([x_center - width/2, x_center + width/2], [m, m], '-', ...
        'Color', color, 'LineWidth', 2);
    return;
end
[f, xi] = ksdensity(data);
f = f / max(f) * (width / 2);
patch([x_center + f, fliplr(x_center - f)], [xi, fliplr(xi)], color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.75);
m = median(data);
plot([x_center - width/3, x_center + width/3], [m, m], 'k-', 'LineWidth', 1.2);
end