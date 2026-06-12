%{
Filename: geometric_metrics_pc1_diagnostic.m
Author:   Zitong Wang
Date:     2026-05-26

Description:
    Feasibility diagnostic for the PC1-vs-signal-direction angle analysis.
    That analysis hinges on each stimulus's within-cloud first principal
    axis (PC1) being a well-defined, dominant direction. With only 10 trials
    per stimulus in N-dimensional space (within-cloud covariance rank <= 9),
    PC1 is trustworthy only if the noise is strongly anisotropic (a large
    eigengap, e.g. a shared population-gain mode). If the spectrum is flat,
    PC1 is mostly sampling noise and the downstream angles would be junk.

    This script reports, per stimulus and per rendering condition:
        - PC1 variance fraction : lambda_1 / sum(lambda), from the SVD of
          the centered 10xN cloud.
        - the full normalized eigenvalue spectrum (scree), top R = 9 PCs.

    To make the fraction interpretable it is compared against a within-
    neuron trial-shuffle null: each neuron's 10 values are independently
    permuted across trials, which destroys cross-neuron noise correlations
    while preserving every neuron's marginal variance. The null therefore
    captures the PC1 dominance expected from sampling + variance
    heterogeneity ALONE. Observed PC1 fraction above the null indicates a
    genuine correlated (shared) noise mode -- i.e. a PC1 worth trusting.
    A uniform reference (1/R) marks the fully-isotropic, no-sampling floor.

    Normalization (`normalize_mode`): 'raw' or 'zscore' (per-neuron z-score
    on the 500xN trial matrix, matching the decoding pipeline). 'pt' is a
    configuration-level normalization of the 50 means and has no within-
    cloud meaning, so it is treated as 'raw' here.

Inputs:
    neuronal_data/<monkey>/<vp>/<monkey>_<vp>_allstim.mat
        variable 'three_stim_array' : 1x3 cell {ac, ec, ex} of
                                       [trials x neurons] spike data.

Outputs:
    results/geometric_metrics_outputs/<monkey>/<vp>/
        pc1_diagnostic_<normalize_mode>_results.mat
    results/figures/simple_geometrics/pc1_angle/<monkey>_<vp>/
        pc1_variance_diagnostic_<normalize_mode>.{fig,png}
%}

%% Initialization
clc; clear;

%% Configuration
monkey         = 'KO';        % 'FR' or 'KO'
vp             = 'V1';        % 'V1' or 'V2'
timewindow     = [330 630];   % spike-count window (ms)
normalize_mode = 'zscore';    % 'raw' | 'zscore'  (pt treated as raw)
shuffle_repeat = 200;         % # within-neuron trial shuffles for the null
rng(1);

%% Paths
data_file = fullfile('..','..','neuronal_data', monkey, vp, ...
    sprintf('%s_%s_allstim.mat', monkey, vp));
result_root = fullfile('..','..','results','geometric_metrics_outputs', monkey, vp);
fig_root    = fullfile('..','..','results','figures','simple_geometrics','pc1_angle', ...
    sprintf('%s_%s', monkey, vp));
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
nstim  = max(labels);
ntrial = sum(labels == 1);          % trials per stimulus (10)
N      = size(data_trial{1}, 2);    % # neurons
R      = min(ntrial - 1, N);        % # non-trivial PCs after centering (9)

%% Normalization (trial level)
data_proc = cell(1, nconds);
for i = 1:nconds
    if strcmp(normalize_mode, 'zscore')
        data_proc{i} = zscore(data_trial{i});
    else
        data_proc{i} = data_trial{i};   % raw (and pt fallback)
    end
end

%% Per-stimulus index lists
idx_by_stim = cell(1, nconds);
for i = 1:nconds
    idx_by_stim{i} = cell(nstim, 1);
    for s = 1:nstim
        idx_by_stim{i}{s} = find(labels == s);
    end
end

%% Observed spectra
fprintf('%s %s | normalize=%s : computing per-stimulus spectra\n', monkey, vp, normalize_mode);
pc1_frac = zeros(nstim, nconds);          % lambda_1 / sum(lambda)
spec_obs = zeros(nstim, R, nconds);       % full normalized spectrum (top R)
for i = 1:nconds
    for s = 1:nstim
        Xs = data_proc{i}(idx_by_stim{i}{s}, :);
        frac = cloud_spectrum(Xs, R);
        spec_obs(s,:,i) = frac;
        pc1_frac(s,i)   = frac(1);
    end
end

%% Null spectra (within-neuron trial shuffle)
fprintf('Null (within-neuron shuffle, R=%d)\n', shuffle_repeat);
null_pc1_frac  = zeros(nstim, nconds, shuffle_repeat);
null_spec_sum  = zeros(R, nconds);        % accumulate mean null spectrum
for i = 1:nconds
    for s = 1:nstim
        Xs = data_proc{i}(idx_by_stim{i}{s}, :);
        for r = 1:shuffle_repeat
            Xsh = shuffle_columns(Xs);
            frac = cloud_spectrum(Xsh, R);
            null_pc1_frac(s,i,r) = frac(1);
            null_spec_sum(:,i)   = null_spec_sum(:,i) + frac;
        end
    end
end
null_spec_mean = null_spec_sum / (nstim * shuffle_repeat);   % [R x nconds]

%% Summaries
uniform_ref = 1 / R;
summary = struct();
for i = 1:nconds
    cn = cond_names{i};
    obs_v  = pc1_frac(:,i);
    null_v = reshape(null_pc1_frac(:,i,:), [], 1);
    null_hi_per_stim = prctile(null_pc1_frac(:,i,:), 97.5, 3);   % per-stimulus null ceiling
    summary.(cn).obs_mean      = mean(obs_v);
    summary.(cn).null_mean     = mean(null_v);
    summary.(cn).frac_above    = mean(obs_v > null_hi_per_stim); % frac of stimuli above own null 97.5
    summary.(cn).eigengap_med  = median(spec_obs(:,1,i) ./ max(spec_obs(:,2,i), eps));
end

%% Save
% save(fullfile(result_root, sprintf('pc1_diagnostic_%s_results.mat', normalize_mode)), ...
%     'pc1_frac','spec_obs','null_pc1_frac','null_spec_mean','summary', ...
%     'uniform_ref','normalize_mode','shuffle_repeat','timewindow', ...
%     'cond_names','R','N','ntrial');

%% Visualization
fprintf('Generating figure...\n');
plot_pc1_diagnostic(pc1_frac, spec_obs, null_pc1_frac, null_spec_mean, uniform_ref, ...
    cond_names, normalize_mode, fig_root, sprintf('pc1_variance_diagnostic_%s', normalize_mode));

%% Print summary
fprintf('\n=== %s %s [normalize=%s] PC1 variance fraction (R=%d PCs, uniform=%.3f) ===\n', ...
    monkey, vp, normalize_mode, R, uniform_ref);
fprintf('%-5s | %-12s | %-12s | %-16s | %-12s\n', ...
    'Cond', 'obs mean', 'null mean', 'frac>null97.5', 'med L1/L2');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:nconds
    cn = cond_names{i};
    fprintf('%-5s | %-12.3f | %-12.3f | %-16.2f | %-12.2f\n', ...
        upper(cn), summary.(cn).obs_mean, summary.(cn).null_mean, ...
        summary.(cn).frac_above, summary.(cn).eigengap_med);
end


%% ====================== Local functions ======================

function frac = cloud_spectrum(Xs, R)
% Normalized variance fraction of the top R PCs of a centered cloud.
% Xs: [n x N]. Returns R x 1 (padded with zeros if rank < R).
Xc = Xs - mean(Xs, 1);
sv = svd(Xc, 'econ');
lam = sv.^2;
tot = sum(lam);
frac = zeros(R, 1);
if tot > 0
    m = min(R, numel(lam));
    frac(1:m) = lam(1:m) / tot;
end
end

function Xsh = shuffle_columns(Xs)
% Independently permute each column (neuron) across rows (trials).
% Preserves per-neuron marginal variance, destroys cross-neuron structure.
[n, p] = size(Xs);
[~, ord] = sort(rand(n, p), 1);
lin = ord + (0:p-1) * n;
Xsh = Xs(lin);
end

function plot_pc1_diagnostic(pc1_frac, spec_obs, null_pc1_frac, null_spec_mean, ...
                             uniform_ref, cond_names, normalize_mode, fig_root, fname)
% Row 1: per-stimulus PC1 variance fraction (sorted) vs null band + uniform.
% Row 2: scree (per-PC variance fraction, mean across stimuli) obs vs null.
nconds = numel(cond_names);
R = size(spec_obs, 2);
nstim = size(pc1_frac, 1);

col_obs  = [0.20 0.40 0.80];
col_null = [0.60 0.60 0.60];

fig = figure('Position', [100 100 1500 760], 'Color', 'w');

% ---- Row 1: PC1 fraction per stimulus ----
ymax1 = max(pc1_frac(:)) * 1.1;
for i = 1:nconds
    subplot(2, nconds, i); hold on;
    obs_sorted = sort(pc1_frac(:,i), 'descend');

    null_v  = reshape(null_pc1_frac(:,i,:), [], 1);
    nlo = prctile(null_v, 2.5);
    nhi = prctile(null_v, 97.5);
    nmn = mean(null_v);

    % pooled null band (shaded horizontal region)
    fill([1 nstim nstim 1], [nlo nlo nhi nhi], col_null, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.25, 'HandleVisibility', 'off');
    h_null = plot([1 nstim], [nmn nmn], '-', 'Color', col_null, 'LineWidth', 1.4);
    h_unif = plot([1 nstim], [uniform_ref uniform_ref], 'k:', 'LineWidth', 1.0);
    h_obs  = plot(1:nstim, obs_sorted, '-o', 'Color', col_obs, ...
        'LineWidth', 1.4, 'MarkerSize', 3, 'MarkerFaceColor', col_obs);

    n_above = sum(pc1_frac(:,i) > prctile(null_pc1_frac(:,i,:), 97.5, 3));
    xlim([1 nstim]); ylim([0 ymax1]);
    xlabel('stimulus (sorted)'); ylabel('PC1 variance fraction');
    title(sprintf('%s   (%d/%d > null)', upper(cond_names{i}), n_above, nstim));
    set(gca, 'Box', 'on');
    if i == 1
        legend([h_obs, h_null, h_unif], {'observed', 'null mean', '1/R'}, ...
            'Location', 'northeast', 'Box', 'off');
    end
end

% ---- Row 2: scree (variance fraction per PC) ----
for i = 1:nconds
    subplot(2, nconds, nconds + i); hold on;
    obs_m  = mean(spec_obs(:,:,i), 1);
    obs_lo = prctile(spec_obs(:,:,i), 2.5, 1);
    obs_hi = prctile(spec_obs(:,:,i), 97.5, 1);
    xpc = 1:R;

    fill([xpc, fliplr(xpc)], [obs_lo, fliplr(obs_hi)], col_obs, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.15, 'HandleVisibility', 'off');
    h_obs  = plot(xpc, obs_m, '-o', 'Color', col_obs, 'LineWidth', 1.4, ...
        'MarkerSize', 3, 'MarkerFaceColor', col_obs);
    h_null = plot(xpc, null_spec_mean(:,i), '--s', 'Color', col_null, ...
        'LineWidth', 1.4, 'MarkerSize', 3);
    h_unif = plot([1 R], [uniform_ref uniform_ref], 'k:', 'LineWidth', 1.0);

    xlim([1 R]); ylim([0 1]);
    xticks(1:R);
    xlabel('PC index'); ylabel('variance fraction');
    title(sprintf('%s scree', upper(cond_names{i})));
    set(gca, 'Box', 'on');
    if i == 1
        legend([h_obs, h_null, h_unif], {'observed', 'null', '1/R'}, ...
            'Location', 'northeast', 'Box', 'off');
    end
end

sgtitle(sprintf('Per-stimulus PC1 variance fraction  [%s]', normalize_mode));
% saveas(fig, fullfile(fig_root, [fname '.fig']));
% saveas(fig, fullfile(fig_root, [fname '.png']));
end
