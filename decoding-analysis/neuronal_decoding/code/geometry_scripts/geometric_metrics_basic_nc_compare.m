%{
Filename: geometric_metrics_basic_nc_compare.m
Author:   Zitong Wang
Date:     2026-06-02

Description:
    Noise-correlation control for the geometry analysis. Reproduces the
    distribution figures from geometric_metrics_basic.m (per rendering pair:
    Null | <i> self | <j> self | Obs violins), for BOTH the centroid cosine
    and the pairwise distance, computed twice -- once on the original trials
    and once after breaking noise correlations with the deterministic affine
    trial-cycling shift used in the decoding analyses
    (build_affine_shift_config + apply_trial_perturbation, mode 'affine') --
    and overlays the two for comparison.

    Rationale. The geometry is built from trial-AVERAGED responses, and the
    affine shift only reorders trials within each (neuron, stimulus) block,
    so the full-sample mean (hence the trial-averaged correlation, the black
    triangle) is invariant to machine precision. Any effect of noise
    correlations can therefore appear only in the bootstrap distributions
    (violin position/width), via the joint bootstrap fluctuation of the
    per-neuron mean estimates.

    Both metrics are computed in the same bootstrap pass (shared resamples),
    and the RNG is reset to the same seed before each version, so the
    original and NC-removed runs see identical resample indices, ceiling
    draws, and shuffle permutations -- the only difference is the
    deterministic affine reordering of the underlying trials.

    Color convention: darker shade = original, lighter shade = NC removed.

Inputs:
    neuronal_data/<monkey>/<vp>/<monkey>_<vp>_allstim.mat  (three_stim_array)
    neuronal_data/<monkey>/<vp>/ac/*.mat                   (session bad_channel)

Outputs:
    results/figures/simple_geometrics/nc_control/<monkey>_<vp>/
        cos_distributions_nc_<normalize_mode>.{fig,png}
        dist_distributions_nc_<normalize_mode>.{fig,png}
    results/geometric_metrics_outputs/<monkey>/<vp>/
        nc_compare_<normalize_mode>_results.mat
%}

%% Initialization
clc; clear;

%% Configuration
monkey              = 'FR';        % 'FR' or 'KO'
vp                  = 'V1';        % 'V1' or 'V2'
timewindow          = [330 630];   % spike-count window (ms)
normalize_mode      = 'zscore';    % 'raw' | 'zscore' | 'pt'
trial_sample_repeat = 100;         % observed bootstrap resamples
ceiling_repeat      = 100;         % bootstrap-pair iterations per condition
shuffle_repeat      = 1000;        % shuffle-null iterations
seed                = 1;           % shared seed for the original/NC-removed comparison

%% Paths
data_root  = fullfile('..','..','neuronal_data');
data_file  = fullfile(data_root, monkey, vp, sprintf('%s_%s_allstim.mat', monkey, vp));
result_root = fullfile('..','..','results','geometric_metrics_outputs', monkey, vp);
fig_root    = fullfile('..','..','results','figures','simple_geometrics','nc_control', ...
    sprintf('%s_%s', monkey, vp));
if ~exist(result_root,'dir'); mkdir(result_root); end
if ~exist(fig_root,'dir');    mkdir(fig_root);    end

%% Load data
tmp = load(data_file, 'three_stim_array');
spike_data = tmp.three_stim_array;
nconds = numel(spike_data);

data_trial = cell(1, nconds);
label_trial = cell(1, nconds);
for i = 1:nconds
    [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, timewindow);
end
labels = label_trial{1};
N = size(data_trial{1}, 2);

cond_names = {'ac','ec','ex'};
pair_list  = {[1 2], [1 3], [2 3]};
pair_names = {'acec', 'acex', 'ecex'};

%% Build affine shift config (per-neuron pattern by session membership)
opts = build_affine_shift_config(monkey, vp, N, data_root);

%% Two trial-data versions: original and affine NC-removed
% Affine shift applied to the raw loader output (decoding-pipeline drop-in
% point), then normalization applied to each version.
data_orig = data_trial;
data_pert = cell(1, nconds);
for i = 1:nconds
    data_pert{i} = apply_trial_perturbation(data_trial{i}, labels, 'affine', opts);
end

if strcmp(normalize_mode, 'zscore')
    for i = 1:nconds
        data_orig{i} = zscore(data_orig{i});
        data_pert{i} = zscore(data_pert{i});
    end
end

%% Compute distance + cosine distributions for each version (identical RNG)
params = struct('normalize_mode', normalize_mode, 'R_obs', trial_sample_repeat, ...
    'R_ceil', ceiling_repeat, 'R_null', shuffle_repeat);

fprintf('%s %s | normalize=%s : original\n', monkey, vp, normalize_mode);
O = compute_distributions(data_orig, labels, pair_list, params, seed);
fprintf('%s %s | normalize=%s : NC removed (affine)\n', monkey, vp, normalize_mode);
P = compute_distributions(data_pert, labels, pair_list, params, seed);

%% Save
% save(fullfile(result_root, sprintf('nc_compare_%s_results.mat', normalize_mode)), ...
%     'O','P','normalize_mode','trial_sample_repeat','ceiling_repeat','shuffle_repeat', ...
%     'seed','timewindow','cond_names','pair_names','pair_list');

%% Figures (one per metric)
fprintf('Generating figures...\n');
plot_compare(O.cos, P.cos, pair_list, pair_names, cond_names, ...
    sprintf('Centroid cosine — Pearson r: original vs noise-correlation removed  [%s]', normalize_mode), ...
    fig_root, sprintf('cos_distributions_nc_%s', normalize_mode));

plot_compare(O.dist, P.dist, pair_list, pair_names, cond_names, ...
    sprintf('Distance — Pearson r: original vs noise-correlation removed  [%s]', normalize_mode), ...
    fig_root, sprintf('dist_distributions_nc_%s', normalize_mode));

%% Print summary
fprintf('\n=== %s %s [normalize=%s] cross-rendering Pearson r (orig / NC removed) ===\n', monkey, vp, normalize_mode);
fprintf('%-6s | %-26s | %-26s\n', 'Pair', 'COS obs / trial-avg', 'DIST obs / trial-avg');
for p = 1:numel(pair_list)
    fprintf('%-6s | %.3f,%.3f / %.3f,%.3f | %.3f,%.3f / %.3f,%.3f\n', pair_names{p}, ...
        mean(O.cos.obs(:,p)),  mean(P.cos.obs(:,p)),  O.cos.avg(p),  P.cos.avg(p), ...
        mean(O.dist.obs(:,p)), mean(P.dist.obs(:,p)), O.dist.avg(p), P.dist.avg(p));
end


%% ====================== Local functions ======================

function out = compute_distributions(data_proc, labels, pair_list, params, seed)
% Distance + cosine observed/null/ceiling distributions + trial-averaged
% point estimates, computed in one shared bootstrap pass. RNG reset to
% `seed` so paired calls share an identical random sequence.
rng(seed);
nconds  = numel(data_proc);
n_pairs = numel(pair_list);
nstim   = max(labels);
ptmode  = strcmp(params.normalize_mode, 'pt');

% trial-averaged point estimates (invariant to within-block trial reorder)
D_avg = cell(1, nconds); C_avg = cell(1, nconds);
for i = 1:nconds
    M = trial_mean(data_proc{i}, labels);
    if ptmode, M = pt_normalize(M); end
    [D_avg{i}, C_avg{i}] = geom_metrics(M);
end
dist_avg = zeros(1, n_pairs); cos_avg = zeros(1, n_pairs);
for p = 1:n_pairs
    i = pair_list{p}(1); j = pair_list{p}(2);
    dist_avg(p) = pcorr(D_avg{i}, D_avg{j});
    cos_avg(p)  = pcorr(C_avg{i}, C_avg{j});
end

% observed distribution (bootstrap)
dist_obs = zeros(params.R_obs, n_pairs); cos_obs = zeros(params.R_obs, n_pairs);
for r = 1:params.R_obs
    Dr = cell(1, nconds); Cr = cell(1, nconds);
    for i = 1:nconds
        M = trial_mean_bootstrap(data_proc{i}, labels);
        if ptmode, M = pt_normalize(M); end
        [Dr{i}, Cr{i}] = geom_metrics(M);
    end
    for p = 1:n_pairs
        i = pair_list{p}(1); j = pair_list{p}(2);
        dist_obs(r,p) = pcorr(Dr{i}, Dr{j});
        cos_obs(r,p)  = pcorr(Cr{i}, Cr{j});
    end
end

% ceiling (two independent bootstrap resamples within a condition)
dist_ceil = zeros(params.R_ceil, nconds); cos_ceil = zeros(params.R_ceil, nconds);
for r = 1:params.R_ceil
    for i = 1:nconds
        M1 = trial_mean_bootstrap(data_proc{i}, labels);
        M2 = trial_mean_bootstrap(data_proc{i}, labels);
        if ptmode, M1 = pt_normalize(M1); M2 = pt_normalize(M2); end
        [D1, C1] = geom_metrics(M1);
        [D2, C2] = geom_metrics(M2);
        dist_ceil(r,i) = pcorr(D1, D2);
        cos_ceil(r,i)  = pcorr(C1, C2);
    end
end

% null (stimulus-label shuffle nested in bootstrap; one perm per pair applied to both metrics)
dist_null = zeros(params.R_null, n_pairs); cos_null = zeros(params.R_null, n_pairs);
for r = 1:params.R_null
    Dr = cell(1, nconds); Cr = cell(1, nconds);
    for i = 1:nconds
        M = trial_mean_bootstrap(data_proc{i}, labels);
        if ptmode, M = pt_normalize(M); end
        [Dr{i}, Cr{i}] = geom_metrics(M);
    end
    for p = 1:n_pairs
        i = pair_list{p}(1); j = pair_list{p}(2);
        perm = randperm(nstim);
        dist_null(r,p) = pcorr(Dr{i}, Dr{j}(perm, perm));
        cos_null(r,p)  = pcorr(Cr{i}, Cr{j}(perm, perm));
    end
end

out.cos  = struct('obs', cos_obs,  'null', cos_null,  'ceil', cos_ceil,  'avg', cos_avg);
out.dist = struct('obs', dist_obs, 'null', dist_null, 'ceil', dist_ceil, 'avg', dist_avg);
end

function M = trial_mean(X, labels)
nstim = max(labels);
M = zeros(nstim, size(X,2));
for s = 1:nstim
    M(s,:) = mean(X(labels==s, :), 1);
end
end

function M = trial_mean_bootstrap(X, labels)
nstim = max(labels);
M = zeros(nstim, size(X,2));
for s = 1:nstim
    idx = find(labels == s);
    boot_idx = idx(randi(numel(idx), numel(idx), 1));
    M(s,:) = mean(X(boot_idx, :), 1);
end
end

function M = pt_normalize(M)
c  = mean(M, 1);
M  = M - c;
fn = norm(M, 'fro');
if fn > 0, M = M / fn; end
end

function [D, Cos] = geom_metrics(M)
% D   : pairwise Euclidean distances.
% Cos : centroid-anchored cosine similarity (diagonal 1).
nstim = size(M, 1);
D = squareform(pdist(M, 'euclidean'));
c  = mean(M, 1);
Mc = M - c;
nrm = sqrt(sum(Mc.^2, 2));
U = Mc ./ max(nrm, eps);
Cos = U * U';
Cos = max(min(Cos, 1), -1);
Cos(1:nstim+1:end) = 1;
end

function r = pcorr(M1, M2)
% Pearson r over the 1225 unique upper-triangular off-diagonal pairs.
n = size(M1, 1);
mask = triu(true(n), 1);
r = corr(M1(mask), M2(mask), 'type', 'Pearson');
end

function plot_compare(Om, Pm, pair_list, pair_names, cond_names, ttl, fig_root, fname)
% Per rendering pair: Null | <i> self | <j> self | Obs, original (darker,
% left) vs NC-removed (lighter, right). Om/Pm carry .obs .null .ceil .avg.
% Original and NC-removed share matched bootstrap iterations (same RNG seed),
% so a PAIRED Wilcoxon signed-rank test compares them per category; a
% significance marker is placed in the upper middle of each violin pair.
fig = figure('Position', [100 100 1600 520], 'Color', 'w');
n_pairs = numel(pair_list);

c_null_o = [0.55 0.55 0.55]; c_null_p = [0.82 0.82 0.82];
c_ceil_o = [0.40 0.70 0.50]; c_ceil_p = [0.73 0.89 0.80];
c_obs_o  = [0.25 0.45 0.80]; c_obs_p  = [0.63 0.77 0.94];

dx = 0.19; w = 0.34;
for p = 1:n_pairs
    i = pair_list{p}(1); j = pair_list{p}(2);
    subplot(1, n_pairs, p); hold on;

    catO = {Om.null(:,p), Om.ceil(:,i), Om.ceil(:,j), Om.obs(:,p)};
    catP = {Pm.null(:,p), Pm.ceil(:,i), Pm.ceil(:,j), Pm.obs(:,p)};
    colO = {c_null_o, c_ceil_o, c_ceil_o, c_obs_o};
    colP = {c_null_p, c_ceil_p, c_ceil_p, c_obs_p};

    for c = 1:4
        violin_at(catO{c}, c-dx, w, colO{c});
        violin_at(catP{c}, c+dx, w, colP{c});

        % paired test on matched bootstrap iterations + significance marker
        d = catO{c} - catP{c};
        if all(d == 0)
            pval = 1;
        else
            pval = signrank(catO{c}, catP{c});
        end
        ytop = min(max([catO{c}; catP{c}]) + 0.05, 1.02);
        if pval < 0.05
            text(c, ytop, sig_marker(pval), 'HorizontalAlignment', 'center', ...
                'FontSize', 13, 'FontWeight', 'bold');
        else
            text(c, ytop, 'n.s.', 'HorizontalAlignment', 'center', ...
                'FontSize', 8, 'Color', [0.55 0.55 0.55]);
        end
    end

    % trial-averaged point estimates (coincide: invariant to the shift)
    plot(4-dx, Om.avg(p), 'k^', 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
    plot(4+dx, Pm.avg(p), 'k^', 'MarkerSize', 7, 'MarkerFaceColor', [0.5 0.5 0.5], 'HandleVisibility', 'off');

    yline(0, 'k:', 'LineWidth', 0.5);
    xlim([0.4 4.6]); ylim([-0.3 1.05]);
    xticks(1:4);
    xticklabels({'Null', [upper(cond_names{i}) ' self'], [upper(cond_names{j}) ' self'], 'Obs'});
    ylabel('Pearson r');
    title(sprintf('%s vs %s', upper(cond_names{i}), upper(cond_names{j})));
    set(gca, 'Box', 'on');

    if p == 1
        hO = patch(NaN, NaN, [0.45 0.45 0.45], 'EdgeColor', 'none');
        hP = patch(NaN, NaN, [0.82 0.82 0.82], 'EdgeColor', 'none');
        legend([hO hP], {'original', 'NC removed (affine)'}, 'Location', 'southeast', 'Box', 'off');
    end
end
sgtitle(ttl);
% saveas(fig, fullfile(fig_root, [fname '.fig']));
saveas(fig, fullfile(fig_root, [fname '.png']));
end

function s = sig_marker(p)
% Conventional significance tiers.
if     p < 0.001, s = '***';
elseif p < 0.01,  s = '**';
elseif p < 0.05,  s = '*';
else,             s = 'n.s.';
end
end

function violin_at(data, x_center, width, color)
% Manual violin (kernel density), centered at x_center, with median bar.
data = data(:);
if numel(unique(data)) < 3 || std(data) == 0
    m = mean(data);
    plot([x_center - width/2, x_center + width/2], [m, m], '-', ...
        'Color', color, 'LineWidth', 2, 'HandleVisibility', 'off');
    return;
end
[f, xi] = ksdensity(data);
f = f / max(f) * (width / 2);
patch([x_center + f, fliplr(x_center - f)], [xi, fliplr(xi)], color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8, 'HandleVisibility', 'off');
m = median(data);
plot([x_center - width/3, x_center + width/3], [m, m], 'k-', 'LineWidth', 1.1, 'HandleVisibility', 'off');
end
