%{
Filename: geometric_metrics_pc1_angle.m
Author:   Zitong Wang
Date:     2026-05-27

Description:
    Computes the angle between each stimulus's within-cloud first principal
    axis (PC1) and the signal directions to other stimuli, then asks:

      (a) WITHIN a rendering -- what is the distribution of these angles?
          Do clouds preferentially align with signal directions
          (information-limiting noise), or sit mostly orthogonal? Reference:
          in high-D two unrelated directions are nearly orthogonal, so a
          random-orientation null concentrates near 90 degrees; systematic
          deviation toward 0 degrees indicates alignment.

      (b) ACROSS renderings -- is the per-pair angle preserved across
          AC/EC/EX? This is the cue-invariance analog of the centroid-cosine
          analysis, applied to the noise-signal geometry.

    Quantity. For stimuli i and j:
        alignment A(i,j) = |cos(theta)| = |<PC1_i, u_ij>|,
        where u_ij = (mu_j - mu_i)/||mu_j - mu_i|| and PC1_i is the top
        principal axis of cloud i. PC1 has an arbitrary sign, so the acute
        angle is used (alignment in [0,1]; angle = acos(alignment) in
        [0,90] deg). The matrix is ASYMMETRIC -- A(i,j) uses cloud i's PC1
        and the i->j direction -- so all 2450 ordered off-diagonal entries
        are used.

    IMPORTANT -- spurious-alignment bias and the split-half fix.
        If PC1_i and mu_i are estimated from the SAME trials, the mean's
        estimation error tends to lie along the high-variance axis (~PC1_i),
        and since u_ij contains -mu_i, the signal direction is pulled along
        +/-PC1_i, biasing the alignment UPWARD -- i.e. in the direction of
        the (a) effect. The random-orientation null does not capture this.
        Default `use_split` = true therefore estimates PC1_i from one
        disjoint half of the trials and mu_i from the other half, removing
        the bias (at the cost of 5-trial estimates). Set `use_split` = false
        to reproduce the biased all-trials version for comparison.

    PC1 is noisy here (10 trials per cloud; see
    geometric_metrics_pc1_diagnostic.m), so this is exploratory and the
    nulls are essential for interpretation.

    Nulls.
        (a) random-orientation: replace each PC1 with a random unit vector
            in N-D, keep the real signal directions. Concentrates near 90.
        (b) stimulus-correspondence shuffle: permute stimulus IDs of one
            rendering's matrix before correlating, nested in the split
            iterations so the noise structure matches the observed.

    Normalization (`normalize_mode`): 'raw' or 'zscore' (per-neuron z-score
    on the 500xN trial matrix). 'pt' is a configuration-level op with no
    within-cloud meaning and is treated as 'raw'.

Inputs:
    neuronal_data/<monkey>/<vp>/<monkey>_<vp>_allstim.mat
        variable 'three_stim_array' : 1x3 cell {ac, ec, ex} of
                                       [trials x neurons] spike data.

Outputs:
    results/geometric_metrics_outputs/<monkey>/<vp>/
        pc1_angle_<normalize_mode>_<split_tag>_results.mat
    results/figures/simple_geometrics/pc1_angle/<monkey>_<vp>/
        angle_within_<normalize_mode>_<split_tag>.{fig,png}
        angle_across_<normalize_mode>_<split_tag>.{fig,png}
%}

%% Initialization
clc; clear;
close all
%% Configuration
monkey         = 'FR';        % 'FR' or 'KO'
vp             = 'V2';        % 'V1' or 'V2'
timewindow     = [330 630];   % spike-count window (ms)
normalize_mode = 'raw';    % 'raw' | 'zscore'
use_split      = false;        % true: PC1 & mean from disjoint halves (bias-free); false: all trials (biased)
n_iter         = 100;         % # split (or bootstrap) iterations
n_null         = 5;           % # shuffles per iteration for (b) cross-rendering null
rng(1);

if use_split, split_tag = 'split'; else, split_tag = 'full'; end

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
N      = size(data_trial{1}, 2);

%% Normalization (trial level)
data_proc = cell(1, nconds);
for i = 1:nconds
    if strcmp(normalize_mode, 'zscore')
        data_proc{i} = zscore(data_trial{i});
    else
        data_proc{i} = data_trial{i};
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

%% Rendering pairs and off-diagonal mask
pair_list  = {[1 2], [1 3], [2 3]};
pair_names = {'acec', 'acex', 'ecex'};
n_pairs    = numel(pair_list);
offmask    = ~eye(nstim);

%% Iterate
fprintf('%s %s | normalize=%s | split=%d | n_iter=%d\n', ...
    monkey, vp, normalize_mode, use_split, n_iter);

A_sum = cell(1, nconds);   % mean observed alignment matrix per rendering
Nl_sum = cell(1, nconds);  % mean random-orientation null matrix per rendering
for i = 1:nconds
    A_sum{i}  = zeros(nstim);
    Nl_sum{i} = zeros(nstim);
end

obs_meanalign  = zeros(n_iter, nconds);   % (a) mean |cos| over pairs, observed
null_meanalign = zeros(n_iter, nconds);   % (a) mean |cos| over pairs, random-orientation null

obs_corr_r   = zeros(n_iter, n_pairs);    % (b) cross-rendering Pearson, observed
obs_corr_rho = zeros(n_iter, n_pairs);
null_corr_r  = zeros(n_iter, n_null, n_pairs);   % (b) shuffle null Pearson
null_corr_rho= zeros(n_iter, n_null, n_pairs);

for it = 1:n_iter
    A_it = cell(1, nconds);
    for i = 1:nconds
        PC1 = zeros(nstim, N);
        M   = zeros(nstim, N);
        for s = 1:nstim
            Xs = data_proc{i}(idx_by_stim{i}{s}, :);
            [PC1(s,:), M(s,:)] = cloud_pc1_mean(Xs, use_split);
        end
        A = alignment_matrix(PC1, M, offmask);
        A_it{i}   = A;
        A_sum{i}  = A_sum{i} + A;

        % (a) random-orientation null: random unit PC1, real signal dirs
        Rdir = randn(nstim, N);
        Rdir = Rdir ./ sqrt(sum(Rdir.^2, 2));
        Nmat = alignment_matrix(Rdir, M, offmask);
        Nl_sum{i} = Nl_sum{i} + Nmat;

        obs_meanalign(it,i)  = mean(A(offmask));
        null_meanalign(it,i) = mean(Nmat(offmask));
    end

    % (b) cross-rendering correlations + nested shuffle null
    for p = 1:n_pairs
        i = pair_list{p}(1); j = pair_list{p}(2);
        vi = A_it{i}(offmask);
        vj = A_it{j}(offmask);
        obs_corr_r(it,p)   = corr(vi, vj, 'type', 'Pearson');
        obs_corr_rho(it,p) = corr(vi, vj, 'type', 'Spearman');
        for sh = 1:n_null
            perm = randperm(nstim);
            Aj_sh = A_it{j}(perm, perm);
            vjs = Aj_sh(offmask);
            null_corr_r(it,sh,p)   = corr(vi, vjs, 'type', 'Pearson');
            null_corr_rho(it,sh,p) = corr(vi, vjs, 'type', 'Spearman');
        end
    end

    if mod(it, 20) == 0, fprintf('  iter %d/%d\n', it, n_iter); end
end

A_mean  = cellfun(@(x) x / n_iter, A_sum,  'UniformOutput', false);
Nl_mean = cellfun(@(x) x / n_iter, Nl_sum, 'UniformOutput', false);

%% Summaries
within_summary = struct();
for i = 1:nconds
    cn = cond_names{i};
    within_summary.(cn).obs_mean_cos   = mean(obs_meanalign(:,i));
    within_summary.(cn).null_mean_cos  = mean(null_meanalign(:,i));
    within_summary.(cn).obs_mean_ang   = acosd(min(max(mean(obs_meanalign(:,i)),0),1));
    within_summary.(cn).null_mean_ang  = acosd(min(max(mean(null_meanalign(:,i)),0),1));
    within_summary.(cn).obs_ci         = prctile(obs_meanalign(:,i), [2.5 97.5]);
end

across_summary = struct();
for p = 1:n_pairs
    pn = pair_names{p};
    nr  = reshape(null_corr_r(:,:,p), [], 1);
    nrh = reshape(null_corr_rho(:,:,p), [], 1);
    across_summary.(pn).obs_r_mean   = mean(obs_corr_r(:,p));
    across_summary.(pn).obs_r_ci     = prctile(obs_corr_r(:,p), [2.5 97.5]);
    across_summary.(pn).obs_rho_mean = mean(obs_corr_rho(:,p));
    across_summary.(pn).null_r_mean  = mean(nr);
    across_summary.(pn).null_r_ci    = prctile(nr, [2.5 97.5]);
    across_summary.(pn).null_rho_mean= mean(nrh);
end

%% Save
% save(fullfile(result_root, sprintf('pc1_angle_%s_%s_results.mat', normalize_mode, split_tag)), ...
%     'A_mean','Nl_mean','obs_meanalign','null_meanalign', ...
%     'obs_corr_r','obs_corr_rho','null_corr_r','null_corr_rho', ...
%     'within_summary','across_summary', ...
%     'normalize_mode','use_split','split_tag','n_iter','n_null','timewindow', ...
%     'cond_names','pair_names','pair_list','nstim','N');

%% Visualization
fprintf('Generating figures...\n');
plot_within(A_mean, Nl_mean, offmask, cond_names, normalize_mode, split_tag, ...
    fig_root, sprintf('angle_within_%s_%s', normalize_mode, split_tag));

plot_across(obs_corr_r, null_corr_r, A_mean, offmask, pair_list, pair_names, cond_names, ...
    normalize_mode, split_tag, fig_root, sprintf('angle_across_%s_%s', normalize_mode, split_tag));

%% Print summary
fprintf('\n=== %s %s [normalize=%s, split=%s] ===\n', monkey, vp, normalize_mode, split_tag);
fprintf('\n(a) Within-rendering: mean |cos| and mean angle (deg), observed vs random-orientation null\n');
fprintf('%-5s | %-22s | %-22s\n', 'Cond', 'observed |cos| / ang', 'null |cos| / ang');
for i = 1:nconds
    cn = cond_names{i};
    fprintf('%-5s | %.3f / %5.1f deg        | %.3f / %5.1f deg\n', ...
        upper(cn), within_summary.(cn).obs_mean_cos, within_summary.(cn).obs_mean_ang, ...
        within_summary.(cn).null_mean_cos, within_summary.(cn).null_mean_ang);
end

fprintf('\n(b) Across-rendering: Pearson r, observed mean [95%% CI] vs shuffle null\n');
for p = 1:n_pairs
    pn = pair_names{p};
    fprintf('%-6s | obs %+.3f [%+.3f, %+.3f]  | null %+.3f [%+.3f, %+.3f]  (rho obs %+.3f)\n', ...
        pn, across_summary.(pn).obs_r_mean, across_summary.(pn).obs_r_ci(1), across_summary.(pn).obs_r_ci(2), ...
        across_summary.(pn).null_r_mean, across_summary.(pn).null_r_ci(1), across_summary.(pn).null_r_ci(2), ...
        across_summary.(pn).obs_rho_mean);
end


%% ====================== Local functions ======================

function [pc1, m] = cloud_pc1_mean(Xs, use_split)
% Return PC1 direction (1xN unit) and mean (1xN) of a cloud.
% use_split: PC1 from one disjoint half, mean from the other (bias-free).
% else: PC1 and mean both from a bootstrap resample of all trials.
n = size(Xs, 1);
if use_split
    p  = randperm(n);
    h  = floor(n/2);
    hA = p(1:h);
    hB = p(h+1:end);
    pc1 = top_pc(Xs(hA,:));
    m   = mean(Xs(hB,:), 1);
else
    bidx = randi(n, n, 1);
    Xb   = Xs(bidx,:);
    pc1  = top_pc(Xb);
    m    = mean(Xb, 1);
end
end

function v = top_pc(X)
% Top principal axis (1xN unit row) of the centered cloud X.
Xc = X - mean(X, 1);
[~, ~, V] = svd(Xc, 'econ');
v = V(:,1)';
end

function A = alignment_matrix(PC1, M, offmask)
% A(i,j) = |cos| between PC1_i and (M_j - M_i). Diagonal NaN.
nstim = size(M, 1);
A = zeros(nstim);
for i = 1:nstim
    D = M - M(i,:);
    nrm = sqrt(sum(D.^2, 2));
    nrm(i) = 1;                       % avoid divide-by-zero at i==j
    Dn = D ./ nrm;
    A(i,:) = abs(Dn * PC1(i,:)');
end
A(~offmask) = NaN;                    % diagonal undefined
end

function plot_within(A_mean, Nl_mean, offmask, cond_names, normalize_mode, split_tag, fig_root, fname)
% (a) Per rendering: distribution of per-pair angles (deg), observed vs
% random-orientation null. 90 deg = orthogonal reference.
nconds = numel(A_mean);
fig = figure('Position', [100 100 1500 430], 'Color', 'w');
col_obs  = [0.20 0.40 0.80];
col_null = [0.60 0.60 0.60];
edges = 0:5:90;

for i = 1:nconds
    subplot(1, nconds, i); hold on;
    a_obs  = A_mean{i}(offmask);
    a_null = Nl_mean{i}(offmask);
    ang_obs  = acosd(min(max(a_obs, 0), 1));
    ang_null = acosd(min(max(a_null, 0), 1));

    histogram(ang_null, edges, 'Normalization', 'probability', ...
        'FaceColor', col_null, 'EdgeColor', 'none', 'FaceAlpha', 0.45);
    histogram(ang_obs, edges, 'Normalization', 'probability', ...
        'FaceColor', col_obs, 'EdgeColor', 'none', 'FaceAlpha', 0.55);

    xline(90, 'k:', 'LineWidth', 1.0);
    xline(median(ang_obs),  '-',  'Color', col_obs,  'LineWidth', 1.4);
    xline(median(ang_null), '--', 'Color', col_null, 'LineWidth', 1.4);

    xlim([0 92]); xlabel('angle PC1 \angle signal (deg)'); ylabel('proportion');
    title(sprintf('%s  (med %.1f vs %.1f)', upper(cond_names{i}), ...
        median(ang_obs), median(ang_null)));
    set(gca, 'Box', 'on');
    if i == 1
        legend({'null (random orient.)', 'observed'}, 'Location', 'northwest', 'Box', 'off');
    end
end
sgtitle(sprintf('(a) Within-rendering PC1\\angle signal angles  [%s, %s]', normalize_mode, split_tag));
% saveas(fig, fullfile(fig_root, [fname '.fig']));
saveas(fig, fullfile(fig_root, [fname '.png']));
end

function plot_across(obs_corr_r, null_corr_r, A_mean, offmask, pair_list, pair_names, cond_names, normalize_mode, split_tag, fig_root, fname)
% (b) Per rendering pair: observed cross-rendering correlation distribution
% vs shuffle null (violins), with the mean-matrix point estimate marked.
n_pairs = numel(pair_list);
fig = figure('Position', [100 100 1400 470], 'Color', 'w');
col_obs  = [0.20 0.40 0.80];
col_null = [0.60 0.60 0.60];

% shared y-limits
allv = [obs_corr_r(:); null_corr_r(:)];
ylims = [min(allv) - 0.05, max(allv) + 0.05];

for p = 1:n_pairs
    i = pair_list{p}(1); j = pair_list{p}(2);
    subplot(1, n_pairs, p); hold on;

    violin_at(reshape(null_corr_r(:,:,p), [], 1), 1, 0.6, col_null);
    violin_at(obs_corr_r(:,p),                    2, 0.6, col_obs);

    % mean-matrix point estimate (less attenuated than the per-iter mean)
    pe = corr(A_mean{i}(offmask), A_mean{j}(offmask), 'type', 'Pearson');
    plot(2, pe, 'k^', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

    yline(0, 'k:', 'LineWidth', 0.5);
    xlim([0.4 2.6]); ylim(ylims);
    xticks([1 2]); xticklabels({'null', 'observed'});
    ylabel('cross-rendering Pearson r');
    title(sprintf('%s vs %s', upper(cond_names{i}), upper(cond_names{j})));
    set(gca, 'Box', 'on');
end
sgtitle(sprintf('(b) Cross-rendering preservation of PC1\\angle signal alignment  [%s, %s]', normalize_mode, split_tag));
% saveas(fig, fullfile(fig_root, [fname '.fig']));
saveas(fig, fullfile(fig_root, [fname '.png']));
end

function violin_at(data, x_center, width, color)
% Manual violin (kernel density), centered at x_center, with median bar.
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
