%{
% FILENAME: procrustes_decoding_simulation.m
% AUTHOR:   Zitong Wang
% DATE:     2026-07-01
%
% DESCRIPTION:
%   Shared driver for the Procrustes transfer-decoding benchmark on synthetic
%   data. One run = R_pop independent synthetic populations of a chosen model.
%   For each population it generates trial data, runs Procrustes transfer
%   decoding for every ordered cue pair, and records self-decoding, PT, the
%   shuffle null, and the model's signal participation ratio. The spread across
%   populations is the error bar (population + trial + decode variability
%   combined -- one population = one complete synthetic experiment).
%
%   The model is selected by a function handle (cfg.model_fn) with the contract
%       [data_trial, labels, diagnostics] = model_fn(gen_cfg)
%   so Models 1/2/3/5/6 share this driver unchanged; only the generator differs.
%   n_cue is inferred from the generator output, so 2-cue (Model 3) and 3-cue
%   (Model 5) runs both work.
%
%   Decoding math is identical to procrustes_decoding_basic.m (z-score per cue,
%   MATLAB `procrustes` similarity transform, ECOC SVM). Two deliberate changes:
%     - No neuron sampling: all N synthetic neurons are used (N is set directly,
%       so the real-data neuron-subsampling loop has no analogue here).
%     - The target ECOC model and the self-decoding cross-validation are fit
%       ONCE per (population, cue-pair). They do not depend on the stim1 trial
%       holdout, so only the transfer estimate is averaged over holdouts.
%
% OUTPUT:
%   results/decoding_outputs/<model_name>/<variant>/simulation_results.mat
%     containing struct `results` with cfg, seeds, per-population metrics,
%     and cross-population aggregates. Datasets are NOT saved -- any population
%     is exactly reproducible from cfg + its seed.
%
%   Per-population metric columns (matching procrustes_decoding_basic.m):
%     1 self (genAcc)   2 no-transform   3 PT (transformed)   4 rand null
%     5 scale-only      6 rotation-only  7 translation-only   8 chance control
%}

%% INITIALIZATION
clc; clear;

%% CONFIGURATION
addpath(fullfile('..','generation_scripts'))
cfg = struct();

% --- model selection ---
cfg.model_name   = 'model1';                      % tags the output folder
cfg.model_fn     = @generate_model1_trial_data;   % generator handle
cfg.noise        = 'poisson';                     % label only (tags variant)

% --- generator parameters (passed through to model_fn) ---
cfg.N            = 50;        % neurons (all used; open parameter for later sweep)
cfg.weight_scale = 10;        % SD of readout weights (SNR knob)
cfg.rho          = 0;          % cross-cue weight correlation (0 = independent)
cfg.baseline     = "auto";     % positivity margin that tracks weight_scale

% --- run parameters ---
cfg.R_pop        = 20;        % number of populations (outer loop = error bar)
cfg.base_seed    = 1000;       % population seeds = base_seed + (1:R_pop)
cfg.trial_repeat = 10;         % trial-holdout resamples per (population, pair)
cfg.n_folds      = 10;         % k-fold for self-decoding (5 is a fine speedup)

%% OUTPUT LOCATION
% Layered, variant-tagged so runs never overwrite (rho / N / noise vary).
variant   = sprintf('rho%g_N%d_ws%g_%s', cfg.rho, cfg.N, cfg.weight_scale, cfg.noise);
save_path = fullfile('..','..','results','decoding_outputs', cfg.model_name, variant);
if ~exist(save_path, 'dir'), mkdir(save_path); end

%% RUN
tic;
results         = run_simulation(cfg);
results.variant = variant;
toc;

%% SAVE
save(fullfile(save_path, 'simulation_results.mat'), 'results');
fprintf('Saved %s [%s], R_pop=%d, to %s\n', ...
    cfg.model_name, variant, cfg.R_pop, save_path);

%% SUMMARY
print_summary(results);


%% ======================================================================
%% Local functions
%% ======================================================================
function results = run_simulation(cfg)
% Outer loop over independent populations (parfor). Each iteration is one
% complete synthetic experiment, reproducible from its seed alone.
R        = cfg.R_pop;
seeds    = cfg.base_seed + (1:R).';
model_fn = cfg.model_fn;                 % local handle for parfor broadcast
pop      = cell(R, 1);                   % sliced cell output (parfor-safe)

parfor r = 1:R
    gen_cfg = build_gen_cfg(cfg, seeds(r));
    [data_trial, labels, gdiag] = model_fn(gen_cfg);
    pop{r} = decode_population(data_trial, labels, gdiag, cfg);
end

results = aggregate_results(pop, cfg, seeds);
end

% ----------------------------------------------------------------------
function out = decode_population(data_trial, labels, gdiag, cfg)
% Procrustes transfer decoding for every ordered cue pair of one population.
n_cue = numel(data_trial);

% z-score each cue once (per-neuron), reused across all pairs.
Z = cell(1, n_cue);
for c = 1:n_cue
    Z{c} = zscore(data_trial{c});
end

pairs   = ordered_pairs(n_cue);          % [n_pairs x 2]: (stim1=transformed, stim2=target)
n_pairs = size(pairs, 1);
names   = {'ac','ec','ex'};
pair_metrics = zeros(n_pairs, 8);
pair_labels  = cell(n_pairs, 1);

for p = 1:n_pairs
    s1 = pairs(p,1);  s2 = pairs(p,2);
    stim1_data = Z{s1};                  % comparison (to be transformed)
    stim2_data = Z{s2};                  % target (decoder trained here)
    pair_labels{p} = sprintf('%s->%s', names{s1}, names{s2});

    % --- target model + self-decoding: fit ONCE (holdout-independent) ---
    stim2_model = fitcecoc(stim2_data, labels);
    genAcc      = 1 - kfoldLoss(crossval(stim2_model, 'KFold', cfg.n_folds));
    stim2_mean  = take_average(stim2_data, 10);        % [50 x N] target means

    % --- transfer: average over stim1 trial holdouts ---
    T   = cfg.trial_repeat;
    acc = zeros(T, 7);   % [no-tf, PT, rand, scale, rot, trans, chance]
    for t = 1:T
        valid_idx = sample_trial(1);                    % 1 test trial/stim (50 rows)
        train_idx = setdiff((1:500).', valid_idx);      % 9 trials/stim (450 rows)

        stim1_test = stim1_data(valid_idx, :);
        test_label = labels(valid_idx);
        stim1_mean = take_average(stim1_data(train_idx, :), 9);

        % similarity transform mapping stim1 means -> stim2 means
        [~, ~, tf] = procrustes(stim2_mean, stim1_mean);
        transformed = tf.b * stim1_test * tf.T + tf.c;

        % shuffled-correspondence transform (PT null)
        sh = randperm(50);
        [~, ~, tf0] = procrustes(stim2_mean, stim1_mean(sh, :));
        toy = tf0.b * stim1_test * tf0.T + tf0.c;

        % single-component partial transforms
        scale_only = tf.b * stim1_test;
        rot_only   = stim1_test * tf.T;
        trans_only = stim1_test + tf.c;

        pred_s1 = predict(stim2_model, stim1_test);
        acc(t,1) = mean(pred_s1                              == test_label);
        acc(t,2) = mean(predict(stim2_model, transformed)   == test_label);
        acc(t,3) = mean(predict(stim2_model, toy)           == test_label);
        acc(t,4) = mean(predict(stim2_model, scale_only)    == test_label);
        acc(t,5) = mean(predict(stim2_model, rot_only)      == test_label);
        acc(t,6) = mean(predict(stim2_model, trans_only)    == test_label);
        acc(t,7) = mean(randperm(50).'                      == pred_s1);   % chance control
    end
    acc_mean = mean(acc, 1);
    pair_metrics(p,:) = [genAcc, acc_mean];              % 8 columns
end

out = struct();
out.pair_metrics = pair_metrics;         % [n_pairs x 8]
out.pair_labels  = pair_labels;          % {n_pairs x 1}
out.signal_PR    = gdiag.signal_PR;      % [1 x n_cue]
out.mean_count   = gdiag.mean_count;     % scalar
out.seed         = gdiag.cfg.seed;       % scalar
end

% ----------------------------------------------------------------------
function results = aggregate_results(pop, cfg, seeds)
% Stack per-population results and compute cross-population statistics.
R        = numel(pop);
n_pairs  = size(pop{1}.pair_metrics, 1);
n_metric = size(pop{1}.pair_metrics, 2);
n_cue    = numel(pop{1}.signal_PR);

PM = zeros(R, n_pairs, n_metric);        % [R x n_pairs x 8]
PR = zeros(R, n_cue);
MC = zeros(R, 1);
for r = 1:R
    PM(r,:,:) = pop{r}.pair_metrics;
    PR(r,:)   = pop{r}.signal_PR;
    MC(r)     = pop{r}.mean_count;
end

per_pop = struct();
per_pop.pair_metrics = PM;
per_pop.signal_PR    = PR;
per_pop.mean_count   = MC;

agg = struct();
agg.mean    = reshape(mean(PM, 1),        [n_pairs, n_metric]);
agg.sem     = reshape(std(PM, 0, 1),      [n_pairs, n_metric]) / sqrt(R);
agg.p2p5    = reshape(prctile(PM, 2.5, 1),[n_pairs, n_metric]);
agg.p97p5   = reshape(prctile(PM,97.5, 1),[n_pairs, n_metric]);
agg.PR_mean = mean(PR, 1);
agg.PR_sem  = std(PR, 0, 1) / sqrt(R);

results = struct();
results.cfg          = cfg;
results.seeds        = seeds;
results.metric_names = {'self','no_transform','PT','rand_null', ...
                        'scale_only','rotation_only','translation_only','chance'};
results.pair_labels  = pop{1}.pair_labels;
results.per_pop      = per_pop;
results.agg          = agg;
end

% ----------------------------------------------------------------------
function g = build_gen_cfg(cfg, seed)
% Pass only generator-relevant fields; each generator fills its own defaults.
g = struct();
g.N            = cfg.N;
g.weight_scale = cfg.weight_scale;
g.rho          = cfg.rho;
g.baseline     = cfg.baseline;
if isfield(cfg, 'baseline_k'), g.baseline_k = cfg.baseline_k; end
g.seed         = seed;
end

% ----------------------------------------------------------------------
function pairs = ordered_pairs(n)
% All ordered cue pairs (i,j), i~=j. n=2 -> [1 2; 2 1].
pairs = zeros(n*(n-1), 2);
k = 0;
for i = 1:n
    for j = 1:n
        if i ~= j
            k = k + 1;
            pairs(k,:) = [i, j];
        end
    end
end
end

% ----------------------------------------------------------------------
function stim_data_trial_averged = take_average(stim_data, number_of_average)
% Mean over each stimulus's trial block. Copied from procrustes_decoding_basic.m.
[~, neuron] = size(stim_data);
stim_data_trial_averged = zeros(50, neuron);
for i = 1:50
    temp = stim_data(i*number_of_average-(number_of_average-1):i*number_of_average, :);
    stim_data_trial_averged(i,:) = mean(temp, 1);
end
end

% ----------------------------------------------------------------------
function sample_trial_idx = sample_trial(trialnum_outof_ten)
% Random holdout of trialnum_outof_ten trials per stimulus (blocked layout).
% Copied from procrustes_decoding_basic.m.
sample_trial_matrix = zeros(50, 10);
for i = 1:50
    sample_trial_matrix(i,:) = randperm(10);
end
sample_trial_4condition = sample_trial_matrix(:, 1:trialnum_outof_ten);
base = (0:10:490).';
sample_trial_idx = sample_trial_4condition + base;
sample_trial_idx = reshape(sample_trial_idx, [], 1);
sample_trial_idx = sort(sample_trial_idx);
end

% ----------------------------------------------------------------------
function print_summary(results)
% Console summary: signal PR and the key metrics per cue pair.
agg    = results.agg;
labels = results.pair_labels;
R      = size(results.per_pop.pair_metrics, 1);
fprintf('\n=== %s [%s]: mean over %d populations ===\n', ...
    results.cfg.model_name, results.variant, R);
fprintf('signal PR per cue: %s\n', num2str(agg.PR_mean, '%.3f '));
for p = 1:numel(labels)
    fprintf('%-10s  self=%.3f  PT=%.3f rot-PT=%.3f shuffle-pt=%.3f no-tf=%.3f  null=%.3f\n', ...
        labels{p}, agg.mean(p,1), agg.mean(p,3), agg.mean(p,6), agg.mean(p,4), agg.mean(p,2), agg.mean(p,4));
end
fprintf('(chance = %.3f)\n', 1/50);
end
