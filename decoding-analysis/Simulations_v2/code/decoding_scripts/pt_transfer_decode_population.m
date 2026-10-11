%{
% FILENAME: pt_transfer_decode_population.m   (Simulations_v2)
% AUTHOR:   Zitong Wang
% DATE:     2026-09-17
%
% DESCRIPTION:
%   Procrustes transfer decoding of ONE synthetic population (all N neurons,
%   every ordered cue pair) with the self-consistent correspondence-shuffled
%   control. This is the per-population core called by
%   procrustes_decoding_simulation_sweep.m; it is a separate function so it can
%   be tested on its own and reused unchanged by future drivers (Model 5/6,
%   PC sweeps, noise-correlation injection).
%
%   Decoding math follows the real-data pipeline
%   (neuronal_decoding/code/decoding_scripts/procrustes_decoding_basic_incremental_pr.m
%   for the full space; ..._incremental_pc.m for the PC-constrained case):
%     - z-score each cue per neuron over all of its trials;
%     - stim2 (target): ECOC linear SVM on all stim2 trials; self-decoding =
%       1 - k-fold CV loss; class means from all n_trials trials;
%     - stim1 (transformed): one held-out trial per stimulus; class means from
%       the remaining n_trials-1 training trials;
%     - MATLAB `procrustes` (similarity: rotation T, scale b, translation c)
%       fit from the stim1 training means onto the stim2 means; the held-out
%       stim1 trials are mapped with the full transform and with each single
%       component, and decoded by the stim2 model against their true labels;
%     - CONTROL (the corrected one): one permutation pi of the stimuli is
%       applied to BOTH the stim1 training-mean matrix and the held-out stim1
%       test trials; a rotation-only map T_pi is fit from the permuted means
%       onto the stim2 means and applied to the permuted test trials; accuracy
%       is scored against the TARGET labels 1:n_stim. This is the accuracy a
%       free rotation can manufacture once the true correspondence has been
%       destroyed, i.e. the overfitting floor. Its matched true-correspondence
%       quantity is the rotation-only PT (column 6). The same predictions
%       scored against the permuted identity labels pi (column 9) must sit at
%       chance, which verifies that the correspondence really was destroyed.
%       (The v1 driver's "rand null" fit the map from shuffled means but
%       applied it to UNshuffled trials, which lands at chance regardless of
%       overfitting; that column is gone.)
%     - optional PC constraint (opts.pc_dims = k): each cue is projected into
%       the top-k of its own mean-PCA basis (stim2's from all trials, stim1's
%       from training trials only, leakage-clean) before the SVM and the
%       Procrustes fit, so the rotation has k(k-1)/2 free parameters instead
%       of N(N-1)/2.
%   As in the v1 driver, the stim2 model and its self-decoding are fit ONCE
%   per pair (they do not depend on the stim1 holdout); only the transfer
%   estimates are averaged over opts.trial_repeat holdouts.
%
% INPUTS:
%   data_trial : 1 x n_cue cell of [n_stim*n_trials x N] counts, blocked by stimulus
%   labels     : [n_stim*n_trials x 1] stimulus labels for that row order
%   opts       : struct (all optional)
%                  trial_repeat  holdout resamples per pair          (default 10)
%                  n_folds       k-fold for self-decoding             (default 10)
%                  pc_dims       [] = full space, or k = retained PCs (default [])
%                  cue_names     cue labels                           (default {'ac','ec','ex'})
%
% OUTPUT: struct `out` with
%   pair_labels  {n_pairs x 1}  e.g. 'ac->ec' (stim1 -> stim2)
%   pair_metrics [n_pairs x 9]  accuracies, columns:
%      1 self (genAcc)      2 no_transform      3 PT_full (b*X*T+c)
%      4 control            5 scale_only        6 PT_rotation (X*T)   <- matched to 4
%      7 translation_only   8 chance            9 control_vs_identity (~chance)
%   metric_names {1 x 9}
%   pair_pdist   [n_pairs x 2]  in-sample Procrustes residual d of the true fit
%                               and of the shuffled-correspondence fit
%   PR_emp       [1 x n_cue]    PR of the z-scored trial-averaged means (the
%                               real-data definition, trial noise included)
%   k_eff        retained PCs actually used ([] in the full space)
%}

function out = pt_transfer_decode_population(data_trial, labels, opts)

%% Options
if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'trial_repeat') || isempty(opts.trial_repeat), opts.trial_repeat = 10; end
if ~isfield(opts, 'n_folds')      || isempty(opts.n_folds),      opts.n_folds      = 10; end
if ~isfield(opts, 'pc_dims'),                                     opts.pc_dims      = []; end
if ~isfield(opts, 'cue_names')    || isempty(opts.cue_names),    opts.cue_names    = {'ac','ec','ex'}; end

n_cue    = numel(data_trial);
n_stim   = max(labels);
n_rows   = numel(labels);
n_trials = n_rows / n_stim;
T        = opts.trial_repeat;
use_pc   = ~isempty(opts.pc_dims);

%% z-score each cue once (per neuron, over all of its trials)
Z = cell(1, n_cue);
for c = 1:n_cue
    Z{c} = zscore(data_trial{c});
end

%% Empirical PR of each cue's z-scored trial-averaged means (real-data definition)
PR_emp = zeros(1, n_cue);
for c = 1:n_cue
    M = class_means(Z{c}, labels, n_stim);
    PR_emp(c) = participation_ratio(M - mean(M, 1));
end

%% Every ordered cue pair: stim1 = transformed (comparison), stim2 = target (trained)
pairs        = ordered_pairs(n_cue);
n_pairs      = size(pairs, 1);
pair_metrics = zeros(n_pairs, 9);
pair_pdist   = zeros(n_pairs, 2);
pair_labels  = cell(n_pairs, 1);
k_eff        = [];

for p = 1:n_pairs
    s1 = pairs(p, 1);  s2 = pairs(p, 2);
    X1 = Z{s1};  X2 = Z{s2};
    pair_labels{p} = sprintf('%s->%s', opts.cue_names{s1}, opts.cue_names{s2});

    % ---- target (stim2): means, optional PC basis, ECOC model, self-decoding: once per pair ----
    M2 = class_means(X2, labels, n_stim);
    if use_pc
        k_eff = min([opts.pc_dims, size(X2, 2), n_stim - 1]);
        [mu2, V2] = mean_basis(M2, k_eff);
        X2 = (X2 - mu2) * V2;
        M2 = (M2 - mu2) * V2;
    end
    model2 = fitcecoc(X2, labels);
    genAcc = 1 - kfoldLoss(crossval(model2, 'KFold', opts.n_folds));

    % ---- transfer: average over stim1 trial holdouts ----
    acc = zeros(T, 8);
    d   = zeros(T, 2);
    for t = 1:T
        % one held-out trial per stimulus; rows ordered by stimulus, so labels are 1:n_stim
        test_idx   = (0:n_stim-1).' * n_trials + randi(n_trials, n_stim, 1);
        train_idx  = setdiff((1:n_rows).', test_idx);
        X1_test    = X1(test_idx, :);
        test_label = labels(test_idx);
        M1 = class_means(X1(train_idx, :), labels(train_idx), n_stim);
        if use_pc
            [mu1, V1] = mean_basis(M1, k_eff);          % stim1 basis from TRAINING trials only
            X1_test = (X1_test - mu1) * V1;
            M1      = (M1 - mu1) * V1;
        end

        % true-correspondence similarity transform: stim1 training means -> stim2 means
        [d_true, ~, tf] = procrustes(M2, M1);
        full_tf = tf.b * X1_test * tf.T + tf.c;
        scale_o = tf.b * X1_test;
        rot_o   = X1_test * tf.T;
        trans_o = X1_test + tf.c;

        % corrected control: the SAME permutation on the means AND the test trials, rotation only
        pi_ = randperm(n_stim);
        [d_sh, ~, tf0] = procrustes(M2, M1(pi_, :));
        ctrl           = X1_test(pi_, :) * tf0.T;
        identity_label = test_label(pi_);

        pred_notf = predict(model2, X1_test);
        pred_ctrl = predict(model2, ctrl);
        acc(t, 1) = mean(pred_notf                  == test_label);      % no transform
        acc(t, 2) = mean(predict(model2, full_tf)   == test_label);      % PT, full similarity transform
        acc(t, 3) = mean(pred_ctrl                  == test_label);      % control vs target labels
        acc(t, 4) = mean(predict(model2, scale_o)   == test_label);      % scale only
        acc(t, 5) = mean(predict(model2, rot_o)     == test_label);      % rotation only (headline)
        acc(t, 6) = mean(predict(model2, trans_o)   == test_label);      % translation only
        acc(t, 7) = mean(randperm(n_stim).'         == pred_notf);       % chance control
        acc(t, 8) = mean(pred_ctrl                  == identity_label);  % control vs identity labels
        d(t, :)   = [d_true, d_sh];
    end
    pair_metrics(p, :) = [genAcc, mean(acc, 1)];
    pair_pdist(p, :)   = mean(d, 1);
end

out = struct();
out.pair_labels  = pair_labels;
out.pair_metrics = pair_metrics;
out.metric_names = {'self', 'no_transform', 'PT_full', 'control', 'scale_only', ...
                    'PT_rotation', 'translation_only', 'chance', 'control_vs_identity'};
out.pair_pdist   = pair_pdist;
out.PR_emp       = PR_emp;
out.k_eff        = k_eff;

end  % main function


%% ----------------------------------------------------------------------
%% Local functions
%% ----------------------------------------------------------------------
function M = class_means(X, lab, n_stim)
% Per-stimulus mean of the rows of X. Works for any trial count or row order.
M = zeros(n_stim, size(X, 2));
for s = 1:n_stim
    M(s, :) = mean(X(lab == s, :), 1);
end
end

function [mu, V] = mean_basis(M, k)
% Top-k principal directions of the centered condition-mean manifold
% (separate per cue, as in procrustes_decoding_basic_incremental_pc.m).
mu = mean(M, 1);
[~, ~, V] = svd(M - mu, 'econ');
V = V(:, 1:k);
end

function pr = participation_ratio(Mc)
% PR = (sum lambda)^2 / sum(lambda^2) over the covariance spectrum of the
% centered [n_stim x N] mean matrix, via SVD (rank-safe when N > n_stim).
s  = svd(Mc, 'econ');
ev = s.^2;
ev = ev(ev > 0);
pr = sum(ev)^2 / sum(ev.^2);
end

function pairs = ordered_pairs(n)
% All ordered cue pairs (i, j), i ~= j. n = 2 -> [1 2; 2 1].
pairs = zeros(n * (n - 1), 2);
k = 0;
for i = 1:n
    for j = 1:n
        if i ~= j
            k = k + 1;
            pairs(k, :) = [i, j];
        end
    end
end
end
