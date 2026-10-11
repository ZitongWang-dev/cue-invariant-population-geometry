function out = pt_transfer_core(Xs, Xt, labels, opts)
%{
Filename: pt_transfer_core.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Procrustes transfer decoding of ONE source -> target comparison, with the
    corrected (self-consistent) correspondence-shuffled control. Every driver in
    neuron_decoding_rev calls this function, so the control, the trial split and
    the PC projection are identical for cue transfer, within-area halves and
    across-area decoding.

    Source = the population whose held-out trials are transformed and tested
    (the "comparison" / stim1 of the old scripts). Target = the population the
    decoder is trained on (the "target" / stim2). Both matrices share one row
    order (stimulus blocks of n_rep repeats).

    Trial split: a 10-fold PARTITION of rows. Within each stimulus the n_rep row
    indices are shuffled once and the k-th of them goes to fold k, so every trial
    is held out exactly once and each fold holds one trial per stimulus. The same
    partition is used on both sides and always moves whole rows, so the
    within-session noise correlations of the pseudo-population are kept.

    For fold f (training rows = all rows outside fold f):
      1. Target: class means, (PC mode) PCA basis, and an ECOC linear SVM
         (fitcecoc defaults, as the published pipeline: one-vs-one linear SVMs;
         support vectors are not stored, which does not change linear-SVM
         predictions), all from training rows.
      2. Source: class means from training rows; (PC mode) its own PCA basis.
      3. Procrustes (MATLAB, similarity fit) from the source means onto the
         target means. PT_rot applies the rotation T only (the headline, as in
         the paper); PT_full applies b*X*T + c.
      4. Self-decoding = the fold model on the TARGET's fold-f trials, so the
         10 fold models give exactly a 10-fold CV, and the target model, means
         and basis never contain the trials they are tested on. For halves (same
         rows on both sides) this also means the target half never trained on
         the physical trial the source half is tested on.
      5. Transfer = the same model on the SOURCE's fold-f trials:
           nopt     no transform (full space with equal neuron counts; in PC
                    mode only when opts.same_neurons, projected onto the
                    TARGET basis). For different neurons it is a sanity check
                    that should sit at chance.
           pt_rot   test * T
           pt_full  b * test * T + c
           ctrl     CORRECTED CONTROL, n_perm draws per fold: one permutation p
                    of the stimuli is applied to BOTH the source training means
                    and the source test trials, a rotation-only fit from the
                    permuted means onto the target means is applied to the
                    permuted test trials, and predictions are scored against the
                    TARGET labels. This is what a rotation with the same degrees
                    of freedom achieves once the true correspondence is gone.
           ctrl_id  the same predictions scored against the permuted identities
                    (must sit at ~1/n_stim; checks that the shuffle really
                    destroyed the correspondence).
      6. In-sample residual shape distances of the true fit on the means:
         d_rot (rotation only, b = 1) and d_full (MATLAB's d), both with the
         target as denominator.

    PC mode (opts.k given): each side is projected onto the top-k principal
    directions of ITS OWN centered class-mean matrix (training rows only), as in
    neuronal_decoding/.../procrustes_decoding_basic_incremental_pc.m, so the
    rotation has k(k-1)/2 free parameters instead of N(N-1)/2.

Inputs:
    Xs, Xt  : [n_rows x Ns], [n_rows x Nt] spike counts (or z-scored data),
              rows in stimulus blocks; Ns == Nt is required in full space
    labels  : [n_rows x 1] stimulus labels, blocked 1..n_stim
    opts    : struct (all optional)
                .k             []: full neuron space (default); k: PC mode
                .n_perm        control permutations per fold (default 10)
                .seed          seed of the fold partition and permutations (default 1)
                .same_neurons  true when Xs and Xt are the same neurons (cue
                               transfer); enables non-PT in PC mode (default false)
                .zscore        z-score each side per neuron over all rows (default
                               true, as the published pipeline)

Output: struct out
    self, nopt, pt_rot, pt_full, d_rot, d_full   [n_folds x 1]
    ctrl, ctrl_id                                [n_folds x n_perm]
    fold      [n_rows x 1] fold of every row
    k         [] or the k used;  n_folds, n_perm, seed
%}

%% Options
if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'k'),                                         opts.k = []; end
if ~isfield(opts, 'n_perm')       || isempty(opts.n_perm),       opts.n_perm = 10; end
if ~isfield(opts, 'seed')         || isempty(opts.seed),         opts.seed = 1; end
if ~isfield(opts, 'same_neurons') || isempty(opts.same_neurons), opts.same_neurons = false; end
if ~isfield(opts, 'zscore')       || isempty(opts.zscore),       opts.zscore = true; end

labels = labels(:);
n_rows = numel(labels);
n_stim = max(labels);
n_rep  = n_rows / n_stim;
assert(size(Xs, 1) == n_rows && size(Xt, 1) == n_rows, 'pt_transfer_core: Xs, Xt and labels need the same rows.');
assert(n_rep == round(n_rep) && isequal(labels, kron((1:n_stim)', ones(n_rep, 1))), ...
    'pt_transfer_core: labels must come in equal, consecutive stimulus blocks.');

Ns = size(Xs, 2);  Nt = size(Xt, 2);
use_pc = ~isempty(opts.k);
if use_pc
    k = opts.k;
    assert(k >= 1 && k == round(k) && k <= min([n_stim - 1, Ns, Nt]), ...
        'pt_transfer_core: k = %g is outside 1..min(n_stim-1, Ns, Nt) = %d.', k, min([n_stim - 1, Ns, Nt]));
else
    k = [];
    assert(Ns == Nt, 'pt_transfer_core: full-space transfer needs equal neuron counts (got %d and %d).', Ns, Nt);
end
if opts.same_neurons
    assert(Ns == Nt, 'pt_transfer_core: same_neurons = true but the neuron counts differ.');
end
do_nopt = (~use_pc && Ns == Nt) || (use_pc && opts.same_neurons);

if opts.zscore
    Xs = zscore(Xs);
    Xt = zscore(Xt);
end

%% Fold partition (whole rows, shared by source and target)
rs = RandStream('twister', 'Seed', opts.seed);
fold = zeros(n_rows, 1);
for s = 1:n_stim
    rows = find(labels == s);
    fold(rows(randperm(rs, n_rep))) = 1:n_rep;
end
F = n_rep;
M = opts.n_perm;

%% Decode fold by fold
svm_template = templateSVM('SaveSupportVectors', false);
self    = nan(F, 1);  nopt   = nan(F, 1);
pt_rot  = nan(F, 1);  pt_full = nan(F, 1);
d_rot   = nan(F, 1);  d_full  = nan(F, 1);
ctrl    = nan(F, M);  ctrl_id = nan(F, M);

for f = 1:F
    tr  = fold ~= f;
    te  = find(fold == f);
    [~, o] = sort(labels(te));
    te  = te(o);                  % one trial per stimulus, ordered 1..n_stim
    ytr = labels(tr);
    yte = labels(te);             % = (1:n_stim)'

    Mt = class_means(Xt(tr, :), ytr, n_stim);
    Ms = class_means(Xs(tr, :), ytr, n_stim);

    if use_pc
        [mut, Vt] = mean_basis(Mt, k);
        [mus, Vs] = mean_basis(Ms, k);
        Ttr = (Xt(tr, :) - mut) * Vt;     % target training trials, target basis
        Tte = (Xt(te, :) - mut) * Vt;     % target held-out trials (self-decoding)
        Ste = (Xs(te, :) - mus) * Vs;     % source held-out trials, source basis
        Mtk = (Mt - mut) * Vt;
        Msk = (Ms - mus) * Vs;
        if do_nopt, Nte = (Xs(te, :) - mut) * Vt; end   % source trials in the target basis
    else
        Ttr = Xt(tr, :);  Tte = Xt(te, :);  Ste = Xs(te, :);
        Mtk = Mt;         Msk = Ms;
        if do_nopt, Nte = Ste; end
    end

    % ECOC linear SVM (fitcecoc defaults; support vectors are not stored, which
    % leaves linear-SVM predictions unchanged and saves ~30% of the fit time)
    mdl = fitcecoc(Ttr, ytr, 'Learners', svm_template);

    % true correspondence
    [d_full(f), ~, tf] = procrustes(Mtk, Msk);
    A0 = Mtk - mean(Mtk, 1);
    B0 = Msk - mean(Msk, 1);
    d_rot(f) = sum((A0 - B0 * tf.T).^2, 'all') / sum(A0.^2, 'all');

    % corrected control: the same permutation on the source means AND the source test trials
    perms = zeros(M, n_stim);
    ctrl_sets = cell(1, M);
    for m = 1:M
        p = randperm(rs, n_stim);
        perms(m, :) = p;
        [~, ~, tp] = procrustes(Mtk, Msk(p, :));
        ctrl_sets{m} = Ste(p, :) * tp.T;
    end

    % one predict call for every test set of this fold (each row is predicted
    % independently, so batching changes nothing but the run time)
    sets = [{Tte}, {Ste * tf.T}, {tf.b * Ste * tf.T + tf.c(1, :)}, ctrl_sets];
    if do_nopt, sets = [sets, {Nte}]; end
    pred = reshape(predict(mdl, vertcat(sets{:})), n_stim, []);

    self(f)    = mean(pred(:, 1) == yte);
    pt_rot(f)  = mean(pred(:, 2) == yte);
    pt_full(f) = mean(pred(:, 3) == yte);
    for m = 1:M
        ctrl(f, m)    = mean(pred(:, 3 + m) == yte);
        ctrl_id(f, m) = mean(pred(:, 3 + m) == perms(m, :)');
    end
    if do_nopt, nopt(f) = mean(pred(:, end) == yte); end
end

out = struct('self', self, 'nopt', nopt, 'pt_rot', pt_rot, 'pt_full', pt_full, ...
    'ctrl', ctrl, 'ctrl_id', ctrl_id, 'd_rot', d_rot, 'd_full', d_full, ...
    'fold', int8(fold), 'k', k, 'n_folds', F, 'n_perm', M, 'seed', opts.seed);
end

%% ---------------------------------------------------------------------------
function M = class_means(X, y, n_stim)
% [n_stim x N] mean response per stimulus
M = zeros(n_stim, size(X, 2));
for s = 1:n_stim
    M(s, :) = mean(X(y == s, :), 1);
end
end

function [mu, V] = mean_basis(M, k)
% Top-k principal directions of the centered class-mean matrix (same
% construction as procrustes_decoding_basic_incremental_pc.m).
mu = mean(M, 1);
[~, ~, V] = svd(M - mu, 'econ');
V = V(:, 1:k);
end
