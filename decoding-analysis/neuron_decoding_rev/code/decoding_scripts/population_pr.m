function pr = population_pr(X, labels, opts)
%{
Filename: population_pr.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Participation ratio (PR) of one population's trial-averaged signal manifold,
    in two versions:

      raw : the definition used throughout the rebuttal
            (procrustes_decoding_basic_incremental_pr.m): trial-average the
            z-scored trials to [n_stim x N], center across stimuli, and take
            PR = (sum lambda)^2 / sum(lambda^2) over the covariance spectrum.
            Trial noise in the 10-trial means adds near-isotropic variance, so
            this PR rises when the signal-to-noise ratio falls.
      cv  : noise-corrected PR by cross-validated PCA (Stringer et al., 2019).
            The repeats of each stimulus are split into two halves; the variance
            of half A's means along half A's principal directions is estimated
            from the covariance with half B's means, which carry independent
            noise. Averaged over both roles and opts.n_split random splits,
            negative estimates set to 0.

    Also returns the cumulative variance fraction of the raw spectrum, so the
    variance captured by the top k PCs is pr.cumvar(k).

Inputs:
    X       : [n_rows x N] trials, normally z-scored per neuron (as for decoding)
    labels  : [n_rows x 1] stimulus labels, blocked 1..n_stim
    opts    : optional struct
                .n_split  random half-splits for the cv estimate (default 20)
                .seed     seed of those splits (default 1)

Output: struct pr
    raw, cv          participation ratios
    cumvar           [1 x r] cumulative variance fraction of the raw spectrum
    cv_cumvar        [1 x r] the same for the noise-corrected spectrum
    n_neurons
%}

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'n_split') || isempty(opts.n_split), opts.n_split = 20; end
if ~isfield(opts, 'seed')    || isempty(opts.seed),    opts.seed = 1; end

labels = labels(:);
n_stim = max(labels);
n_rep  = numel(labels) / n_stim;
N      = size(X, 2);

%% raw PR (rebuttal definition)
Mc  = class_means(X, labels, n_stim);
Mc  = Mc - mean(Mc, 1);
lam = svd(Mc, 'econ').^2;
lam = lam(lam > 0);
pr.raw    = sum(lam)^2 / sum(lam.^2);
pr.cumvar = (cumsum(lam) / sum(lam))';

%% noise-corrected PR (cross-validated PCA)
r   = min(n_stim - 1, N);
h   = floor(n_rep / 2);
acc = zeros(1, r);
rs  = RandStream('twister', 'Seed', opts.seed);
for t = 1:opts.n_split
    inA = false(numel(labels), 1);
    for s = 1:n_stim
        rows = find(labels == s);
        inA(rows(randperm(rs, n_rep, h))) = true;
    end
    A = class_means(X(inA, :),  labels(inA),  n_stim);  A = A - mean(A, 1);
    B = class_means(X(~inA, :), labels(~inA), n_stim);  B = B - mean(B, 1);
    [~, ~, VA] = svd(A, 'econ');  VA = VA(:, 1:r);
    [~, ~, VB] = svd(B, 'econ');  VB = VB(:, 1:r);
    lamA = sum((A * VA) .* (B * VA), 1);   % signal variance along A's directions
    lamB = sum((B * VB) .* (A * VB), 1);   % and along B's
    acc  = acc + (lamA + lamB) / 2;
end
lam_cv = acc / opts.n_split;
lam_cv(lam_cv < 0) = 0;
pr.cv        = sum(lam_cv)^2 / sum(lam_cv.^2);
pr.cv_cumvar = cumsum(lam_cv) / sum(lam_cv);
pr.n_neurons = N;
end

function M = class_means(X, y, n_stim)
M = zeros(n_stim, size(X, 2));
for s = 1:n_stim
    M(s, :) = mean(X(y == s, :), 1);
end
end
