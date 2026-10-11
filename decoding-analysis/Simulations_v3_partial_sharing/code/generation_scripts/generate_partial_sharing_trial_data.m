%{
% FILENAME: generate_partial_sharing_trial_data.m   (Simulations_v3_partial_sharing)
% AUTHOR:   Zitong Wang
% DATE:     2026-10-07
%
% DESCRIPTION:
%   Partially shared linear mixed selectivity: the reviewer's "linear
%   combination of contour shape parameters (the same across renderings) and
%   rendering style/appearance parameters (different)". Every cue has the same
%   number K of latent dimensions; K_s of them are SHARED across cues (shape)
%   and the remaining K - K_s are drawn independently per cue (appearance).
%   The amount of shared structure is
%       alpha = K_s / K,
%   so alpha = 0 is a null with no shared structure at the same dimensionality
%   as the recordings, and alpha = 1 is Simulations_v2 Model 3 (fully shared
%   K-dimensional latent, independent per-cue readout).
%
%   Fixed-K design: each cue's latent is whitened as a whole (identity
%   covariance over its K dimensions) whatever alpha is, so every within-cue
%   quantity (each neuron's tuning depth, the signal spectrum, the
%   participation ratio, and hence self-decoding) has the same distribution at
%   every alpha. Only the relation BETWEEN cues changes. A parametrization that
%   weights a shared and a cue-specific subspace by sqrt(alpha) and
%   sqrt(1-alpha) would instead change the per-cue PR with alpha
%   (PR = 1/(alpha^2/D_s + (1-alpha)^2/d_c)), and the correspondence-shuffled
%   floor depends on dimensionality, so sharing and floor would be confounded.
%
%   Generative model (neuron i, stimulus s, cue c, trial t):
%       r_i, cv_i      mean count and tuning depth, as in Simulations_v2 Models 1 and 3
%       x_s            [1 x K] latent draw of stimulus s (cue 1's latent)
%       z_{s,1}        = whiten(x_s)
%       z_{s,c}        = whiten([x_s(1:K_s), y_{s,c}]),  y_{s,c} ~ N(0, I_{K-K_s}) for c >= 2
%                        Whitening is Cholesky (Gram-Schmidt in column order), so the
%                        first K_s whitened columns are the SAME in every cue, and the
%                        last K - K_s are orthogonalized against them and independent
%                        across cues. Each cue's latent has zero mean, identity covariance.
%       u_{i,c}        unit readout direction in the K-dim latent space, random per
%                      neuron and cue (rho-mixed across cues, as in v2 Model 3)
%       lambda_{i,s,c} = r_i + r_i * cv_i * u_{i,c} . z_{s,c}
%       count          ~ Poisson(lambda)
%   Positivity: each neuron's depth is capped at the largest value that keeps
%   its own rates non-negative in every cue (cap_rule 'exact'), and at cv_max.
%
%   Pairing across alpha and with Simulations_v2. The RNG is consumed in the
%   order r, cv, x, readouts (exactly as v2 Model 3), and only then the
%   cue-specific columns y (none when K_s = K). Hence:
%     - n_shared = n_latent reproduces generate_model3_trial_data.m
%       (n_latent = 5, cap_rule 'exact') bit for bit, i.e. the
%       sweep_N_S_model3_D5 populations;
%     - at every alpha a seed gives the same neurons (r, cv), the same shared
%       latent and the same readouts, and cue 1's latent is identical at
%       every alpha. Differences between alphas at a seed come from the
%       cue-specific columns, the depth cap they imply and the Poisson draws.
%
%   rho: cross-cue readout correlation, as in v2 Model 3 (0 = independent
%   readouts, the default; 1 = identical readouts, i.e. each neuron keeps the
%   same filter in every cue).
%
% OUTPUT CONTRACT (same as the Simulations_v2 generators, used by
%   Simulations_v2/code/decoding_scripts/pt_transfer_decode_population.m):
%   data_trial : 1 x n_cue cell; data_trial{c} is [n_stim*n_trials x N] counts,
%                blocked by stimulus (rows (s-1)*n_trials+1 : s*n_trials).
%   labels     : [n_stim*n_trials x 1] stimulus labels for that row order.
%   diagnostics: latents, weights, means, per-neuron r and cv, noise-free
%                signal PR per cue, cross-cue latent covariance spectrum,
%                positivity checks and the cfg used.
%
% USAGE:
%   [data_trial, labels, diag] = generate_partial_sharing_trial_data();     % alpha = 1 (K = 5)
%   cfg = struct('N', 50, 'n_stim', 30, 'n_latent', 5, 'n_shared', 2, 'seed', 7);   % alpha = 0.4
%   [data_trial, labels, diag] = generate_partial_sharing_trial_data(cfg);
%}

function [data_trial, labels, diagnostics] = generate_partial_sharing_trial_data(cfg)

%% Configuration (defaults; override any field via the cfg input struct)
if nargin < 1 || isempty(cfg), cfg = struct(); end
cfg = merge_defaults(cfg, default_config());
rng(cfg.seed, cfg.rng_algorithm);
N      = cfg.N;
n_stim = cfg.n_stim;
n_cue  = cfg.n_cue;
L      = cfg.n_latent;
Ks     = cfg.n_shared;
Kc     = L - Ks;                                % cue-specific dimensions per cue
assert(Ks >= 0 && Ks <= L && Ks == round(Ks), 'n_shared must be an integer in 0..n_latent (got %g).', Ks);
assert(n_stim > L, 'Need n_stim > n_latent (got %d <= %d).', n_stim, L);

%% Per-neuron statistics (same draws, in the same order, as Simulations_v2 Models 1 and 3)
r  = cfg.rate_median * exp(cfg.rate_log_sd * randn(1, N));   % [1 x N] mean count
cv = cfg.cv_median   * exp(cfg.cv_log_sd   * randn(1, N));   % [1 x N] tuning depth (CV)

%% Latent draw (cue 1's latent; its first K_s columns are the shared shape latent)
X = randn(n_stim, L);

%% Per-cue readout directions (unit vectors in latent space; rho-mixed across cues)
G_shared = randn(N, L);
U = cell(1, n_cue);
for c = 1:n_cue
    G = sqrt(cfg.rho) * G_shared + sqrt(1 - cfg.rho) * randn(N, L);
    U{c} = G ./ sqrt(sum(G.^2, 2));                     % [N x L] unit directions
end

%% Per-cue latents: shared columns + independent cue-specific columns, whitened per cue
Z = cell(1, n_cue);
Z{1} = whiten(X);
for c = 2:n_cue
    if Kc > 0
        Z{c} = whiten([X(:, 1:Ks), randn(n_stim, Kc)]);
        % Cholesky whitening leaves the leading K_s columns equal to cue 1's up to
        % rounding; make the sharing exact
        assert(Ks == 0 || max(abs(Z{c}(:, 1:Ks) - Z{1}(:, 1:Ks)), [], 'all') < 1e-10, ...
            'Shared latent columns differ across cues.');
        Z{c}(:, 1:Ks) = Z{1}(:, 1:Ks);
    else
        Z{c} = Z{1};
    end
end

%% Tuning-depth cap that keeps every rate positive without an offset or a floor
% lowest rate of neuron i = r_i * (1 - cv_i * m_i),  m_i = max_{s,c} (-u_{i,c}' z_{s,c})
switch lower(cfg.cap_rule)
    case 'exact'                                         % each neuron's own worst stimulus
        m = zeros(1, N);
        for c = 1:n_cue
            m = max(m, max(-(Z{c} * U{c}.'), [], 1));
        end
    case 'radius'                                        % conservative bound (v2's first run)
        m = zeros(1, N);
        for c = 1:n_cue
            m = max(m, max(sqrt(sum(Z{c}.^2, 2))));
        end
    otherwise
        error('Unknown cap_rule ''%s'' (use ''exact'' or ''radius'').', cfg.cap_rule);
end
cv_cap   = min(cfg.cv_max, (1 - cfg.min_rate_frac) ./ max(m, eps));   % [1 x N]
n_capped = sum(cv > cv_cap);
cv = min(cv, cv_cap);

%% Mean rates
W      = cell(1, n_cue);
lambda = cell(1, n_cue);
n_neg  = 0;
for c = 1:n_cue
    W{c} = (r .* cv).' .* U{c};                         % [N x L], row norm = r_i*cv_i
    lambda{c} = r + Z{c} * W{c}.';                      % [n_stim x N]
    n_neg = n_neg + sum(lambda{c}(:) < -1e-9);
    lambda{c} = max(lambda{c}, 0);                      % removes rounding error only
end
if n_neg > 0
    warning('generate_partial_sharing_trial_data: %d negative rates were floored (should be 0 by construction).', n_neg);
end

%% Poisson trials (blocked by stimulus)
data_trial = cell(1, n_cue);
for c = 1:n_cue
    data_trial{c} = poissrnd(repelem(lambda{c}, cfg.n_trials, 1));
end
labels = repelem((1:n_stim).', cfg.n_trials);

%% Diagnostics
diagnostics = struct();
diagnostics.Z              = Z;                            % 1 x n_cue, each [n_stim x L] whitened latents
diagnostics.W              = W;                            % 1 x n_cue, each [N x L]
diagnostics.lambda         = lambda;                       % 1 x n_cue, each [n_stim x N] (noise-free means)
diagnostics.rate           = r;                            % [1 x N]
diagnostics.cv             = cv;                           % [1 x N] (after cap)
diagnostics.cv_realized    = cellfun(@(L_) std(L_, 0, 1) ./ mean(L_, 1), lambda, 'UniformOutput', false);
diagnostics.signal_PR      = cellfun(@signal_participation_ratio, lambda);
diagnostics.mean_count     = mean(cellfun(@(x) mean(x(:)), data_trial));
diagnostics.min_lambda     = min(cellfun(@(L_) min(L_(:)), lambda));
diagnostics.frac_floored   = n_neg / (n_cue * n_stim * N);
diagnostics.frac_cv_capped = n_capped / N;
diagnostics.cv_cap         = cv_cap;                       % [1 x N] per-neuron depth cap
diagnostics.n_latent       = L;
diagnostics.n_shared       = Ks;
diagnostics.alpha          = Ks / L;
% singular values of the cue-1 / cue-2 latent cross-covariance: K_s ones (shared)
% and K - K_s chance-level values (cue-specific, ~1/sqrt(n_stim))
diagnostics.latent_cross_sv = svd(Z{1}.' * Z{min(2, n_cue)} / (n_stim - 1)).';
diagnostics.cfg            = cfg;

end  % main function


%% ----------------------------------------------------------------------
%% Local functions
%% ----------------------------------------------------------------------
function cfg = default_config()
% Defaults: Simulations_v2 Model 3 with a 5-D latent (sweep_N_S_model3_D5), fully shared.
cfg = struct();
cfg.N             = 100;        % number of neurons (all used by the driver)
cfg.n_cue         = 2;          % renderings: cue 1 = "ac", cue 2 = "ec"
cfg.n_stim        = 50;         % stimuli (swept by the driver)
cfg.n_trials      = 10;         % trials per stimulus (as recorded)
cfg.n_latent      = 5;          % K: latent dimensions per cue (fixed across alpha)
cfg.n_shared      = [];         % K_s: shared latent dimensions ([] = n_latent, i.e. alpha = 1)
cfg.rate_median   = 9.4;        % median mean count per neuron (330-630 ms window; calibrate_from_real_data.m)
cfg.rate_log_sd   = 0.66;       % log-SD of mean count across neurons
cfg.cv_median     = 0.14;       % median tuning depth (noise-corrected CV of the tuning curve)
cfg.cv_log_sd     = 0.65;       % log-SD of tuning depth across neurons
cfg.cv_max        = 0.35;       % fixed cap on tuning depth (same as v2 Models 1 and 3)
cfg.min_rate_frac = 0;          % lowest allowed rate as a fraction of r_i (0 = rate may touch 0)
cfg.cap_rule      = 'exact';    % depth cap: 'exact' per-neuron bound, or 'radius'
cfg.rho           = 0;          % cross-cue readout correlation (0 = independent)
cfg.seed          = 1;          % per-population RNG seed
cfg.rng_algorithm = 'twister';  % global-stream algorithm
end

function cfg = merge_defaults(cfg, def)
% Overlay defaults for any field the user did not provide.
f = fieldnames(def);
for k = 1:numel(f)
    if ~isfield(cfg, f{k}) || isempty(cfg.(f{k}))
        cfg.(f{k}) = def.(f{k});
    end
end
if isempty(cfg.n_shared), cfg.n_shared = cfg.n_latent; end
end

function Z = whiten(X)
% Zero mean and identity covariance across stimuli (Cholesky: column j of the
% result depends only on columns 1..j of X). Same operations as v2 Model 3.
Z = X - mean(X, 1);
Z = Z / chol(cov(Z));
end

function pr = signal_participation_ratio(mu)
% Participation ratio of the noise-free signal covariance across stimuli.
mu_c = mu - mean(mu, 1);
s    = svd(mu_c, 'econ');
ev   = s.^2;
pr   = sum(ev)^2 / sum(ev.^2);
end
