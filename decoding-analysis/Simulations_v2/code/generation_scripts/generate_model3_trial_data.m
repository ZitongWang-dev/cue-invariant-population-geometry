%{
% FILENAME: generate_model3_trial_data.m   (Simulations_v2)
% AUTHOR:   Zitong Wang
% DATE:     2026-09-17
%
% DESCRIPTION:
%   Model 3 (linear mixed selectivity) generator, version 2. Both cues share
%   one fixed low-dimensional stimulus geometry (two latent shape parameters
%   per stimulus, the cue-invariant contour geometry) that each neuron reads
%   out LINEARLY with its own weights, drawn independently per cue. The
%   reviewer's "more interesting" case: mean representations that are a linear
%   combination of contour-shape parameters shared across renderings.
%
%   Generative model (neuron i, stimulus s, cue c, trial t):
%       z_s            shared latent, [n_stim x n_latent], whitened across stimuli
%                      (zero mean, identity covariance), drawn once per population
%       r_i    = rate_median * exp(rate_log_sd * xi_i)        mean count     (as Model 1)
%       cv_i   = cv_median   * exp(cv_log_sd   * eta_i)       tuning depth   (as Model 1)
%       u_{i,c}        unit readout direction in latent space, random per neuron
%                      and cue (rho-mixed across cues, see below)
%       w_{i,c} = r_i * cv_i * u_{i,c}
%       lambda_{i,s,c} = r_i + w_{i,c} . z_s
%       count_{i,s,c,t} ~ Poisson(lambda_{i,s,c})
%   Because z is whitened and |u| = 1, the SD of neuron i's tuning curve across
%   stimuli is exactly r_i*cv_i in every cue, i.e. the same per-neuron tuning
%   depth and mean count as Model 1 with the same seed. The two models then
%   differ ONLY in geometry: rank-2 shared latents here versus independent
%   full-rank random means there.
%
%   What changed from Simulations/ (v1), and why:
%     v1 used a fixed 10 x 5 factorial grid of latents (only valid for 50
%     stimuli), Gaussian weights of SD `weight_scale`, and one large constant
%     baseline (ceil(5.5 * weight_scale * M_grid)) for positivity, which made
%     tuning shallow and forced absurd mean counts (~140) to reach a data-like
%     self-decoding. v2 (a) draws the latent configuration randomly so any
%     n_stim works (a grid layout is still available for n_stim = 50 parity),
%     (b) sets each neuron's mean count and tuning depth directly from the
%     recorded distributions (calibrate_from_real_data.m), and (c) keeps rates
%     positive by bounding each neuron's tuning depth instead of adding an
%     offset. Writing lambda_{i,s,c} = r_i * (1 + cv_i * u_{i,c}' z_s), the
%     lowest rate of neuron i is r_i * (1 - cv_i * m_i) with
%     m_i = max_{s,c} (-u_{i,c}' z_s), so
%         cv_i <= min( cv_max , (1 - min_rate_frac) / m_i )
%     keeps every rate >= r_i * min_rate_frac with no floor and no
%     rectification (rectification would curve the manifold and inflate PR,
%     turning Model 3 into a nonlinear model). cv_max is the same fixed cap as
%     Model 1, and in practice it is the binding one (m_i ~ 2.3-2.8 for 50
%     whitened stimuli), so Model 3 neurons keep Model 1's single-neuron
%     statistics at any latent dimension. The capped fraction is reported.
%     cap_rule = 'radius' reproduces the earlier, more conservative bound
%     cv_i <= 1 / max_s |z_s|, which tightens as n_latent grows (it capped 15%
%     of neurons at 0.28 for n_latent = 5) and was used for the first sweep_N_S
%     run (n_latent = 2).
%
%   rho: cross-cue readout correlation. The readout directions are
%   u_c = normalize(sqrt(rho)*g_shared + sqrt(1-rho)*g_indep_c), g ~ N(0, I).
%   rho = 0 (default): independent readouts per cue, the Model 3 base case.
%   rho = 1: identical readouts (identical means) across cues, trivial transfer.
%
%   Expected behaviour: signal PR = n_latent (2) exactly (noise-free), shared
%   structure across cues, so transfer decoding should beat the
%   correspondence-shuffled control whenever the rotation is not free enough
%   to overfit (N small relative to n_stim), and the control should approach
%   PT as N grows past n_stim.
%
% OUTPUT CONTRACT (shared with generate_model1_trial_data.m and the driver):
%   data_trial : 1 x n_cue cell; data_trial{c} is [n_stim*n_trials x N] counts,
%                blocked by stimulus (rows (s-1)*n_trials+1 : s*n_trials).
%   labels     : [n_stim*n_trials x 1] stimulus labels for that row order.
%   diagnostics: struct with the latents, weights, means, per-neuron r and cv,
%                noise-free signal PR per cue, mean count, positivity checks
%                and the cfg used.
%
% USAGE:
%   [data_trial, labels, diag] = generate_model3_trial_data();            % defaults
%   cfg = struct('N', 50, 'n_stim', 20, 'seed', 7);
%   [data_trial, labels, diag] = generate_model3_trial_data(cfg);
%   cfg = struct('latent_layout', 'grid', 'n_levels', [10 5]);           % v1 grid, 50 stimuli
%}

function [data_trial, labels, diagnostics] = generate_model3_trial_data(cfg)

%% Configuration (defaults; override any field via the cfg input struct)
if nargin < 1 || isempty(cfg), cfg = struct(); end
cfg = merge_defaults(cfg, default_config());
rng(cfg.seed, cfg.rng_algorithm);
N      = cfg.N;
n_stim = cfg.n_stim;
L      = cfg.n_latent;

%% Per-neuron statistics (same draws, in the same order, as Model 1)
r  = cfg.rate_median * exp(cfg.rate_log_sd * randn(1, N));   % [1 x N] mean count
cv = cfg.cv_median   * exp(cfg.cv_log_sd   * randn(1, N));   % [1 x N] tuning depth (CV)

%% Shared latent geometry (one configuration per population, common to all cues)
switch lower(cfg.latent_layout)
    case 'random'
        Z = randn(n_stim, L);
    case 'grid'
        assert(numel(cfg.n_levels) == L && prod(cfg.n_levels) == n_stim, ...
            'latent_layout=''grid'' needs prod(n_levels) == n_stim (got %s vs %d).', ...
            mat2str(cfg.n_levels), n_stim);
        axes_ = arrayfun(@(n) linspace(-1, 1, n), cfg.n_levels, 'UniformOutput', false);
        G = cell(1, L);
        [G{:}] = ndgrid(axes_{:});
        Z = cell2mat(cellfun(@(g) g(:), G, 'UniformOutput', false));
    otherwise
        error('Unknown latent_layout ''%s'' (use ''random'' or ''grid'').', cfg.latent_layout);
end
assert(n_stim > L, 'Need n_stim > n_latent (got %d <= %d).', n_stim, L);
Z = Z - mean(Z, 1);
Z = Z / chol(cov(Z));                        % whiten: cov(Z) = I, so tuning SD = |w| exactly

%% Per-cue readout directions (unit vectors in latent space; rho-mixed across cues)
% Drawn before the depth cap (the cap needs them); the RNG order is unchanged.
G_shared = randn(N, L);
U = cell(1, cfg.n_cue);
for c = 1:cfg.n_cue
    G = sqrt(cfg.rho) * G_shared + sqrt(1 - cfg.rho) * randn(N, L);
    U{c} = G ./ sqrt(sum(G.^2, 2));                     % [N x L] unit directions
end

%% Tuning-depth cap that keeps every rate positive without an offset or a floor
% lowest rate of neuron i = r_i * (1 - cv_i * m_i),  m_i = max_{s,c} (-u_{i,c}' z_s)
switch lower(cfg.cap_rule)
    case 'exact'                                         % each neuron's own worst stimulus
        m = zeros(1, N);
        for c = 1:cfg.n_cue
            m = max(m, max(-(Z * U{c}.'), [], 1));
        end
    case 'radius'                                        % earlier conservative bound (first sweep_N_S run)
        m = repmat(max(sqrt(sum(Z.^2, 2))), 1, N);
    otherwise
        error('Unknown cap_rule ''%s'' (use ''exact'' or ''radius'').', cfg.cap_rule);
end
cv_cap   = min(cfg.cv_max, (1 - cfg.min_rate_frac) ./ max(m, eps));   % [1 x N]
n_capped = sum(cv > cv_cap);
cv = min(cv, cv_cap);

%% Mean rates
W      = cell(1, cfg.n_cue);
lambda = cell(1, cfg.n_cue);
n_neg  = 0;
for c = 1:cfg.n_cue
    W{c} = (r .* cv).' .* U{c};                         % [N x L], row norm = r_i*cv_i
    lambda{c} = r + Z * W{c}.';                         % [n_stim x N]
    n_neg = n_neg + sum(lambda{c}(:) < -1e-9);
    lambda{c} = max(lambda{c}, 0);                      % removes rounding error only
end
if n_neg > 0
    warning('generate_model3_trial_data: %d negative rates were floored (should be 0 by construction).', n_neg);
end

%% Poisson trials (blocked by stimulus)
data_trial = cell(1, cfg.n_cue);
for c = 1:cfg.n_cue
    data_trial{c} = poissrnd(repelem(lambda{c}, cfg.n_trials, 1));
end
labels = repelem((1:n_stim).', cfg.n_trials);

%% Diagnostics
diagnostics = struct();
diagnostics.Z             = Z;                             % [n_stim x L] shared whitened latents
diagnostics.W             = W;                             % 1 x n_cue, each [N x L]
diagnostics.lambda        = lambda;                        % 1 x n_cue, each [n_stim x N] (noise-free means)
diagnostics.rate          = r;                             % [1 x N]
diagnostics.cv            = cv;                            % [1 x N] (after cap)
diagnostics.cv_realized   = cellfun(@(L_) std(L_, 0, 1) ./ mean(L_, 1), lambda, 'UniformOutput', false);
diagnostics.signal_PR     = cellfun(@signal_participation_ratio, lambda);   % ~ n_latent
diagnostics.mean_count    = mean(cellfun(@(x) mean(x(:)), data_trial));
diagnostics.min_lambda    = min(cellfun(@(L_) min(L_(:)), lambda));
diagnostics.frac_floored  = n_neg / (cfg.n_cue * n_stim * N);
diagnostics.frac_cv_capped = n_capped / N;
diagnostics.cv_cap        = cv_cap;                        % [1 x N] per-neuron depth cap
diagnostics.n_latent      = L;
diagnostics.cfg           = cfg;

end  % main function


%% ----------------------------------------------------------------------
%% Local functions
%% ----------------------------------------------------------------------
function cfg = default_config()
% Default Model 3 (v2) configuration. All fields overridable via the cfg input.
cfg = struct();
cfg.N             = 100;        % number of neurons (all used by the driver)
cfg.n_cue         = 2;          % renderings: cue 1 = "ac", cue 2 = "ec"
cfg.n_stim        = 50;         % stimuli (swept by the driver)
cfg.n_trials      = 10;         % trials per stimulus (as recorded)
cfg.n_latent      = 2;          % shared latent shape parameters (rank of the signal)
cfg.latent_layout = 'random';   % 'random' (any n_stim) or 'grid' (needs n_levels, v1 parity)
cfg.n_levels      = [10 5];     % grid layout only: levels per latent axis, prod == n_stim
cfg.rate_median   = 9.4;        % median mean count per neuron (330-630 ms window; calibrate_from_real_data.m)
cfg.rate_log_sd   = 0.66;       % log-SD of mean count across neurons
cfg.cv_median     = 0.14;       % median tuning depth (noise-corrected CV of the tuning curve)
cfg.cv_log_sd     = 0.65;       % log-SD of tuning depth across neurons
cfg.cv_max        = 0.35;       % fixed cap on tuning depth (same as Model 1)
cfg.min_rate_frac = 0;          % lowest allowed rate as a fraction of r_i (0 = rate may touch 0)
cfg.cap_rule      = 'exact';    % depth cap: 'exact' per-neuron bound, or 'radius' (first sweep_N_S run)
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
end

function pr = signal_participation_ratio(mu)
% Participation ratio of the noise-free signal covariance across stimuli.
% mu : [n_stim x N] mean rates, centered across stimuli. Rank-L -> pr = L.
mu_c = mu - mean(mu, 1);
s    = svd(mu_c, 'econ');
ev   = s.^2;
pr   = sum(ev)^2 / sum(ev.^2);
end
