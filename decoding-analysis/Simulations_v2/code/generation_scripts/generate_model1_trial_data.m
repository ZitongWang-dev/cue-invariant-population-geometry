%{
% FILENAME: generate_model1_trial_data.m   (Simulations_v2)
% AUTHOR:   Zitong Wang
% DATE:     2026-09-17
%
% DESCRIPTION:
%   Model 1 (null) generator, version 2. Each neuron's mean response to each
%   stimulus, in each cue, is an INDEPENDENT random draw: the reviewer's
%   "simplest case", with no structure shared across cues by construction.
%
%   What changed from Simulations/ (v1), and why:
%     v1 drew Gaussian deviations of SD `weight_scale` around one large constant
%     baseline (ceil(5.5 * weight_scale)) chosen to keep every rate positive.
%     The offset was large relative to the signal, so tuning was shallow, and
%     the operating point was set indirectly (count and depth both moved with
%     weight_scale). v2 sets the single-neuron statistics DIRECTLY and matches
%     them to the recorded populations (see calibrate_from_real_data.m):
%       r_i   mean spike count of neuron i in the window     (log-normal across neurons)
%       cv_i  tuning depth of neuron i = SD of its mean response across stimuli
%             divided by its mean response                    (log-normal across neurons)
%     and draws the per-stimulus means from a Gamma distribution with exactly
%     that mean and CV:  lambda_{i,s,c} ~ Gamma(shape = 1/cv_i^2, scale = r_i*cv_i^2).
%     Gamma means are strictly positive, so no positivity offset, no floor and
%     no rectification are ever needed, and the tuning curves are right-skewed
%     like real ones. Poisson trials are then drawn around each mean.
%
%   Generative model (neuron i, stimulus s, cue c, trial t):
%       r_i    = rate_median * exp(rate_log_sd * xi_i)        xi  ~ N(0,1)
%       cv_i   = cv_median   * exp(cv_log_sd   * eta_i)       eta ~ N(0,1), capped at cv_max
%       lambda_{i,s,c} ~ Gamma(1/cv_i^2, r_i*cv_i^2)          i.i.d. over s and c
%       count_{i,s,c,t} ~ Poisson(lambda_{i,s,c})
%   r_i and cv_i are properties of the neuron and are shared across cues; the
%   per-stimulus means are not. After the pipeline's per-neuron z-scoring the
%   two cues are therefore independent random configurations.
%
%   rho (interface parity with Model 3): fraction of each mean shared across
%   cues, implemented through Gamma additivity, lambda_c = G_shared + G_indep_c
%   with shapes rho/cv^2 and (1-rho)/cv^2, so every lambda_c keeps exactly the
%   same Gamma(1/cv^2, r*cv^2) marginal and the cross-cue correlation is rho.
%   rho = 0 (default) is the true null; rho = 1 gives identical means across
%   cues (trivial transfer).
%
%   The RNG is consumed in the order r, cv, then the model-specific draws, in
%   both v2 generators, so Model 1 and Model 3 populations with the same seed
%   share the same neurons (same r_i, cv_i) and differ only in geometry.
%
%   Defaults are the pooled empirical operating point of the four recorded
%   populations (FR/KO x V1/V2, 330-630 ms, all three cues; see
%   calibrate_from_real_data.m): median count 9.4, log-SD 0.66 across
%   neurons; median noise-corrected tuning CV 0.14, log-SD 0.65; Fano
%   1.08 (Poisson adequate). cv_max caps the log-normal tail at the same
%   value in both models (Model 3 needs a cap for positivity; applying the same
%   one here keeps the single-neuron SNR distributions identical).
%
% OUTPUT CONTRACT (shared with generate_model3_trial_data.m and the driver):
%   data_trial : 1 x n_cue cell; data_trial{c} is [n_stim*n_trials x N] counts,
%                blocked by stimulus (rows (s-1)*n_trials+1 : s*n_trials).
%   labels     : [n_stim*n_trials x 1] stimulus labels for that row order.
%   diagnostics: struct with the drawn means, per-neuron r and cv, noise-free
%                signal PR per cue, mean count, and the cfg used.
%
% USAGE:
%   [data_trial, labels, diag] = generate_model1_trial_data();            % defaults
%   cfg = struct('N', 50, 'n_stim', 20, 'seed', 7);
%   [data_trial, labels, diag] = generate_model1_trial_data(cfg);
%}

function [data_trial, labels, diagnostics] = generate_model1_trial_data(cfg)

%% Configuration (defaults; override any field via the cfg input struct)
if nargin < 1 || isempty(cfg), cfg = struct(); end
cfg = merge_defaults(cfg, default_config());
rng(cfg.seed, cfg.rng_algorithm);
N      = cfg.N;
n_stim = cfg.n_stim;

%% Per-neuron statistics (properties of the neuron: shared across cues)
r  = cfg.rate_median * exp(cfg.rate_log_sd * randn(1, N));   % [1 x N] mean count
cv = cfg.cv_median   * exp(cfg.cv_log_sd   * randn(1, N));   % [1 x N] tuning depth (CV)
n_capped = sum(cv > cfg.cv_max);
cv = min(cv, cfg.cv_max);
shape = 1 ./ cv.^2;             % Gamma shape: mean = shape*scale, CV = 1/sqrt(shape)
scale = r .* cv.^2;             % Gamma scale

%% Per-cue random means (Gamma, i.i.d. over stimuli; rho-shared across cues)
lambda = cell(1, cfg.n_cue);
if cfg.rho > 0
    G_shared = gamrnd(repmat(cfg.rho * shape, n_stim, 1), repmat(scale, n_stim, 1));
else
    G_shared = zeros(n_stim, N);
end
for c = 1:cfg.n_cue
    if cfg.rho < 1
        G_indep = gamrnd(repmat((1 - cfg.rho) * shape, n_stim, 1), repmat(scale, n_stim, 1));
    else
        G_indep = zeros(n_stim, N);
    end
    lambda{c} = G_shared + G_indep;             % [n_stim x N], strictly positive
end

%% Poisson trials (blocked by stimulus)
data_trial = cell(1, cfg.n_cue);
for c = 1:cfg.n_cue
    data_trial{c} = poissrnd(repelem(lambda{c}, cfg.n_trials, 1));
end
labels = repelem((1:n_stim).', cfg.n_trials);

%% Diagnostics
diagnostics = struct();
diagnostics.lambda        = lambda;                        % 1 x n_cue, each [n_stim x N] (noise-free means)
diagnostics.rate          = r;                             % [1 x N] drawn mean count
diagnostics.cv            = cv;                            % [1 x N] drawn tuning depth (after cap)
diagnostics.cv_realized   = cellfun(@(L) std(L, 0, 1) ./ mean(L, 1), lambda, 'UniformOutput', false);
diagnostics.signal_PR     = cellfun(@signal_participation_ratio, lambda);   % HIGH (near full rank)
diagnostics.mean_count    = mean(cellfun(@(x) mean(x(:)), data_trial));
diagnostics.min_lambda    = min(cellfun(@(L) min(L(:)), lambda));
diagnostics.frac_floored  = 0;                             % Gamma means are positive by construction
diagnostics.frac_cv_capped = n_capped / N;
diagnostics.cfg           = cfg;

end  % main function


%% ----------------------------------------------------------------------
%% Local functions
%% ----------------------------------------------------------------------
function cfg = default_config()
% Default Model 1 (v2) configuration. All fields overridable via the cfg input.
cfg = struct();
cfg.N             = 100;        % number of neurons (all used by the driver)
cfg.n_cue         = 2;          % renderings: cue 1 = "ac", cue 2 = "ec"
cfg.n_stim        = 50;         % stimuli (swept by the driver)
cfg.n_trials      = 10;         % trials per stimulus (as recorded)
cfg.rate_median   = 9.4;        % median mean count per neuron (330-630 ms window; calibrate_from_real_data.m)
cfg.rate_log_sd   = 0.66;       % log-SD of mean count across neurons
cfg.cv_median     = 0.14;       % median tuning depth (noise-corrected CV of the tuning curve)
cfg.cv_log_sd     = 0.65;       % log-SD of tuning depth across neurons
cfg.cv_max        = 0.35;       % cap on tuning depth (same in both models; see header)
cfg.rho           = 0;          % fraction of each mean shared across cues (0 = independent)
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
% mu : [n_stim x N] mean rates, centered across stimuli. Full-rank random -> HIGH.
mu_c = mu - mean(mu, 1);
s    = svd(mu_c, 'econ');
ev   = s.^2;
pr   = sum(ev)^2 / sum(ev.^2);
end
