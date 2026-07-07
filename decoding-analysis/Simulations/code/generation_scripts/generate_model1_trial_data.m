%{
% FILENAME: generate_model1_trial_data.m
% AUTHOR:   Zitong Wang
% DATE:     2026-07-01
%
% DESCRIPTION:
%   Generator for the Model 1 (null) benchmark of the Procrustes transfer-
%   decoding pipeline. Each neuron's mean response to each stimulus, in each
%   rendering, is an INDEPENDENT random draw -- there is no shared latent
%   structure across renderings. This is the reviewer's "simplest case": mean
%   representations chosen randomly and independently for each stimulus class
%   and each rendering style.
%
%   Generative model (per neuron i, stimulus s, cue c):
%       lambda_{i,s,c} = baseline_i + signal_{i,s,c}
%       count_{i,s,c,t} ~ Poisson(lambda_{i,s,c}),   t = 1..n_trials
%   with signal_{i,s,c} ~ N(0, weight_scale^2) drawn independently per element.
%
%   Expected result: self-decoding is high (50 random patterns in N-D are
%   near-orthogonal, so trivially separable), but PT is at CHANCE -- the two
%   renderings' point clouds are independent random configurations, so no
%   similarity transform maps stimulus s of one cue onto stimulus s of the
%   other. This confirms that transfer requires shared structure; PT working
%   anywhere else in the family is therefore not trivial. signal_PR is HIGH
%   (near full rank, ~30-35 for 50 stimuli in 100 neurons) -- the opposite
%   extreme from Model 3's PR ~ 2.
%
%   Interface parity with the shared driver / Model 3:
%   - weight_scale : per-element signal SD. (Kept as "weight_scale" so the
%     driver's build_gen_cfg passes it unchanged; Model 1 has no readout
%     weights -- it is simply the magnitude of the random means.)
%   - rho          : fraction of the random means SHARED across cues.
%     rho = 0 (default) = independent means = the true Model 1 (PT at chance).
%     rho -> 1 = identical means across cues (trivial transfer). Interface only.
%   - baseline     : additive, constant across stimuli (pure translation,
%     removed by z-scoring and by the Procrustes translation term). 'auto' by
%     default. Because the signal is i.i.d. with per-element SD = weight_scale,
%     the auto margin uses M_eff = 1 (no grid factor, unlike Model 3).
%
%   Rectification is deliberately avoided (see Model 3): a large-enough baseline
%   keeps lambda >= 0 without half-wave rectifying; a floor is a safety net only
%   and WARNS if it engages.
%
% OUTPUT CONTRACT (identical to generate_model3_trial_data.m):
%   data_trial : 1 x n_cue cell; data_trial{c} is a [500 x N] double of spike
%                counts, blocked by stimulus in groups of n_trials rows.
%   labels     : [500 x 1] stimulus labels matching the blocked row order.
%   diagnostics: struct with the drawn means, achieved signal PR per cue, mean
%                count, and positivity checks.
%
% USAGE:
%   [data_trial, labels, diag] = generate_model1_trial_data();          % defaults
%   cfg.seed = 42; cfg.weight_scale = 3;
%   [data_trial, labels, diag] = generate_model1_trial_data(cfg);
%}

function [data_trial, labels, diagnostics] = generate_model1_trial_data(cfg)

%% Configuration (defaults; override any field via the cfg input struct)
if nargin < 1 || isempty(cfg), cfg = struct(); end
cfg = merge_defaults(cfg, default_config());

% One global-stream seed per population (see Model 3): reproducible per seed.
rng(cfg.seed, cfg.rng_algorithm);

n_stim = cfg.n_stim;

%% Per-cue random mean signal (i.i.d. -- the null: no shared structure)
% signal_c = sqrt(rho)*S_shared + sqrt(1-rho)*S_indep_c
%   rho = 0 : independent random means per cue (true Model 1) -> PT at chance.
%   rho -> 1: identical means across cues (trivial transfer). Interface parity.
S_shared = cfg.weight_scale * randn(n_stim, cfg.N);
signal = cell(1, cfg.n_cue);
for c = 1:cfg.n_cue
    S_indep   = cfg.weight_scale * randn(n_stim, cfg.N);
    signal{c} = sqrt(cfg.rho) * S_shared + sqrt(1 - cfg.rho) * S_indep;   % [n_stim x N]
end

%% Baseline (constant across stimuli -> pure translation of the cloud)
% i.i.d. signal: per-element SD = weight_scale, so M_eff = 1 (no grid factor).
M_eff     = 1;
safe_base = cfg.baseline_k * cfg.weight_scale * M_eff;   % recommended margin

if (ischar(cfg.baseline) || isstring(cfg.baseline)) && strcmpi(cfg.baseline, 'auto')
    baseline = ceil(safe_base);
else
    baseline = cfg.baseline;
    if ~(isscalar(baseline) || numel(baseline) == cfg.N)
        error('cfg.baseline must be a scalar, a length-%d vector, or ''auto''.', cfg.N);
    end
    % A-priori check: a manual baseline below the safe margin risks flooring.
    if min(baseline) < safe_base
        warning(['generate_model1_trial_data: cfg.baseline (min %.3g) is below the ' ...
            'recommended margin %.3g (= %.1f * weight_scale) and may floor ' ...
            'negative rates. Consider cfg.baseline = ''auto'' or >= %g.'], ...
            min(baseline), safe_base, cfg.baseline_k, ceil(safe_base));
    end
end
baseline = baseline(:).';                    % 1x1 or 1xN row for broadcasting

%% Build rates, check positivity, draw Poisson trials
data_trial = cell(1, cfg.n_cue);
lambda_all = cell(1, cfg.n_cue);
n_neg = 0; n_tot = 0;
for c = 1:cfg.n_cue
    lambda = signal{c} + baseline;           % broadcast baseline over stimuli
    n_neg  = n_neg + sum(lambda(:) < 0);
    n_tot  = n_tot + numel(lambda);
    lambda = max(lambda, cfg.lambda_floor);  % safety net only (see warning)
    lambda_all{c} = lambda;

    lambda_expanded = repelem(lambda, cfg.n_trials, 1);   % [n_stim*n_trials x N]
    data_trial{c}   = poissrnd(lambda_expanded);
end

labels = repelem((1:n_stim).', cfg.n_trials);            % [500 x 1]

%% Positivity warning (baseline too low => rectification)
frac_floored = n_neg / n_tot;
if frac_floored > cfg.floor_warn_tol
    min_signal   = min(cellfun(@(s) min(s(:)), signal));
    rec_baseline = ceil(-min_signal + cfg.floor_margin_sd * cfg.weight_scale);
    warning(['generate_model1_trial_data: %.3f%% of rates were negative and ' ...
        'floored. Raise cfg.baseline (suggest >= %g) or lower cfg.weight_scale.'], ...
        100 * frac_floored, rec_baseline);
end

%% Diagnostics
diagnostics = struct();
diagnostics.signal       = signal;                       % 1 x n_cue, each [n_stim x N] (means)
diagnostics.lambda       = lambda_all;                   % 1 x n_cue, each [n_stim x N]
diagnostics.signal_PR    = cellfun(@signal_participation_ratio, lambda_all);  % HIGH (~30-35)
diagnostics.mean_count   = mean(cellfun(@(x) mean(x(:)), data_trial));
diagnostics.min_lambda   = min(cellfun(@(L) min(L(:)), lambda_all));
diagnostics.frac_floored = frac_floored;
diagnostics.cfg          = cfg;

end  % main function


%% ----------------------------------------------------------------------
%% Local functions
%% ----------------------------------------------------------------------
function cfg = default_config()
% Default Model 1 configuration. All fields overridable via the cfg input.
cfg = struct();
cfg.N              = 100;         % number of neurons (all used)
cfg.n_cue          = 2;          % renderings: cue 1 = "ac", cue 2 = "ec"
cfg.n_stim         = 50;         % stimuli
cfg.n_trials       = 10;         % trials per stimulus (pipeline expects 10)
cfg.weight_scale   = 3;          % per-element signal SD (placeholder; calibrate to target self)
cfg.rho            = 0;          % fraction of means shared across cues (0 = independent)
cfg.baseline       = 'auto';     % 'auto' = ceil(baseline_k*weight_scale); or scalar / [N x 1]
cfg.baseline_k     = 5.5;        % margin factor for 'auto' (extreme-value z* + slack)
cfg.lambda_floor   = 0;          % safety floor on rates (should never engage)
cfg.floor_warn_tol = 1e-3;       % warn if >0.1% of rates get floored
cfg.floor_margin_sd = 5;         % SDs of margin used in the baseline suggestion
cfg.seed           = 1;          % per-population RNG seed
cfg.rng_algorithm  = 'twister';  % global-stream algorithm
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
% mu : [n_stim x N] mean rates. Baseline is constant across stimuli and removed
% by the centering, so PR(lambda) = PR(signal). Full-rank random -> HIGH.
mu_c = mu - mean(mu, 1);          % center across stimuli
s    = svd(mu_c, 'econ');         % singular values
ev   = s.^2;                      % covariance eigenvalues (up to scale)
pr   = sum(ev)^2 / sum(ev.^2);
end
