%{
% FILENAME: generate_model3_trial_data.m
% AUTHOR:   Zitong Wang
% DATE:     2026-07-01
%
% DESCRIPTION:
%   Generator for the Model 3 (linear mixed selectivity) benchmark of the
%   Procrustes transfer-decoding pipeline. Produces synthetic trial-by-trial
%   spike counts for TWO renderings (cue 1 = "ac", cue 2 = "ec") that share a
%   common, fixed low-dimensional stimulus geometry but are read out by
%   independent (or, optionally, correlated) linear weights.
%
%   Generative model (per neuron i, stimulus s, cue c):
%       lambda_{i,s,c} = baseline_i + w1_{i,c} * z1_s + w2_{i,c} * z2_s
%       count_{i,s,c,t} ~ Poisson(lambda_{i,s,c}),   t = 1..n_trials
%
%   - z1_s, z2_s : two standardized linear shape parameters, SHARED across
%     cues (the cue-invariant contour geometry). Defined once on a fixed
%     10 x 5 full-factorial grid; deterministic, identical across populations
%     and models.
%   - W_c = [w1 w2] : [N x 2] readout weights drawn per cue. Signed,
%     zero-mean Gaussian. Independent across cues by default (rho = 0);
%     cfg.rho in (0,1] induces element-wise cross-cue correlation for later
%     experiments (rho = 1 -> shared readout / trivial transfer).
%   - baseline : additive, constant across stimuli, so it is a pure
%     translation of the cloud (removed by z-scoring and by the Procrustes
%     translation term). It sets the mean count / SNR without touching the
%     rank-2 signal geometry. Defaults to 'auto' = ceil(baseline_k * weight_scale
%     * M_grid): a single fixed constant per run whose positivity margin tracks
%     weight_scale as you calibrate SNR. Also accepts a manual scalar, or an
%     [N x 1] vector for per-neuron rates.
%
%   Rectification is deliberately AVOIDED: a large-enough baseline keeps
%   lambda >= 0 without half-wave rectifying the signal, which would curve the
%   manifold and inflate the participation ratio (turning Model 3 into a
%   Model-6-like model). A floor is applied only as a safety net, and the
%   function WARNS if it engages, signalling that baseline is too low relative
%   to weight_scale.
%
% OUTPUT CONTRACT (matches the decoding helpers, e.g. the pair_pcaloader /
% take_average / sample_trial functions in procrustes_decoding_basic.m):
%   data_trial : 1 x n_cue cell. data_trial{c} is a [500 x N] double of spike
%                counts, blocked by stimulus in groups of n_trials rows
%                (rows (s-1)*10+1 : s*10 are the 10 trials of stimulus s).
%                Consumed directly by pair_pcaloader -- no cell2mat and no
%                multiclass_svmloader_PT step (that loader is only for the raw
%                ms-level neural recordings).
%   labels     : [500 x 1] stimulus labels matching the blocked row order.
%   diagnostics: struct with the fixed geometry, drawn weights, achieved
%                signal participation ratio per cue, mean count, and
%                positivity checks (for the isolated PR / calibration test).
%
% USAGE:
%   [data_trial, labels, diag] = generate_model3_trial_data();          % defaults
%   cfg.seed = 42; cfg.weight_scale = 0.8;
%   [data_trial, labels, diag] = generate_model3_trial_data(cfg);
%}

function [data_trial, labels, diagnostics] = generate_model3_trial_data(cfg)

%% Configuration (defaults; override any field via the cfg input struct)
if nargin < 1 || isempty(cfg), cfg = struct(); end
cfg = merge_defaults(cfg, default_config());

% Reproducibility: one global-stream seed per population. Called inside the
% driver's parfor with a distinct seed per iteration, this makes every
% population independently reproducible regardless of worker scheduling.
rng(cfg.seed, cfg.rng_algorithm);

%% Fixed, shared stimulus geometry (deterministic; consumes no RNG)
% Two linear shape axes on a full-factorial grid -> n_stim points, then
% standardize each axis to zero mean / unit variance across stimuli so the two
% latents contribute comparable signal variance (target signal PR ~ 2).
[g1, g2] = ndgrid(linspace(-1, 1, cfg.n_levels(1)), ...
                  linspace(-1, 1, cfg.n_levels(2)));
Phi = [g1(:), g2(:)];                       % [n_stim x 2]
Phi = (Phi - mean(Phi, 1)) ./ std(Phi, 0, 1);
n_stim = size(Phi, 1);
assert(n_stim == cfg.n_stim, ...
    'Grid n_levels %s yields %d stimuli, expected %d.', ...
    mat2str(cfg.n_levels), n_stim, cfg.n_stim);

%% Per-cue linear readout weights (signed, zero-mean Gaussian)
% W_c = sqrt(rho)*W_shared + sqrt(1-rho)*W_indep_c
%   rho = 0 (default): independent readouts per cue (the Model 3 base case).
%   rho -> 1         : shared readout (trivial transfer). Interface for later.
W_shared = cfg.weight_scale * randn(cfg.N, 2);
W = cell(1, cfg.n_cue);
for c = 1:cfg.n_cue
    W_indep = cfg.weight_scale * randn(cfg.N, 2);
    W{c} = sqrt(cfg.rho) * W_shared + sqrt(1 - cfg.rho) * W_indep;
end

%% Baseline (constant across stimuli -> pure translation of the cloud)
% M_grid = norm of the most extreme (corner) standardized grid point; the
% worst-case negative signal scales as ~ z* * weight_scale * M_grid.
M_grid    = max(sqrt(sum(Phi.^2, 2)));       % ~1.95 for the 10x5 grid
safe_base = cfg.baseline_k * cfg.weight_scale * M_grid;   % recommended margin

if (ischar(cfg.baseline) || isstring(cfg.baseline)) && strcmpi(cfg.baseline, 'auto')
    % One fixed constant for the whole run. Tracks weight_scale between runs so
    % the positivity margin stays correct during SNR calibration, while staying
    % identical across populations within a run (constant SNR).
    baseline = ceil(safe_base);
else
    baseline = cfg.baseline;
    if ~(isscalar(baseline) || numel(baseline) == cfg.N)
        error('cfg.baseline must be a scalar, a length-%d vector, or ''auto''.', cfg.N);
    end
    % A-priori check: a manual baseline below the safe margin risks flooring.
    % For a per-neuron vector, the smallest entry is the one most at risk.
    if min(baseline) < safe_base
        warning(['generate_model3_trial_data: cfg.baseline (min %.3g) is below the ' ...
            'recommended margin %.3g (= %.1f * weight_scale * M_grid) and may ' ...
            'floor negative rates. Consider cfg.baseline = ''auto'' or >= %g.'], ...
            min(baseline), safe_base, cfg.baseline_k, ceil(safe_base));
    end
end
baseline = baseline(:).';                    % 1x1 or 1xN row for broadcasting

%% Build rates, check positivity, draw Poisson trials
data_trial = cell(1, cfg.n_cue);
lambda_all = cell(1, cfg.n_cue);
n_neg = 0; n_tot = 0;
for c = 1:cfg.n_cue
    signal = Phi * W{c}.';                   % [n_stim x N] noise-free means
    lambda = signal + baseline;              % broadcast baseline over stimuli
    n_neg  = n_neg + sum(lambda(:) < 0);
    n_tot  = n_tot + numel(lambda);
    lambda = max(lambda, cfg.lambda_floor);  % safety net only (see warning)
    lambda_all{c} = lambda;

    % Vectorized Poisson draws: expand each stimulus row to n_trials rows,
    % preserving the [500 x N] blocked layout the decoding helpers expect.
    lambda_expanded = repelem(lambda, cfg.n_trials, 1);   % [n_stim*n_trials x N]
    data_trial{c} = poissrnd(lambda_expanded);
end

labels = repelem((1:n_stim).', cfg.n_trials);            % [500 x 1]

%% Positivity warning (baseline too low => rectification corrupts Model 3)
frac_floored = n_neg / n_tot;
if frac_floored > cfg.floor_warn_tol
    % Suggest a baseline that clears the most negative signal with margin.
    min_signal   = min(cellfun(@(w) min(min(Phi * w.')), W));
    rec_baseline = ceil(-min_signal + cfg.floor_margin_sd * sqrt(2) * cfg.weight_scale);
    warning(['generate_model3_trial_data: %.3f%% of rates were negative and ' ...
        'floored. This half-wave rectifies the signal and inflates the ' ...
        'participation ratio, corrupting the rank-2 linear model. Raise ' ...
        'cfg.baseline (suggest >= %g) or lower cfg.weight_scale.'], ...
        100 * frac_floored, rec_baseline);
end

%% Diagnostics (for the isolated PR / calibration test in the driver)
diagnostics = struct();
diagnostics.Phi          = Phi;                              % [n_stim x 2] shared geometry
diagnostics.W            = W;                                % 1 x n_cue, each [N x 2]
diagnostics.lambda       = lambda_all;                       % 1 x n_cue, each [n_stim x N]
diagnostics.signal_PR    = cellfun(@signal_participation_ratio, lambda_all);  % ~2 for Model 3
diagnostics.mean_count   = mean(cellfun(@(x) mean(x(:)), data_trial));
diagnostics.min_lambda   = min(cellfun(@(L) min(L(:)), lambda_all));
diagnostics.frac_floored = frac_floored;
diagnostics.cfg          = cfg;

end  % main function


%% ----------------------------------------------------------------------
%% Local functions
%% ----------------------------------------------------------------------
function cfg = default_config()
% Default Model 3 configuration. All fields overridable via the cfg input.
cfg = struct();
cfg.N              = 100;         % number of neurons (open parameter; sweep later)
cfg.n_cue          = 2;          % renderings: cue 1 = "ac", cue 2 = "ec"
cfg.n_stim         = 50;         % stimuli (fixed by the design)
cfg.n_trials       = 10;         % trials per stimulus (fixed; pipeline expects 10)
cfg.n_levels       = [10, 5];    % full-factorial grid -> 10*5 = 50 stimuli
cfg.weight_scale   = 0.7;        % SD of Gaussian readout weights (sets SNR)
cfg.rho            = 0;          % cross-cue weight correlation (0 = independent)
cfg.baseline       = 'auto';     % 'auto' = ceil(baseline_k*weight_scale*M_grid), fixed per run; or a scalar / [N x 1] vector
cfg.baseline_k     = 5.5;        % margin factor for 'auto' (~ extreme-value z* over a run + slack)
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
% mu : [n_stim x N] mean rates. Rank-2 (Model 3) -> pr ~ 2. Baseline is
% constant across stimuli and removed by the centering, so PR(lambda)=PR(signal).
mu_c = mu - mean(mu, 1);          % center across stimuli
s    = svd(mu_c, 'econ');         % singular values
ev   = s.^2;                      % covariance eigenvalues (up to scale)
pr   = sum(ev)^2 / sum(ev.^2);
end
