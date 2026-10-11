%{
% FILENAME: procrustes_decoding_simulation_sweep.m   (Simulations_v2)
% AUTHOR:   Zitong Wang
% DATE:     2026-09-17 (2026-09-23: all results in one structured file)
%
% DESCRIPTION:
%   Benchmark of Procrustes transfer decoding on synthetic populations, swept
%   over the number of neurons N and the number of stimulus classes n_stim,
%   for Model 1 (independent random means: no shared structure) and Model 3
%   (linear mixed selectivity: shared rank-2 latent geometry, independent
%   per-cue readout). One run = every (model, N, n_stim) grid point x R_pop
%   independent synthetic populations; one population = one complete synthetic
%   experiment (generation + decoding), reproducible from its seed alone, and
%   the spread across populations is the error bar.
%
%   Differences from Simulations/code/decoding_scripts/procrustes_decoding_simulation.m (v1):
%     1. CONTROL. v1's "rand null" fit the transform from a shuffled mean
%        matrix but applied it to the UNSHUFFLED test trials, which maps every
%        trial to the wrong class and sits at 1/n_stim whatever the overfitting.
%        v2 uses the self-consistent correspondence-shuffled control of the
%        real-data pipeline (neuronal_decoding_control/..._trialShuffledWithMean.m,
%        neuronal_decoding/..._incremental_pr.m and ..._incremental_pc.m): one
%        permutation shuffles BOTH the stim1 training-mean matrix and the
%        held-out stim1 test trials, a rotation-only map is fit from the
%        permuted means onto the stim2 means, applied to the identically
%        permuted test trials, and scored against the target labels. That is
%        the accuracy the rotation can manufacture with the correspondence
%        destroyed (the overfitting floor); PT is judged by its margin over it.
%        See pt_transfer_decode_population.m.
%     2. SWEEP. N and n_stim are both swept (cfg.N_list, cfg.S_list); v1 fixed
%        n_stim = 50 and ran one N per invocation.
%     3. GENERATORS. v2 generators set single-neuron statistics directly
%        (mean count and tuning depth matched to the recordings, no positivity
%        offset). See generation_scripts/.
%     4. Optional PC-dimension constraint (cfg.pc_dims = k), as in
%        procrustes_decoding_basic_incremental_pc.m. Empty = full space.
%     5. ONE structured results file per sweep (below). It is rewritten after
%        every grid point (atomically: temporary file, then replace), so it is
%        always complete up to the last finished point. Rerunning with the same
%        configuration resumes: points already marked done are skipped. A
%        different configuration under the same sweep_tag stops with an error
%        unless cfg.overwrite = true.
%
%   Quick test / override without editing this file: define a struct
%   `cfg_override` in the workspace before running, e.g.
%       cfg_override = struct('R_pop', 2, 'N_list', [10 100], 'S_list', [5 50], ...
%                             'trial_repeat', 2, 'sweep_tag', 'smoke');
%   Any cfg field can be overridden that way (the script does not `clear` it).
%
% OUTPUT: results/decoding_outputs/<sweep_tag>_results.mat, one variable `sweep`
%   sweep.info          how to read the file
%   sweep.created, .updated
%   sweep.cfg           run configuration (generator names instead of handles)
%   sweep.models        {'model1','model3'}
%   sweep.N_list        [1 x nN]  neuron counts    = dimension 1 of every grid array
%   sweep.S_list        [1 x nS]  stimulus counts  = dimension 2 of every grid array
%   sweep.metric_names  {1 x 9}   accuracy columns = last dimension of the metric arrays:
%                                 1 self  2 no_transform  3 PT_full (b*X*T+c)  4 control
%                                 5 scale_only  6 PT_rotation (X*T, matched to 4)
%                                 7 translation_only  8 chance  9 control_vs_identity (~chance)
%   sweep.pair_labels   {n_pairs x 1} ordered cue pairs, e.g. 'ac->ec' (stim1 -> stim2)
%   sweep.summary       table, one row per (model, N, n_stim): mean and SEM of every
%                       metric, PT_rotation/control ratio, PR, mean count, done flag
%   sweep.<model>       one struct per model (sweep.model1, sweep.model3):
%     .description, .generator
%     .done            [nN x nS]         grid point computed
%     .mean, .sem      [nN x nS x 9]     over populations, of the pair-averaged accuracies
%     .ratio_mean, .ratio_sem [nN x nS]  PT_rotation / control per population, then averaged
%     .ratio_of_means  [nN x nS]         mean PT_rotation / mean control
%     .PR_true         [nN x nS]         noise-free signal PR (generator), mean over populations and cues
%     .PR_emp          [nN x nS]         PR of the z-scored trial-averaged means (real-data definition)
%     .count           [nN x nS]         mean spike count
%     .pdist_mean      [nN x nS x 2]     in-sample Procrustes residual d: true fit, shuffled fit
%     .k_eff           [nN x nS]         retained PCs (NaN in the full neuron space)
%     .runtime_s       [nN x nS]         wall time of the grid point
%     .per_pop         every population (R = cfg.R_pop):
%        .seeds        [nN x nS x R]
%        .pair_metrics [nN x nS x R x n_pairs x 9]
%        .pop_metrics  [nN x nS x R x 9]           pair-averaged
%        .ratio        [nN x nS x R]               PT_rotation / control
%        .pair_pdist   [nN x nS x R x n_pairs x 2]
%        .signal_PR    [nN x nS x R x n_cue]
%        .PR_emp       [nN x nS x R x n_cue]
%        .mean_count   [nN x nS x R]
%   Example: iN = find(sweep.N_list == 100); iS = find(sweep.S_list == 50);
%            sweep.model1.mean(iN, iS, 6)                          % PT_rotation
%            squeeze(sweep.model1.per_pop.pop_metrics(iN, iS, :, 6))   % its 16 populations
%}

%% INITIALIZATION
clc; clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
addpath(script_dir);
addpath(fullfile(script_dir, '..', 'generation_scripts'));
results_root = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
if ~exist(results_root, 'dir'), mkdir(results_root); end

%% CONFIGURATION
cfg = struct();

% --- models (generator handles and descriptions in the same order) ---
cfg.models             = {'model1', 'model3'};
cfg.model_fns          = {@generate_model1_trial_data, @generate_model3_trial_data};
cfg.model_descriptions = {'independent random means (null: no structure shared across cues)', ...
                          'linear mixed selectivity (shared rank-2 latent geometry, independent per-cue readout)'};

% --- sweep axes ---
cfg.N_list = [5 10 20 30 40 50 60 70 80 90 100];   % neurons (all used, no subsampling)
cfg.S_list = [5 10 20 30 40 50];                   % stimulus classes

% --- run parameters ---
cfg.R_pop        = 16;      % populations per grid point (error bar); 16 = two full rounds of 8 workers
cfg.trial_repeat = 10;      % stim1 holdout resamples per (population, pair)
cfg.n_folds      = 10;      % k-fold CV for self-decoding (as in the real pipeline)
cfg.base_seed    = 1000;    % seeds are unique per (model, grid point, population)
cfg.seed_N_list  = [];      % reference grid that numbers the seeds ([] = N_list / S_list). Set it to
cfg.seed_S_list  = [];      %   another run's grid to reuse that run's seeds (same neurons) on a sub-grid
cfg.pc_dims      = [];      % [] = full neuron space; k = top-k mean-PCA per cue
cfg.gen          = struct();% generator overrides passed to every model, e.g. cfg.gen.cv_median = 0.2
cfg.sweep_tag    = 'sweep_N_S';   % names the file: results/decoding_outputs/<tag>_results.mat
cfg.overwrite    = false;   % false = resume (skip grid points already in the results file)

% --- optional overrides from the workspace (quick tests) ---
if exist('cfg_override', 'var') && isstruct(cfg_override)
    fo = fieldnames(cfg_override);
    for i = 1:numel(fo), cfg.(fo{i}) = cfg_override.(fo{i}); end
end

if isempty(cfg.seed_N_list), cfg.seed_N_list = cfg.N_list; end
if isempty(cfg.seed_S_list), cfg.seed_S_list = cfg.S_list; end
assert(all(ismember(cfg.N_list, cfg.seed_N_list)) && all(ismember(cfg.S_list, cfg.seed_S_list)), ...
    'Every N and n_stim of the sweep must appear in seed_N_list / seed_S_list.');

tag = cfg.sweep_tag;
if ~isempty(cfg.pc_dims), tag = sprintf('%s_pc%d', tag, cfg.pc_dims); end
cfg.tag = tag;
results_file = fullfile(results_root, sprintf('%s_results.mat', tag));

%% RESULTS FILE: resume it, or start a new one
nN = numel(cfg.N_list);  nS = numel(cfg.S_list);  n_models = numel(cfg.models);
if exist(results_file, 'file') && ~cfg.overwrite
    loaded = load(results_file, 'sweep');
    sweep  = loaded.sweep;
    check_compatible(sweep.cfg, cfg, results_file);
    n_done = sum(cellfun(@(m) nnz(sweep.(m).done), cfg.models));
    fprintf('Resuming %s (%d of %d grid points already done)\n', results_file, n_done, n_models * nN * nS);
else
    if exist(results_file, 'file')
        fprintf('cfg.overwrite = true: %s will be replaced\n', results_file);
    end
    sweep = init_sweep(cfg);
end

%% RUN
t_all = tic;
opts = struct('trial_repeat', cfg.trial_repeat, 'n_folds', cfg.n_folds, 'pc_dims', cfg.pc_dims);
fprintf('\n===== models %s : %d x %d grid, R_pop = %d, trial_repeat = %d, n_folds = %d, tag = %s =====\n', ...
    strjoin(cfg.models, ', '), nN, nS, cfg.R_pop, cfg.trial_repeat, cfg.n_folds, tag);

% grid outer, model inner: every finished grid point has all models, so a
% partial run is already comparable across models
for iS = 1:nS
    for iN = 1:nN
        N = cfg.N_list(iN);  S = cfg.S_list(iS);
        for im = 1:n_models
            model_name = cfg.models{im};
            if sweep.(model_name).done(iN, iS)
                fprintf('[%s] N=%3d S=%2d : already in the results file, skipped\n', model_name, N, S);
                continue;
            end

            % unique, reproducible seeds per (model, grid point, population),
            % numbered by the point's position in the reference seed grid
            grid_idx = (find(cfg.seed_S_list == S) - 1) * numel(cfg.seed_N_list) + find(cfg.seed_N_list == N);
            seeds    = cfg.base_seed + 100000 * (im - 1) + 1000 * (grid_idx - 1) + (1:cfg.R_pop).';

            t_pt = tic;
            R = cfg.R_pop;  gen_over = cfg.gen;  model_fn = cfg.model_fns{im};
            pop = cell(R, 1);
            parfor r = 1:R
                gen_cfg = gen_over;
                gen_cfg.N = N;  gen_cfg.n_stim = S;  gen_cfg.seed = seeds(r);
                [data_trial, labels, gdiag] = model_fn(gen_cfg);
                dec = pt_transfer_decode_population(data_trial, labels, opts);
                dec.signal_PR  = gdiag.signal_PR;
                dec.mean_count = gdiag.mean_count;
                pop{r} = dec;
            end

            % store in the structured results and rewrite the file (checkpoint)
            if isempty(sweep.metric_names)
                sweep.metric_names = pop{1}.metric_names;
                sweep.pair_labels  = pop{1}.pair_labels;
            end
            sweep.(model_name) = store_point(sweep.(model_name), iN, iS, pop, seeds, toc(t_pt), nN, nS);
            sweep = save_sweep(sweep, results_file);

            M = sweep.(model_name);
            a = squeeze(M.mean(iN, iS, :)).';
            fprintf('[%s] N=%3d S=%2d | self %.3f  noTF %.3f  PT_rot %.3f  ctrl %.3f  ratio %.2f  ctrl_id %.3f  chance %.3f | PR %.1f (emp %.1f)  count %.1f | %.0f s\n', ...
                model_name, N, S, a(1), a(2), a(6), a(4), M.ratio_mean(iN, iS), a(9), a(8), ...
                M.PR_true(iN, iS), M.PR_emp(iN, iS), M.count(iN, iS), M.runtime_s(iN, iS));
        end
    end
end
fprintf('\nAll results in %s\nTotal time: %.1f min\n', results_file, toc(t_all) / 60);


%% ======================================================================
%% Local functions
%% ======================================================================
function sweep = init_sweep(cfg)
% Empty structured results for this configuration. Grid arrays are allocated
% at each model's first finished point (their sizes come from the decoder).
nN = numel(cfg.N_list);  nS = numel(cfg.S_list);
sweep = struct();
sweep.info = sprintf([ ...
    'Procrustes transfer decoding on synthetic populations (Simulations_v2), written by\n' ...
    'procrustes_decoding_simulation_sweep.m and plotted by visualization_simulation_sweep.m.\n' ...
    'Grid arrays are [numel(N_list) x numel(S_list) x ...]: dim 1 = N (neurons), dim 2 = n_stim (classes).\n' ...
    'Accuracy columns (metric_names): 1 self, 2 no_transform, 3 PT_full, 4 control, 5 scale_only,\n' ...
    '6 PT_rotation, 7 translation_only, 8 chance, 9 control_vs_identity. Headline: PT_rotation (6)\n' ...
    'against its matched correspondence-shuffled control (4); ratio = PT_rotation / control per\n' ...
    'population, then averaged. <model>.mean and .sem are over populations of the pair-averaged\n' ...
    'accuracies; <model>.per_pop keeps every population ([nN x nS x R x ...]). summary = one table\n' ...
    'row per model x N x n_stim. done marks computed grid points (file rewritten after each one).']);
sweep.created      = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
sweep.updated      = sweep.created;
sweep.cfg          = saved_cfg(cfg);
sweep.models       = cfg.models;
sweep.N_list       = cfg.N_list;
sweep.S_list       = cfg.S_list;
sweep.metric_names = {};
sweep.pair_labels  = {};
sweep.summary      = table();
for m = 1:numel(cfg.models)
    assert(isvarname(cfg.models{m}), 'Model name ''%s'' must be a valid MATLAB field name.', cfg.models{m});
    desc = cfg.models{m};
    if isfield(cfg, 'model_descriptions') && numel(cfg.model_descriptions) >= m
        desc = cfg.model_descriptions{m};
    end
    M = struct();
    M.description    = desc;
    M.generator      = func2str(cfg.model_fns{m});
    M.done           = false(nN, nS);
    M.mean           = [];
    M.sem            = [];
    M.ratio_mean     = [];
    M.ratio_sem      = [];
    M.ratio_of_means = [];
    M.PR_true        = [];
    M.PR_emp         = [];
    M.count          = [];
    M.pdist_mean     = [];
    M.k_eff          = [];
    M.runtime_s      = [];
    M.per_pop        = struct('seeds', [], 'pair_metrics', [], 'pop_metrics', [], 'ratio', [], ...
                              'pair_pdist', [], 'signal_PR', [], 'PR_emp', [], 'mean_count', []);
    sweep.(cfg.models{m}) = M;
end
end

% ----------------------------------------------------------------------
function M = store_point(M, iN, iS, pop, seeds, runtime, nN, nS)
% Put one finished grid point (R populations) into the model's arrays.
R       = numel(pop);
n_pairs = size(pop{1}.pair_metrics, 1);
n_met   = size(pop{1}.pair_metrics, 2);
n_cue   = numel(pop{1}.PR_emp);

if isempty(M.mean)          % first finished point of this model: allocate
    M.mean           = nan(nN, nS, n_met);
    M.sem            = nan(nN, nS, n_met);
    M.ratio_mean     = nan(nN, nS);
    M.ratio_sem      = nan(nN, nS);
    M.ratio_of_means = nan(nN, nS);
    M.PR_true        = nan(nN, nS);
    M.PR_emp         = nan(nN, nS);
    M.count          = nan(nN, nS);
    M.pdist_mean     = nan(nN, nS, 2);
    M.k_eff          = nan(nN, nS);
    M.runtime_s      = nan(nN, nS);
    M.per_pop.seeds        = nan(nN, nS, R);
    M.per_pop.pair_metrics = nan(nN, nS, R, n_pairs, n_met);
    M.per_pop.pop_metrics  = nan(nN, nS, R, n_met);
    M.per_pop.ratio        = nan(nN, nS, R);
    M.per_pop.pair_pdist   = nan(nN, nS, R, n_pairs, 2);
    M.per_pop.signal_PR    = nan(nN, nS, R, n_cue);
    M.per_pop.PR_emp       = nan(nN, nS, R, n_cue);
    M.per_pop.mean_count   = nan(nN, nS, R);
end

PM  = zeros(R, n_pairs, n_met);
PD  = zeros(R, n_pairs, 2);
PRt = zeros(R, n_cue);
PRe = zeros(R, n_cue);
MC  = zeros(R, 1);
for r = 1:R
    PM(r, :, :) = pop{r}.pair_metrics;
    PD(r, :, :) = pop{r}.pair_pdist;
    PRt(r, :)   = pop{r}.signal_PR;
    PRe(r, :)   = pop{r}.PR_emp;
    MC(r)       = pop{r}.mean_count;
end
pop_metrics = reshape(mean(PM, 2), R, n_met);          % pair-averaged accuracies per population
pop_pdist   = reshape(mean(PD, 2), R, 2);
ratio       = pop_metrics(:, 6) ./ pop_metrics(:, 4);  % PT_rotation / control, per population

M.per_pop.seeds(iN, iS, :)              = seeds;
M.per_pop.pair_metrics(iN, iS, :, :, :) = PM;
M.per_pop.pop_metrics(iN, iS, :, :)     = pop_metrics;
M.per_pop.ratio(iN, iS, :)              = ratio;
M.per_pop.pair_pdist(iN, iS, :, :, :)   = PD;
M.per_pop.signal_PR(iN, iS, :, :)       = PRt;
M.per_pop.PR_emp(iN, iS, :, :)          = PRe;
M.per_pop.mean_count(iN, iS, :)         = MC;

M.mean(iN, iS, :)       = mean(pop_metrics, 1);
M.sem(iN, iS, :)        = std(pop_metrics, 0, 1) / sqrt(R);
M.ratio_mean(iN, iS)    = mean(ratio);
M.ratio_sem(iN, iS)     = std(ratio) / sqrt(R);
M.ratio_of_means(iN, iS) = mean(pop_metrics(:, 6)) / mean(pop_metrics(:, 4));
M.PR_true(iN, iS)       = mean(PRt(:));
M.PR_emp(iN, iS)        = mean(PRe(:));
M.count(iN, iS)         = mean(MC);
M.pdist_mean(iN, iS, :) = mean(pop_pdist, 1);
k = pop{1}.k_eff;  if isempty(k), k = NaN; end
M.k_eff(iN, iS)         = k;
M.runtime_s(iN, iS)     = runtime;
M.done(iN, iS)          = true;
end

% ----------------------------------------------------------------------
function sweep = save_sweep(sweep, results_file)
% Refresh the summary table and rewrite the results file atomically (write a
% temporary file, then replace), so an interruption never truncates it.
sweep.updated = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
sweep.summary = build_summary(sweep);
tmp_file = [results_file(1:end-4) '_partial.mat'];
save(tmp_file, 'sweep');
msg = '';
for attempt = 1:5                     % retry: a scanner may briefly hold the file on Windows
    [ok, msg] = movefile(tmp_file, results_file, 'f');
    if ok, return; end
    pause(2);
end
error('Could not replace %s: %s', results_file, msg);
end

% ----------------------------------------------------------------------
function T = build_summary(sweep)
% One row per (model, N, n_stim), rows ordered by model, then n_stim, then N.
T = table();
if isempty(sweep.metric_names), return; end
nN = numel(sweep.N_list);  nS = numel(sweep.S_list);
[IN, IS] = ndgrid(1:nN, 1:nS);
IN = IN(:);  IS = IS(:);
idx = sub2ind([nN nS], IN, IS);
parts = {};
for m = 1:numel(sweep.models)
    M = sweep.(sweep.models{m});
    if isempty(M.mean), continue; end          % no point of this model finished yet
    t = table();
    t.model         = repmat(string(sweep.models{m}), numel(idx), 1);
    t.N             = reshape(sweep.N_list(IN), [], 1);
    t.n_stim        = reshape(sweep.S_list(IS), [], 1);
    t.N_over_n_stim = t.N ./ t.n_stim;
    t.done          = M.done(idx);
    for k = 1:numel(sweep.metric_names)
        mk = M.mean(:, :, k);  sk = M.sem(:, :, k);
        t.(sweep.metric_names{k})          = mk(idx);
        t.([sweep.metric_names{k} '_sem']) = sk(idx);
    end
    t.ratio_PT_rotation_over_control     = M.ratio_mean(idx);
    t.ratio_PT_rotation_over_control_sem = M.ratio_sem(idx);
    t.PR_true    = M.PR_true(idx);
    t.PR_emp     = M.PR_emp(idx);
    t.mean_count = M.count(idx);
    parts{end + 1} = t; %#ok<AGROW>
end
if ~isempty(parts), T = vertcat(parts{:}); end
end

% ----------------------------------------------------------------------
function c = saved_cfg(cfg)
% cfg as stored in the file: generator names instead of function handles.
c = rmfield(cfg, 'model_fns');
c.generators = cellfun(@func2str, cfg.model_fns, 'UniformOutput', false);
end

% ----------------------------------------------------------------------
function check_compatible(saved, cfg, results_file)
% Refuse to mix grid points computed under different settings in one file.
cur  = saved_cfg(cfg);
if ~isfield(saved, 'seed_N_list'), saved.seed_N_list = saved.N_list; end   % files written before
if ~isfield(saved, 'seed_S_list'), saved.seed_S_list = saved.S_list; end   % the reference seed grid
keys = {'models', 'generators', 'N_list', 'S_list', 'seed_N_list', 'seed_S_list', 'R_pop', ...
        'trial_repeat', 'n_folds', 'base_seed', 'pc_dims', 'gen'};
for k = 1:numel(keys)
    if ~isfield(saved, keys{k}) || ~isequal(saved.(keys{k}), cur.(keys{k}))
        error('sweep:configMismatch', ['%s holds a sweep with a different ''%s''. Set ' ...
              'cfg.overwrite = true to replace it, or choose a new cfg.sweep_tag.'], results_file, keys{k});
    end
end
end
