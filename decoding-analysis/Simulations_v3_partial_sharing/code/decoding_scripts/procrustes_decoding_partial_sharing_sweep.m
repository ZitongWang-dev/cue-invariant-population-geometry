%{
% FILENAME: procrustes_decoding_partial_sharing_sweep.m   (Simulations_v3_partial_sharing)
% AUTHOR:   Zitong Wang
% DATE:     2026-10-07
%
% DESCRIPTION:
%   Procrustes transfer decoding on synthetic populations with PARTIALLY
%   shared structure across cues, swept over the shared fraction alpha, the
%   number of neurons N and the number of stimulus classes n_stim.
%   Generator: generate_partial_sharing_trial_data.m (K latent dimensions per
%   cue, K_s shared, alpha = K_s / K; every within-cue statistic is the same
%   at every alpha, only the cross-cue relation changes).
%
%   Decoding is Simulations_v2's per-population core,
%   Simulations_v2/code/decoding_scripts/pt_transfer_decode_population.m,
%   used unchanged (z-score per cue, MATLAB procrustes on trial-averaged means,
%   one held-out trial per stimulus, ECOC linear SVM, self-consistent
%   correspondence-shuffled control, 9 accuracy columns).
%
%   Seeds do not depend on alpha: at a grid point, population r has the same
%   seed at every alpha (paired design), and seeds are numbered exactly as in
%   Simulations_v2 (base_seed 1000, first model block, reference seed grid
%   N_list x [5 10 20 30 40 50]). Every population therefore has the same
%   neurons as the Model 1 population (sweep_N_S) and the fully shared Model 3
%   population (sweep_N_S_model3_D5) at the same grid point, and alpha = 1 IS
%   the sweep_N_S_model3_D5 run. That slice is imported from the v2 results
%   file (cfg.import_file) after checking that its grid, seeds, run settings
%   and generator settings match; set cfg.import_file = '' to recompute it.
%
%   Order of computation: n_stim from 50 down, N up, alpha innermost, so any
%   finished grid point has every alpha, and the 50-class row (the recorded
%   operating point) finishes first.
%
%   The results file is rewritten after every (N, n_stim, alpha) point
%   (temporary file, then replace). Rerunning the same configuration resumes;
%   a different configuration under the same sweep_tag stops with an error
%   unless cfg.overwrite = true.
%
%   Quick test without editing this file: define `cfg_override` first, e.g.
%       cfg_override = struct('R_pop', 2, 'N_list', [20 60], 'S_list', [10 50], ...
%                             'trial_repeat', 2, 'import_file', '', 'sweep_tag', 'smoke');
%
% OUTPUT: results/decoding_outputs/<sweep_tag>_results.mat, one variable `sweep`
%   sweep.info, .created, .updated, .cfg
%   sweep.N_list [1 x nN], .S_list [1 x nS], .n_shared_list [1 x nA], .alpha_list [1 x nA]
%                        grid arrays are [nN x nS x nA x ...]
%   sweep.metric_names {1 x 9}: 1 self 2 no_transform 3 PT_full 4 control 5 scale_only
%                      6 PT_rotation (matched to 4) 7 translation_only 8 chance 9 control_vs_identity
%   sweep.pair_labels  ordered cue pairs ('ac->ec', 'ec->ac')
%   sweep.done   [nN x nS x nA]  logical
%   sweep.source [nN x nS x nA]  0 = not done, 1 = computed here, 2 = imported (alpha = 1)
%   sweep.mean, .sem [nN x nS x nA x 9]   over populations, of pair-averaged accuracies
%   sweep.ratio_mean, .ratio_sem, .ratio_of_means [nN x nS x nA]   PT_rotation / control
%   sweep.headroom_mean, .headroom_sem, .headroom_of_means [nN x nS x nA]
%                        (PT_rotation - control) / (self - control): the fraction of the
%                        gap between the shuffled floor and self-decoding that PT recovers
%                        (0 = no shared structure, 1 = full transfer). _of_means uses the
%                        mean accuracies (robust where self ~ control); _mean averages the
%                        per-population values.
%   sweep.PR_true, .PR_emp, .count, .frac_cv_capped, .runtime_s [nN x nS x nA]
%   sweep.pdist_mean [nN x nS x nA x 2]   in-sample Procrustes residual: true fit, shuffled fit
%   sweep.per_pop: .seeds [nN x nS x R]; .pair_metrics [nN x nS x nA x R x n_pairs x 9];
%                  .pop_metrics [nN x nS x nA x R x 9]; .ratio, .headroom, .mean_count,
%                  .frac_cv_capped [nN x nS x nA x R]; .pair_pdist [.. x R x n_pairs x 2];
%                  .signal_PR, .PR_emp [.. x R x n_cue]; .latent_cross_sv [.. x R x K]
%                  (latent_cross_sv is NaN for the imported alpha = 1 slice: K ones by construction)
%   sweep.summary  table, one row per alpha x n_stim x N
%   Example: iN = find(sweep.N_list == 100); iS = find(sweep.S_list == 50);
%            squeeze(sweep.mean(iN, iS, :, 6)).'      % PT_rotation at every alpha
%            squeeze(sweep.headroom_of_means(iN, iS, :)).'
%}

%% INITIALIZATION
clc; clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
v3_root = fullfile(script_dir, '..', '..');
v2_root = fullfile(v3_root, '..', 'Simulations_v2');
addpath(fullfile(v3_root, 'code', 'generation_scripts'));
addpath(fullfile(v2_root, 'code', 'decoding_scripts'));          % pt_transfer_decode_population.m
results_root = fullfile(v3_root, 'results', 'decoding_outputs');
if ~exist(results_root, 'dir'), mkdir(results_root); end

%% CONFIGURATION
cfg = struct();

% --- partial-sharing axis: K latent dims per cue, K_s of them shared, alpha = K_s / K ---
cfg.n_latent      = 5;               % K (as sweep_N_S_model3_D5: matches the recorded PR)
cfg.n_shared_list = [0 1 2 3 4 5];   % K_s  ->  alpha = 0, 0.2, ..., 1

% --- sweep axes (as sweep_N_S_model3_D5) ---
cfg.N_list      = [5 10 20 30 40 50 60 70 80 90 100];   % neurons (all used, no subsampling)
cfg.S_list      = [10 20 30 40 50];                    % stimulus classes (> K)
cfg.seed_N_list = [];                                  % reference seed grid ([] = N_list)
cfg.seed_S_list = [5 10 20 30 40 50];                  % v2's grid, so seeds match sweep_N_S

% --- run parameters (as Simulations_v2) ---
cfg.R_pop        = 16;      % populations per grid point
cfg.trial_repeat = 10;      % stim1 holdout resamples per (population, pair)
cfg.n_folds      = 10;      % k-fold CV for self-decoding
cfg.base_seed    = 1000;
cfg.pc_dims      = [];      % [] = full neuron space; k = top-k mean-PCA per cue
cfg.gen          = struct('cap_rule', 'exact');   % generator overrides (n_latent / n_shared set by the driver)
cfg.sweep_tag    = 'sweep_alpha_N_S';
cfg.overwrite    = false;   % false = resume

% --- alpha = 1 slice imported from Simulations_v2 ('' = compute it here) ---
cfg.import_file  = fullfile(v2_root, 'results', 'decoding_outputs', 'sweep_N_S_model3_D5_results.mat');

% --- optional overrides from the workspace (quick tests) ---
if exist('cfg_override', 'var') && isstruct(cfg_override)
    fo = fieldnames(cfg_override);
    for i = 1:numel(fo), cfg.(fo{i}) = cfg_override.(fo{i}); end
end
if isempty(cfg.seed_N_list), cfg.seed_N_list = cfg.N_list; end
if isempty(cfg.seed_S_list), cfg.seed_S_list = cfg.S_list; end
assert(all(ismember(cfg.N_list, cfg.seed_N_list)) && all(ismember(cfg.S_list, cfg.seed_S_list)), ...
    'Every N and n_stim of the sweep must appear in seed_N_list / seed_S_list.');
assert(all(cfg.S_list > cfg.n_latent), 'Every n_stim must exceed n_latent.');
assert(all(cfg.n_shared_list >= 0 & cfg.n_shared_list <= cfg.n_latent), 'n_shared must be in 0..n_latent.');
cfg.alpha_list = cfg.n_shared_list / cfg.n_latent;

tag = cfg.sweep_tag;
if ~isempty(cfg.pc_dims), tag = sprintf('%s_pc%d', tag, cfg.pc_dims); end
cfg.tag = tag;
results_file = fullfile(results_root, sprintf('%s_results.mat', tag));

%% RESULTS FILE: resume it, or start a new one (and import alpha = 1)
nN = numel(cfg.N_list);  nS = numel(cfg.S_list);  nA = numel(cfg.n_shared_list);
if exist(results_file, 'file') && ~cfg.overwrite
    loaded = load(results_file, 'sweep');
    sweep  = loaded.sweep;
    check_compatible(sweep.cfg, cfg, results_file);
    fprintf('Resuming %s (%d of %d points already done)\n', results_file, nnz(sweep.done), nN * nS * nA);
else
    if exist(results_file, 'file')
        fprintf('cfg.overwrite = true: %s will be replaced\n', results_file);
    end
    sweep = init_sweep(cfg);
    if ~isempty(cfg.import_file)
        sweep = import_full_sharing(sweep, cfg);
        sweep = save_sweep(sweep, results_file);
    end
end

%% RUN
t_all = tic;
opts = struct('trial_repeat', cfg.trial_repeat, 'n_folds', cfg.n_folds, 'pc_dims', cfg.pc_dims);
fprintf('\n===== partial sharing: K = %d, alpha = %s : %d x %d grid, R_pop = %d, trial_repeat = %d, n_folds = %d, tag = %s =====\n', ...
    cfg.n_latent, mat2str(cfg.alpha_list, 2), nN, nS, cfg.R_pop, cfg.trial_repeat, cfg.n_folds, tag);

for iS = nS:-1:1
    for iN = 1:nN
        N = cfg.N_list(iN);  S = cfg.S_list(iS);
        seeds = point_seeds(cfg, N, S);
        for iA = 1:nA
            if sweep.done(iN, iS, iA)
                continue;
            end
            Ks = cfg.n_shared_list(iA);

            t_pt = tic;
            R = cfg.R_pop;
            gen_base = cfg.gen;
            gen_base.n_latent = cfg.n_latent;  gen_base.n_shared = Ks;
            gen_base.N = N;  gen_base.n_stim = S;
            pop = cell(R, 1);
            parfor r = 1:R
                gen_cfg = gen_base;
                gen_cfg.seed = seeds(r);
                [data_trial, labels, gdiag] = generate_partial_sharing_trial_data(gen_cfg);
                dec = pt_transfer_decode_population(data_trial, labels, opts);
                dec.signal_PR       = gdiag.signal_PR;
                dec.mean_count      = gdiag.mean_count;
                dec.frac_cv_capped  = gdiag.frac_cv_capped;
                dec.latent_cross_sv = gdiag.latent_cross_sv;
                pop{r} = dec;
            end

            if isempty(sweep.metric_names)
                sweep.metric_names = pop{1}.metric_names;
                sweep.pair_labels  = pop{1}.pair_labels;
            end
            sweep = store_point(sweep, iN, iS, iA, pop, seeds, toc(t_pt), 1);
            sweep = save_sweep(sweep, results_file);
            print_point(sweep, iN, iS, iA);
        end
    end
end
fprintf('\nAll results in %s\nTotal time: %.1f min\n', results_file, toc(t_all) / 60);


%% ======================================================================
%% Local functions
%% ======================================================================
function seeds = point_seeds(cfg, N, S)
% Same numbering as Simulations_v2 (first model block), independent of alpha.
grid_idx = (find(cfg.seed_S_list == S) - 1) * numel(cfg.seed_N_list) + find(cfg.seed_N_list == N);
seeds    = cfg.base_seed + 1000 * (grid_idx - 1) + (1:cfg.R_pop).';
end

% ----------------------------------------------------------------------
function sweep = init_sweep(cfg)
nN = numel(cfg.N_list);  nS = numel(cfg.S_list);  nA = numel(cfg.n_shared_list);
sweep = struct();
sweep.info = sprintf([ ...
    'Procrustes transfer decoding on partially shared synthetic populations (Simulations_v3_partial_sharing),\n' ...
    'written by procrustes_decoding_partial_sharing_sweep.m. Each cue has n_latent = K latent dimensions,\n' ...
    'n_shared = K_s of them shared across cues, alpha = K_s / K (within-cue statistics fixed across alpha).\n' ...
    'Grid arrays are [numel(N_list) x numel(S_list) x numel(alpha_list) x ...].\n' ...
    'Accuracy columns (metric_names): 1 self, 2 no_transform, 3 PT_full, 4 control, 5 scale_only,\n' ...
    '6 PT_rotation, 7 translation_only, 8 chance, 9 control_vs_identity. Headline: PT_rotation (6) vs its\n' ...
    'correspondence-shuffled control (4): ratio = PT/control, headroom = (PT - control)/(self - control).\n' ...
    'source: 1 = computed here, 2 = imported from Simulations_v2 sweep_N_S_model3_D5 (alpha = 1, same seeds).']);
sweep.created       = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
sweep.updated       = sweep.created;
sweep.cfg           = cfg;
sweep.N_list        = cfg.N_list;
sweep.S_list        = cfg.S_list;
sweep.n_shared_list = cfg.n_shared_list;
sweep.alpha_list    = cfg.alpha_list;
sweep.n_latent      = cfg.n_latent;
sweep.metric_names  = {};
sweep.pair_labels   = {};
sweep.summary       = table();
sweep.done          = false(nN, nS, nA);
sweep.source        = zeros(nN, nS, nA, 'uint8');
sweep.per_pop       = struct();
sweep.per_pop.seeds = nan(nN, nS, cfg.R_pop);
end

% ----------------------------------------------------------------------
function sweep = allocate(sweep, R, n_pairs, n_met, n_cue, K)
% Grid arrays, allocated at the first finished point (sizes come from the decoder).
dims = [numel(sweep.N_list), numel(sweep.S_list), numel(sweep.alpha_list)];
for f = {'ratio_mean', 'ratio_sem', 'ratio_of_means', 'headroom_mean', 'headroom_sem', ...
         'headroom_of_means', 'PR_true', 'PR_emp', 'count', 'frac_cv_capped', 'runtime_s'}
    sweep.(f{1}) = nan(dims);
end
sweep.mean       = nan([dims, n_met]);
sweep.sem        = nan([dims, n_met]);
sweep.pdist_mean = nan([dims, 2]);
P = sweep.per_pop;
P.pair_metrics    = nan([dims, R, n_pairs, n_met]);
P.pop_metrics     = nan([dims, R, n_met]);
P.ratio           = nan([dims, R]);
P.headroom        = nan([dims, R]);
P.mean_count      = nan([dims, R]);
P.frac_cv_capped  = nan([dims, R]);
P.pair_pdist      = nan([dims, R, n_pairs, 2]);
P.signal_PR       = nan([dims, R, n_cue]);
P.PR_emp          = nan([dims, R, n_cue]);
P.latent_cross_sv = nan([dims, R, K]);
sweep.per_pop = P;
end

% ----------------------------------------------------------------------
function sweep = store_point(sweep, iN, iS, iA, pop, seeds, runtime, source)
% Put one finished (N, n_stim, alpha) point (R populations) into the grid arrays.
R       = numel(pop);
n_pairs = size(pop{1}.pair_metrics, 1);
n_met   = size(pop{1}.pair_metrics, 2);
n_cue   = numel(pop{1}.PR_emp);
K       = sweep.n_latent;
if ~isfield(sweep, 'mean'), sweep = allocate(sweep, R, n_pairs, n_met, n_cue, K); end

PM = zeros(R, n_pairs, n_met);  PD = zeros(R, n_pairs, 2);
PRt = zeros(R, n_cue);  PRe = zeros(R, n_cue);
MC = nan(R, 1);  FC = nan(R, 1);  XS = nan(R, K);
for r = 1:R
    PM(r, :, :) = pop{r}.pair_metrics;
    PD(r, :, :) = pop{r}.pair_pdist;
    PRt(r, :)   = pop{r}.signal_PR;
    PRe(r, :)   = pop{r}.PR_emp;
    MC(r)       = pop{r}.mean_count;
    if isfield(pop{r}, 'frac_cv_capped'),  FC(r)    = pop{r}.frac_cv_capped;  end
    if isfield(pop{r}, 'latent_cross_sv'), XS(r, :) = pop{r}.latent_cross_sv; end
end
pm       = reshape(mean(PM, 2), R, n_met);               % pair-averaged accuracies per population
pdist    = reshape(mean(PD, 2), R, 2);
ratio    = pm(:, 6) ./ pm(:, 4);
headroom = (pm(:, 6) - pm(:, 4)) ./ (pm(:, 1) - pm(:, 4));
mu       = mean(pm, 1);

P = sweep.per_pop;
P.seeds(iN, iS, :)                  = seeds;
P.pair_metrics(iN, iS, iA, :, :, :) = PM;
P.pop_metrics(iN, iS, iA, :, :)     = pm;
P.ratio(iN, iS, iA, :)              = ratio;
P.headroom(iN, iS, iA, :)           = headroom;
P.mean_count(iN, iS, iA, :)         = MC;
P.frac_cv_capped(iN, iS, iA, :)     = FC;
P.pair_pdist(iN, iS, iA, :, :, :)   = PD;
P.signal_PR(iN, iS, iA, :, :)       = PRt;
P.PR_emp(iN, iS, iA, :, :)          = PRe;
P.latent_cross_sv(iN, iS, iA, :, :) = XS;
sweep.per_pop = P;

sweep.mean(iN, iS, iA, :)         = mu;
sweep.sem(iN, iS, iA, :)          = std(pm, 0, 1) / sqrt(R);
sweep.ratio_mean(iN, iS, iA)      = mean(ratio);
sweep.ratio_sem(iN, iS, iA)       = std(ratio) / sqrt(R);
sweep.ratio_of_means(iN, iS, iA)  = mu(6) / mu(4);
sweep.headroom_mean(iN, iS, iA)   = mean(headroom, 'omitnan');
sweep.headroom_sem(iN, iS, iA)    = std(headroom, 0, 'omitnan') / sqrt(sum(~isnan(headroom)));
sweep.headroom_of_means(iN, iS, iA) = (mu(6) - mu(4)) / (mu(1) - mu(4));
sweep.PR_true(iN, iS, iA)         = mean(PRt(:));
sweep.PR_emp(iN, iS, iA)          = mean(PRe(:));
sweep.count(iN, iS, iA)           = mean(MC);
sweep.frac_cv_capped(iN, iS, iA)  = mean(FC);
sweep.pdist_mean(iN, iS, iA, :)   = mean(pdist, 1);
sweep.runtime_s(iN, iS, iA)       = runtime;
sweep.done(iN, iS, iA)            = true;
sweep.source(iN, iS, iA)          = source;
end

% ----------------------------------------------------------------------
function sweep = import_full_sharing(sweep, cfg)
% Copy the alpha = 1 slice from Simulations_v2 sweep_N_S_model3_D5 after checking
% that it was computed with this grid, these seeds and these settings.
iA = find(cfg.n_shared_list == cfg.n_latent);
if isempty(iA), return; end
assert(exist(cfg.import_file, 'file') == 2, 'Import file not found: %s', cfg.import_file);
src = load(cfg.import_file, 'sweep');  src = src.sweep;
c   = src.cfg;
M   = src.model3;
need = {'R_pop', 'trial_repeat', 'n_folds', 'pc_dims'};
for k = 1:numel(need)
    assert(isequal(c.(need{k}), cfg.(need{k})), 'Import: ''%s'' differs from this run.', need{k});
end
assert(isequal(c.gen.n_latent, cfg.n_latent) && strcmpi(c.gen.cap_rule, cfg.gen.cap_rule), ...
    'Import: generator settings (n_latent, cap_rule) differ from this run.');
extra = setdiff(fieldnames(cfg.gen), {'cap_rule'});
assert(isempty(extra), 'Import: this run overrides generator fields the import did not use (%s).', strjoin(extra, ', '));

n_imp = 0;
for iS = 1:numel(cfg.S_list)
    for iN = 1:numel(cfg.N_list)
        jN = find(src.N_list == cfg.N_list(iN));  jS = find(src.S_list == cfg.S_list(iS));
        if isempty(jN) || isempty(jS) || ~M.done(jN, jS), continue; end
        seeds = point_seeds(cfg, cfg.N_list(iN), cfg.S_list(iS));
        assert(isequal(squeeze(M.per_pop.seeds(jN, jS, :)), seeds), ...
            'Import: seeds differ at N = %d, n_stim = %d.', cfg.N_list(iN), cfg.S_list(iS));
        R = cfg.R_pop;
        pop = cell(R, 1);
        for r = 1:R
            p = struct();
            p.pair_metrics = squeeze(M.per_pop.pair_metrics(jN, jS, r, :, :));
            p.pair_pdist   = squeeze(M.per_pop.pair_pdist(jN, jS, r, :, :));
            p.signal_PR    = squeeze(M.per_pop.signal_PR(jN, jS, r, :)).';
            p.PR_emp       = squeeze(M.per_pop.PR_emp(jN, jS, r, :)).';
            p.mean_count   = M.per_pop.mean_count(jN, jS, r);
            pop{r} = p;
        end
        if isempty(sweep.metric_names)
            sweep.metric_names = src.metric_names;
            sweep.pair_labels  = src.pair_labels;
        end
        sweep = store_point(sweep, iN, iS, iA, pop, seeds, M.runtime_s(jN, jS), 2);
        n_imp = n_imp + 1;
    end
end
fprintf('Imported %d alpha = 1 points from %s\n', n_imp, cfg.import_file);
end

% ----------------------------------------------------------------------
function print_point(sweep, iN, iS, iA)
a = squeeze(sweep.mean(iN, iS, iA, :)).';
fprintf(['[alpha %.1f] N=%3d S=%2d | self %.3f  noTF %.3f  PT_rot %.3f  ctrl %.3f  ratio %.2f  headroom %.2f' ...
         '  ctrl_id %.3f  chance %.3f | PR %.1f (emp %.1f)  count %.1f  capped %.2f | %.0f s\n'], ...
    sweep.alpha_list(iA), sweep.N_list(iN), sweep.S_list(iS), a(1), a(2), a(6), a(4), ...
    sweep.ratio_mean(iN, iS, iA), sweep.headroom_of_means(iN, iS, iA), a(9), a(8), ...
    sweep.PR_true(iN, iS, iA), sweep.PR_emp(iN, iS, iA), sweep.count(iN, iS, iA), ...
    sweep.frac_cv_capped(iN, iS, iA), sweep.runtime_s(iN, iS, iA));
end

% ----------------------------------------------------------------------
function sweep = save_sweep(sweep, results_file)
% Refresh the summary table and rewrite the results file atomically.
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
% One row per (alpha, n_stim, N), ordered by alpha, then n_stim, then N.
T = table();
if ~isfield(sweep, 'mean'), return; end
nN = numel(sweep.N_list);  nS = numel(sweep.S_list);  nA = numel(sweep.alpha_list);
[IN, IS, IA] = ndgrid(1:nN, 1:nS, 1:nA);
idx = sub2ind([nN nS nA], IN(:), IS(:), IA(:));
T.alpha         = reshape(sweep.alpha_list(IA(:)), [], 1);
T.n_shared      = reshape(sweep.n_shared_list(IA(:)), [], 1);
T.n_stim        = reshape(sweep.S_list(IS(:)), [], 1);
T.N             = reshape(sweep.N_list(IN(:)), [], 1);
T.N_over_n_stim = T.N ./ T.n_stim;
T.done          = pick(sweep.done, idx);
T.source        = pick(sweep.source, idx);
for k = 1:numel(sweep.metric_names)
    T.(sweep.metric_names{k})          = pick(sweep.mean(:, :, :, k), idx);
    T.([sweep.metric_names{k} '_sem']) = pick(sweep.sem(:, :, :, k), idx);
end
T.ratio_PT_rotation_over_control     = pick(sweep.ratio_mean, idx);
T.ratio_PT_rotation_over_control_sem = pick(sweep.ratio_sem, idx);
T.ratio_of_means                     = pick(sweep.ratio_of_means, idx);
T.headroom_of_means                  = pick(sweep.headroom_of_means, idx);
T.headroom_mean                      = pick(sweep.headroom_mean, idx);
T.headroom_sem                       = pick(sweep.headroom_sem, idx);
T.PR_true        = pick(sweep.PR_true, idx);
T.PR_emp         = pick(sweep.PR_emp, idx);
T.mean_count     = pick(sweep.count, idx);
T.frac_cv_capped = pick(sweep.frac_cv_capped, idx);
T = sortrows(T, {'alpha', 'n_stim', 'N'});
end

function v = pick(A, idx)
% Column of A's entries at linear grid indices idx (any singleton dimensions).
v = A(:);
v = v(idx);
end

% ----------------------------------------------------------------------
function check_compatible(saved, cfg, results_file)
% Refuse to mix points computed under different settings in one file.
keys = {'n_latent', 'n_shared_list', 'N_list', 'S_list', 'seed_N_list', 'seed_S_list', 'R_pop', ...
        'trial_repeat', 'n_folds', 'base_seed', 'pc_dims', 'gen', 'import_file'};
for k = 1:numel(keys)
    if ~isfield(saved, keys{k}) || ~isequal(saved.(keys{k}), cfg.(keys{k}))
        error('sweep:configMismatch', ['%s holds a sweep with a different ''%s''. Set ' ...
              'cfg.overwrite = true to replace it, or choose a new cfg.sweep_tag.'], results_file, keys{k});
    end
end
end
