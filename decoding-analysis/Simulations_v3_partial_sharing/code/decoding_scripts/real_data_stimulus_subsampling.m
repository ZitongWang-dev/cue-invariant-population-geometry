%{
% FILENAME: real_data_stimulus_subsampling.m   (Simulations_v3_partial_sharing)
% AUTHOR:   Zitong Wang
% DATE:     2026-10-08
%
% DESCRIPTION:
%   Procrustes transfer decoding of the RECORDED populations on random subsets
%   of C of the 50 stimuli, so the recordings can be compared with the
%   synthetic sweeps along the number-of-classes axis. In the partially shared
%   family (generate_partial_sharing_trial_data.m) the headroom
%       (PT_rotation - control) / (self - control)
%   depends on the shared fraction alpha but hardly on C or N, while
%   PT/control rises steeply with C. If that family describes the recordings,
%   their headroom should stay roughly flat as C drops from 50.
%
%   Decoding is exactly the synthetic pipeline: each subset is passed to
%   Simulations_v2/code/decoding_scripts/pt_transfer_decode_population.m
%   (unchanged), which z-scores each cue per neuron over the subset's trials,
%   fits one ECOC linear SVM per target cue (self-decoding by 10-fold CV),
%   fits MATLAB procrustes on trial-averaged means with one held-out trial per
%   stimulus (trial_repeat holdouts), and scores the self-consistent
%   correspondence-shuffled control. All 6 ordered cue pairs (ac, ec, ex).
%   (The recordings' own incremental_pr script z-scores over all 50 stimuli
%   and refits the SVM for every holdout; here the subset is treated like a
%   synthetic population of C stimuli.)
%
%   Two modes:
%     'full_N'          all recorded neurons, C random stimuli
%     'matched_ratio'   N = round(N_full * C / 50) random neurons (redrawn per
%                       subset), so N / C stays at the population's recorded ratio
%   C = 50 uses every stimulus and every neuron in both modes, so it is
%   computed once (n_repeat_full repeats of the trial holdouts and CV folds)
%   and stored in both. C < 50: n_subsets random stimulus subsets.
%   A neuron with no spikes in any trial of a cue's selected stimuli is kept:
%   MATLAB's zscore maps a constant column to zeros, as in the recordings'
%   own pipeline (FR V2 has 7 such neurons over all 50 stimuli). The number of
%   such silent neurons per subset (in any cue) is stored as n_silent.
%
%   Every subset is reproducible from its seed (stimulus draw, neuron draw,
%   holdouts, shuffles, CV folds). Results go to ONE structured file,
%   rewritten after every (population, mode, C) cell; rerunning resumes.
%
% INPUT:  neuronal_decoding/neuronal_data/<M>/<A>/<M>_<A>_allstim.mat ('three_stim_array'),
%         spike counts in [330 630] ms (same loader logic as multiclass_svmloader_PT.m)
% OUTPUT: Simulations_v3_partial_sharing/results/decoding_outputs/real_stimsub_results.mat,
%         one variable `real`:
%   real.populations {1 x nP}  'FR_V1' 'FR_V2' 'KO_V1' 'KO_V2';  real.N_full [1 x nP]
%   real.modes {1 x nM}        'full_N' 'matched_ratio'
%   real.C_list [1 x nC]       stimulus counts;  real.n_rep [1 x nC] subsets (repeats) per C
%   real.metric_names, real.pair_labels   (as the synthetic sweeps; 6 ordered pairs)
%   real.done [nP x nM x nC]
%   real.mean, .sem [nP x nM x nC x 9]       over subsets, of pair-averaged accuracies
%   real.pair_mean  [nP x nM x nC x 6 x 9]   per ordered pair, mean over subsets
%   real.ratio_of_means, .headroom_of_means, .headroom_mean, .headroom_sem, .PR_emp,
%        .N_used [nP x nM x nC]
%   real.per_sub: .pair_metrics [nP x nM x nC x R x 6 x 9]; .pop_metrics [.. x R x 9];
%        .headroom, .ratio, .n_silent, .N_used, .seed [.. x R]; .PR_emp [.. x R x 3];
%        .stim_idx [.. x R x 50]; .neuron_idx [.. x R x max N_full]   (NaN-padded)
%   real.summary  table, one row per population x mode x C
%}

%% INITIALIZATION
clc; clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
v3_root = fullfile(script_dir, '..', '..');
da_root = fullfile(v3_root, '..');
addpath(fullfile(da_root, 'Simulations_v2', 'code', 'decoding_scripts'));   % pt_transfer_decode_population.m
results_root = fullfile(v3_root, 'results', 'decoding_outputs');
if ~exist(results_root, 'dir'), mkdir(results_root); end

%% CONFIGURATION
cfg = struct();
cfg.populations   = {'FR_V1', 'FR_V2', 'KO_V1', 'KO_V2'};
cfg.data_root     = fullfile(da_root, 'neuronal_decoding', 'neuronal_data');
cfg.timewindow    = [330 630];          % ms, as every real-data decoding script
cfg.n_trials      = 10;                 % trials per stimulus (recorded)
cfg.n_stim_total  = 50;
cfg.C_list        = [10 20 30 40 50];   % stimulus classes (as the synthetic sweeps)
cfg.modes         = {'full_N', 'matched_ratio'};
cfg.n_subsets     = 16;                 % random stimulus subsets per C < 50 (two rounds of 8 workers)
cfg.n_repeat_full = 8;                  % repeats at C = 50 (holdouts / CV folds only; one round)
cfg.trial_repeat  = 10;                 % stim1 holdout resamples per (subset, pair)
cfg.n_folds       = 10;                 % k-fold CV for self-decoding
cfg.base_seed     = 500000;
cfg.tag           = 'real_stimsub';
cfg.overwrite     = false;

if exist('cfg_override', 'var') && isstruct(cfg_override)
    fo = fieldnames(cfg_override);
    for i = 1:numel(fo), cfg.(fo{i}) = cfg_override.(fo{i}); end
end
results_file = fullfile(results_root, sprintf('%s_results.mat', cfg.tag));

%% LOAD the recorded spike counts (1 x 3 cell per population: ac, ec, ex)
nP = numel(cfg.populations);
counts = cell(1, nP);
N_full = zeros(1, nP);
for p = 1:nP
    parts = strsplit(cfg.populations{p}, '_');
    f = fullfile(cfg.data_root, parts{1}, parts{2}, sprintf('%s_%s_allstim.mat', parts{1}, parts{2}));
    tmp = load(f, 'three_stim_array');
    counts{p} = cellfun(@(a) spike_counts(a, cfg.timewindow, cfg.n_trials), tmp.three_stim_array, ...
                        'UniformOutput', false);
    N_full(p) = size(counts{p}{1}, 2);
    assert(all(cellfun(@(x) isequal(size(x), [cfg.n_stim_total * cfg.n_trials, N_full(p)]), counts{p})), ...
        'Unexpected data size for %s.', cfg.populations{p});
    fprintf('%s: %d neurons, %d cues\n', cfg.populations{p}, N_full(p), numel(counts{p}));
end
cfg.N_full = N_full;

%% RESULTS FILE: resume or start
nM = numel(cfg.modes);  nC = numel(cfg.C_list);
if exist(results_file, 'file') && ~cfg.overwrite
    loaded = load(results_file, 'real');
    real   = loaded.real;
    check_compatible(real.cfg, cfg, results_file);
    fprintf('Resuming %s (%d of %d cells done)\n', results_file, nnz(real.done), numel(real.done));
else
    real = init_real(cfg);
end

%% RUN: C from 50 down, then population, then mode
opts = struct('trial_repeat', cfg.trial_repeat, 'n_folds', cfg.n_folds, 'pc_dims', []);
t_all = tic;
for iC = nC:-1:1
    C = cfg.C_list(iC);
    for p = 1:nP
        for m = 1:nM
            if real.done(p, m, iC), continue; end
            t_cell = tic;
            if C == cfg.n_stim_total && m > 1           % identical procedure: copy mode 1
                assert(real.done(p, 1, iC), 'C = 50 must be computed in mode 1 first.');
                real = copy_cell(real, p, 1, m, iC);
            else
                if C == cfg.n_stim_total, R = cfg.n_repeat_full; else, R = cfg.n_subsets; end
                if strcmp(cfg.modes{m}, 'matched_ratio')
                    N_use = round(N_full(p) * C / cfg.n_stim_total);
                else
                    N_use = N_full(p);
                end
                seeds = cfg.base_seed + 100000 * (p - 1) + 10000 * (m - 1) + 100 * (iC - 1) + (1:R).';
                X_all = counts{p};  n_tr = cfg.n_trials;  S_tot = cfg.n_stim_total;  Nf = N_full(p);
                sub = cell(R, 1);
                parfor r = 1:R
                    rng(seeds(r), 'twister');
                    stim_idx   = sort(randperm(S_tot, C));
                    neuron_idx = sort(randperm(Nf, N_use));
                    rows = reshape((stim_idx - 1) * n_tr + (1:n_tr).', [], 1);   % blocked by stimulus
                    X = cellfun(@(d) d(rows, neuron_idx), X_all, 'UniformOutput', false);
                    silent = false(1, numel(neuron_idx));
                    for c = 1:numel(X)
                        silent = silent | var(X{c}, 0, 1) == 0;
                    end
                    labels = repelem((1:C).', n_tr);
                    out = pt_transfer_decode_population(X, labels, opts);
                    out.stim_idx   = stim_idx;
                    out.neuron_idx = neuron_idx;
                    out.n_silent   = sum(silent);
                    out.seed       = seeds(r);
                    sub{r} = out;
                end
                if isempty(real.metric_names)
                    real.metric_names = sub{1}.metric_names;
                    real.pair_labels  = sub{1}.pair_labels;
                end
                real = store_cell(real, p, m, iC, sub, toc(t_cell));
            end
            real = save_real(real, results_file);
            a = squeeze(real.mean(p, m, iC, :)).';
            fprintf(['[%s %-13s] C=%2d N=%3d | self %.3f  noTF %.3f  PT_rot %.3f  ctrl %.3f  ratio %.2f  headroom %.2f' ...
                     '  ctrl_id %.3f  chance %.3f | PR %.1f  silent %.1f | %.0f s\n'], ...
                cfg.populations{p}, cfg.modes{m}, C, real.N_used(p, m, iC), a(1), a(2), a(6), a(4), ...
                real.ratio_of_means(p, m, iC), real.headroom_of_means(p, m, iC), a(9), a(8), ...
                real.PR_emp(p, m, iC), mean(real.per_sub.n_silent(p, m, iC, :), 'omitnan'), toc(t_cell));
        end
    end
end
fprintf('\nAll results in %s\nTotal time: %.1f min\n', results_file, toc(t_all) / 60);


%% ======================================================================
%% Local functions
%% ======================================================================
function data = spike_counts(array, timewindow, repeat)
% Spike counts in timewindow, [condition*repeat x neuron], blocked by condition.
% Same logic as neuronal_decoding/code/decoding_scripts/multiclass_svmloader_PT.m.
[neuron_num, condition_num] = size(array);
data = zeros(condition_num * repeat, neuron_num);
for con = 1:condition_num
    raster = cell2mat(array(:, con));                    % [neuron*repeat x time]
    sc = sum(raster(:, timewindow(1):timewindow(2)), 2);
    data((con - 1) * repeat + (1:repeat), :) = reshape(sc, repeat, neuron_num);
end
end

% ----------------------------------------------------------------------
function real = init_real(cfg)
nP = numel(cfg.populations);  nM = numel(cfg.modes);  nC = numel(cfg.C_list);
R  = max(cfg.n_subsets, cfg.n_repeat_full);
real = struct();
real.info = sprintf([ ...
    'Procrustes transfer decoding of the recorded populations on random subsets of C stimuli\n' ...
    '(Simulations_v3_partial_sharing/code/decoding_scripts/real_data_stimulus_subsampling.m), decoded with\n' ...
    'Simulations_v2 pt_transfer_decode_population.m exactly as the synthetic sweeps. Arrays are\n' ...
    '[population x mode x C x ...]. Modes: full_N (all neurons), matched_ratio (N = N_full*C/50).\n' ...
    'Accuracy columns: 1 self, 2 no_transform, 3 PT_full, 4 control, 5 scale_only, 6 PT_rotation,\n' ...
    '7 translation_only, 8 chance, 9 control_vs_identity. headroom = (PT_rotation - control)/(self - control).\n' ...
    'C = 50 is computed once (all stimuli, all neurons) and stored in both modes.']);
real.created      = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
real.updated      = real.created;
real.cfg          = cfg;
real.populations  = cfg.populations;
real.N_full       = cfg.N_full;
real.modes        = cfg.modes;
real.C_list       = cfg.C_list;
real.n_rep        = cfg.n_subsets * ones(1, nC);
real.n_rep(cfg.C_list == cfg.n_stim_total) = cfg.n_repeat_full;
real.metric_names = {};
real.pair_labels  = {};
real.done         = false(nP, nM, nC);
d = [nP, nM, nC];
real.mean = nan([d 9]);  real.sem = nan([d 9]);  real.pair_mean = nan([d 6 9]);
for f = {'ratio_of_means', 'headroom_of_means', 'headroom_mean', 'headroom_sem', 'PR_emp', 'N_used', 'runtime_s'}
    real.(f{1}) = nan(d);
end
P = struct();
P.pair_metrics = nan([d R 6 9]);
P.pop_metrics  = nan([d R 9]);
P.headroom     = nan([d R]);
P.ratio        = nan([d R]);
P.n_silent    = nan([d R]);
P.N_used       = nan([d R]);
P.seed         = nan([d R]);
P.PR_emp       = nan([d R 3]);
P.stim_idx     = nan([d R cfg.n_stim_total]);
P.neuron_idx   = nan([d R max(cfg.N_full)]);
real.per_sub   = P;
real.summary   = table();
end

% ----------------------------------------------------------------------
function real = store_cell(real, p, m, iC, sub, runtime)
R = numel(sub);
PM = zeros(R, 6, 9);  PR = zeros(R, 3);
P = real.per_sub;
for r = 1:R
    PM(r, :, :) = sub{r}.pair_metrics;
    PR(r, :)    = sub{r}.PR_emp;
    P.n_silent(p, m, iC, r) = sub{r}.n_silent;
    P.N_used(p, m, iC, r)    = numel(sub{r}.neuron_idx);
    P.seed(p, m, iC, r)      = sub{r}.seed;
    P.stim_idx(p, m, iC, r, 1:numel(sub{r}.stim_idx))     = sub{r}.stim_idx;
    P.neuron_idx(p, m, iC, r, 1:numel(sub{r}.neuron_idx)) = sub{r}.neuron_idx;
end
pm = reshape(mean(PM, 2), R, 9);                  % pair-averaged accuracies per subset
hr = (pm(:, 6) - pm(:, 4)) ./ (pm(:, 1) - pm(:, 4));
P.pair_metrics(p, m, iC, 1:R, :, :) = PM;
P.pop_metrics(p, m, iC, 1:R, :)     = pm;
P.headroom(p, m, iC, 1:R)           = hr;
P.ratio(p, m, iC, 1:R)              = pm(:, 6) ./ pm(:, 4);
P.PR_emp(p, m, iC, 1:R, :)          = PR;
real.per_sub = P;
mu = mean(pm, 1);
real.mean(p, m, iC, :)          = mu;
real.sem(p, m, iC, :)           = std(pm, 0, 1) / sqrt(R);
real.pair_mean(p, m, iC, :, :)  = mean(PM, 1);
real.ratio_of_means(p, m, iC)   = mu(6) / mu(4);
real.headroom_of_means(p, m, iC) = (mu(6) - mu(4)) / (mu(1) - mu(4));
real.headroom_mean(p, m, iC)    = mean(hr, 'omitnan');
real.headroom_sem(p, m, iC)     = std(hr, 0, 'omitnan') / sqrt(sum(~isnan(hr)));
real.PR_emp(p, m, iC)           = mean(PR(:));
real.N_used(p, m, iC)           = mean(P.N_used(p, m, iC, 1:R));
real.runtime_s(p, m, iC)        = runtime;
real.done(p, m, iC)             = true;
end

% ----------------------------------------------------------------------
function real = copy_cell(real, p, m_from, m_to, iC)
% C = 50: both modes run the identical procedure, so mode m_to reuses m_from.
for f = {'mean', 'sem', 'pair_mean'}
    real.(f{1})(p, m_to, iC, :) = real.(f{1})(p, m_from, iC, :);
end
for f = {'ratio_of_means', 'headroom_of_means', 'headroom_mean', 'headroom_sem', 'PR_emp', 'N_used', 'runtime_s'}
    real.(f{1})(p, m_to, iC) = real.(f{1})(p, m_from, iC);
end
fs = fieldnames(real.per_sub);
for k = 1:numel(fs)
    real.per_sub.(fs{k})(p, m_to, iC, :) = real.per_sub.(fs{k})(p, m_from, iC, :);
end
real.done(p, m_to, iC) = true;
end

% ----------------------------------------------------------------------
function real = save_real(real, results_file)
real.updated = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
real.summary = build_summary(real);
tmp_file = [results_file(1:end-4) '_partial.mat'];
save(tmp_file, 'real');
msg = '';
for attempt = 1:5
    [ok, msg] = movefile(tmp_file, results_file, 'f');
    if ok, return; end
    pause(2);
end
error('Could not replace %s: %s', results_file, msg);
end

% ----------------------------------------------------------------------
function T = build_summary(real)
T = table();
if isempty(real.metric_names), return; end
nP = numel(real.populations);  nM = numel(real.modes);  nC = numel(real.C_list);
[IP, IM, IC] = ndgrid(1:nP, 1:nM, 1:nC);
idx = sub2ind([nP nM nC], IP(:), IM(:), IC(:));
T.population = reshape(string(real.populations(IP(:))), [], 1);
T.mode       = reshape(string(real.modes(IM(:))), [], 1);
T.C          = reshape(real.C_list(IC(:)), [], 1);
T.N_used     = pick(real.N_used, idx);
T.done       = pick(real.done, idx);
for k = 1:numel(real.metric_names)
    T.(real.metric_names{k})          = pick(real.mean(:, :, :, k), idx);
    T.([real.metric_names{k} '_sem']) = pick(real.sem(:, :, :, k), idx);
end
T.ratio_of_means    = pick(real.ratio_of_means, idx);
T.headroom_of_means = pick(real.headroom_of_means, idx);
T.headroom_mean     = pick(real.headroom_mean, idx);
T.headroom_sem      = pick(real.headroom_sem, idx);
T.PR_emp            = pick(real.PR_emp, idx);
T = sortrows(T, {'population', 'mode', 'C'});
end

function v = pick(A, idx)
v = A(:);
v = v(idx);
end

% ----------------------------------------------------------------------
function check_compatible(saved, cfg, results_file)
keys = {'populations', 'timewindow', 'C_list', 'modes', 'n_subsets', 'n_repeat_full', ...
        'trial_repeat', 'n_folds', 'base_seed', 'N_full'};
for k = 1:numel(keys)
    if ~isfield(saved, keys{k}) || ~isequal(saved.(keys{k}), cfg.(keys{k}))
        error('real:configMismatch', ['%s holds a run with a different ''%s''. Set ' ...
              'cfg.overwrite = true to replace it, or choose a new cfg.tag.'], results_file, keys{k});
    end
end
end
