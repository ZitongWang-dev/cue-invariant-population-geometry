%{
Filename: run_gabor_pc_sweep.m
Author:   Zitong Wang
Date:     2026-10-15

Description:
    Cue-transfer decoding of the Gabor-filter cell models as a function of the
    number of retained principal components k, on whole models (default: both,
    simple 192 cells and complex 48 cells), six ordered rendering pairs.

    The PC projection is the one of the neuron and CNN PC sweeps
    (neuron_decoding_rev/run_pc_sweep.m, ANN_decoding_rev/run_ann_pc_sweep.m),
    done inside the shared core pt_transfer_core.m: each side is projected onto
    the top-k principal directions of ITS OWN centered class-mean matrix
    (training trials only), and the rotation is a k x k map with k(k-1)/2 free
    parameters whatever the model size, so model cells, CNN units and neurons
    can be compared at the same rotation freedom (k = 5 in the main
    comparison). Non-PT projects the source trials onto the TARGET basis.

    k grid: cfg.k_base (default 1:10, 15:5:45, 49), with values above the
    model's maximum rank removed and that maximum added (simple: 49 =
    n_stim - 1; complex: 48 cells). The fold partition and the control
    permutations depend on the pair but not on k, so the k curves share their
    random numbers.

Run:
    By default it runs both models. To run one:
        cfg_override = struct('models', {{'simple'}});
        run_gabor_pc_sweep
    Re-running the same configuration resumes (a finished model is skipped).
    Clear cfg_override before running another driver.

Output: results/decoding_outputs/gabor_pc_sweep/<model>.mat, variable `res`:
    res.model, res.population, res.N, res.n_silent
    res.k_list [nK x 1]; res.pair_labels, res.pair_group {6 x 1}
    [nK x 6 pairs x 10 folds]: self, nopt, pt_rot, pt_full, d_rot, d_full
    [nK x 6 x 10 x n_perm]: ctrl, ctrl_id
    res.pr_raw, res.pr_cv [1 x 3 cues]; res.cumvar, res.cv_cumvar {1 x 3}
        (variance captured by the top k PCs of cue c = res.cumvar{c}(k))
    res.summary: one row per k (pooled and per pair) with ratios of pooled
    means (transfer_metrics.m).
%}

%% Initialization
clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
core_dir = fullfile(script_dir, '..', '..', '..', 'neuron_decoding_rev', 'code', 'decoding_scripts');
assert(exist(fullfile(core_dir, 'pt_transfer_core.m'), 'file') == 2, 'Shared core not found in %s', core_dir);
addpath(script_dir);
addpath(core_dir);                  % pt_transfer_core, population_pr, transfer_metrics, progress_queue

%% Configuration
cfg = struct();
cfg.models    = {'simple', 'complex'};   % default: both
cfg.k_base    = [1:10, 15:5:45, 49];
cfg.n_perm    = 10;
cfg.pr_splits = 20;
cfg.pairs     = {'ac','ec'; 'ec','ac'; 'ec','ex'; 'ex','ec'; 'ac','ex'; 'ex','ac'};   % {source, target}
cfg.base_seed = 33000;
cfg.data_root = '';                 % '' = Gabor_decoding/.../sampled_trials (load_gabor_model default)
cfg.out_root  = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag       = 'gabor_pc_sweep';
cfg.overwrite = false;

if exist('cfg_override', 'var') && isstruct(cfg_override)
    fo  = fieldnames(cfg_override);
    bad = setdiff(fo, fieldnames(cfg));
    assert(isempty(bad), 'cfg_override has unknown field(s): %s', strjoin(bad', ', '));
    for i = 1:numel(fo), cfg.(fo{i}) = cfg_override.(fo{i}); end
    fprintf('cfg_override applied: %s\n', strjoin(fo', ', '));
end

out_dir = fullfile(cfg.out_root, cfg.tag);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
if isempty(gcp('nocreate')), parpool('Processes'); end
t_all = tic;

%% Run
nP = size(cfg.pairs, 1);
for im = 1:numel(cfg.models)
    model    = lower(cfg.models{im});
    model_id = gabor_model_id(model);
    pop = load_gabor_model(model, struct('data_root', cfg.data_root));
    F = pop.n_rep;  M = cfg.n_perm;
    k_list = k_grid(cfg.k_base, min([pop.n_stim - 1, pop.N]));
    nK = numel(k_list);
    src = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 1));
    tgt = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 2));

    file = fullfile(out_dir, sprintf('%s.mat', model));
    [~, data_name] = fileparts(pop.file);
    sig  = struct('k_list', k_list, 'n_perm', M, 'pr_splits', cfg.pr_splits, 'pairs', {cfg.pairs}, ...
        'base_seed', cfg.base_seed, 'data_file', data_name);
    res = open_results(file, sig, cfg.overwrite);
    if ~isempty(res) && all(res.done)
        fprintf('%s PC sweep already done, skipped.\n', pop.name);
        continue;
    end
    res = struct();
    res.info        = 'Gabor-model PC sweep, cue transfer, whole model. Field definitions: header of run_gabor_pc_sweep.m.';
    res.population  = pop.name;  res.model = model;  res.N = pop.N;
    res.n_silent    = pop.n_silent;
    res.sig         = sig;  res.cfg = cfg;  res.created = char(datetime('now'));
    res.k_list      = k_list(:);
    res.pair_labels = strcat(cfg.pairs(:, 1), '->', cfg.pairs(:, 2));
    res.pair_group  = pair_group_names(cfg.pairs);
    res.cue_names   = pop.cue_names;
    res.done        = false;
    for c = 1:3
        pr = population_pr(zscore(pop.data{c}), pop.labels, ...
            struct('n_split', cfg.pr_splits, 'seed', cfg.base_seed + model_id * 1e7 + 10 + c));
        res.pr_raw(1, c) = pr.raw;  res.pr_cv(1, c) = pr.cv;
        res.cumvar{1, c} = pr.cumvar;  res.cv_cumvar{1, c} = pr.cv_cumvar;
    end
    fprintf('\n===== %s PC sweep: N = %d, k = %s =====\n', pop.name, pop.N, mat2str(k_list));

    t0 = tic;
    nI = nK * nP;
    data = pop.data;  labels = pop.labels;
    seeds = zeros(nI, 1);
    for ii = 1:nI
        [~, p] = ind2sub([nK, nP], ii);
        seeds(ii) = cfg.base_seed + model_id * 1e7 + p;      % independent of k
    end
    outs = cell(nI, 1);
    fprintf('[%s] %d k values x %d pairs = %d calls\n', pop.name, nK, nP, nI);
    q = progress_queue(sprintf('%s PC sweep', pop.name), nI);
    parfor ii = 1:nI
        [ik, p] = ind2sub([nK, nP], ii);
        o = struct('k', k_list(ik), 'n_perm', M, 'seed', seeds(ii), 'same_neurons', true);
        outs{ii} = pt_transfer_core(data{src(p)}, data{tgt(p)}, labels, o);
        send(q, 1);
    end
    res = unpack(res, outs, [nK, nP], F, M);
    res.done = true;  res.runtime_s = toc(t0);  res.updated = char(datetime('now'));
    groups = [{'pooled'}; unique(res.pair_group, 'stable')];
    sets = [{1:nP}; cellfun(@(g) find(strcmp(res.pair_group, g))', groups(2:end), 'UniformOutput', false)];
    res.summary = sweep_summary(res, k_list, groups, sets);
    save_atomic(file, res);
    print_k5(res, pop.name);
end
fprintf('\nDone in %.1f min. Results in %s\n', toc(t_all) / 60, out_dir);

%% ===========================================================================
function id = gabor_model_id(model)
% fixed index of a model, so seeds do not depend on which models a run includes
id = find(strcmp({'simple', 'complex'}, model));
assert(~isempty(id), 'Unknown Gabor model %s', model);
end

function k_list = k_grid(k_base, k_cap)
% values of k_base below the maximum rank, plus the maximum itself
k_list = unique([k_base(k_base < k_cap), k_cap]);
end

function g = pair_group_names(pairs)
names = containers.Map({'ac|ec', 'ec|ex', 'ac|ex'}, {'EC-AC', 'EC-EX', 'AC-EX'});
g = cell(size(pairs, 1), 1);
for i = 1:size(pairs, 1)
    g{i} = names(strjoin(sort(pairs(i, :)), '|'));
end
end

function res = unpack(res, outs, dims, F, M)
% outs{ii} -> res.<field>(ik, p, :) with [ik, p] = ind2sub(dims, ii)
for fn = {'self', 'nopt', 'pt_rot', 'pt_full', 'd_rot', 'd_full'}
    res.(fn{1}) = nan(dims(1), dims(2), F);
end
res.ctrl    = nan(dims(1), dims(2), F, M);
res.ctrl_id = nan(dims(1), dims(2), F, M);
for ii = 1:numel(outs)
    [ik, p] = ind2sub(dims, ii);
    o = outs{ii};
    for fn = {'self', 'nopt', 'pt_rot', 'pt_full', 'd_rot', 'd_full'}
        res.(fn{1})(ik, p, :) = o.(fn{1});
    end
    res.ctrl(ik, p, :, :)    = reshape(o.ctrl,    [1 1 F M]);
    res.ctrl_id(ik, p, :, :) = reshape(o.ctrl_id, [1 1 F M]);
end
end

function T = sweep_summary(res, k_list, groups, sets)
% k is the first dimension, the pair the second
rows = {};
for ik = 1:numel(k_list)
    for g = 1:numel(groups)
        S = sets{g};
        m = transfer_metrics(res.self(ik, S, :), res.pt_rot(ik, S, :), res.ctrl(ik, S, :, :), res.nopt(ik, S, :));
        rows(end + 1, :) = {k_list(ik), groups{g}, m.self, m.nopt, m.pt, mean(res.pt_full(ik, S, :), 'all', 'omitnan'), ...
            m.ctrl, mean(res.ctrl_id(ik, S, :, :), 'all', 'omitnan'), m.ratio, m.diff, m.headroom}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'k', 'set', 'self', 'nopt', 'pt_rot', 'pt_full', ...
    'ctrl', 'ctrl_id', 'pt_over_ctrl', 'pt_minus_ctrl', 'headroom'});
end

function print_k5(res, name)
ik = find(res.k_list == 5);
if isempty(ik), return; end
m = transfer_metrics(res.self(ik, :, :), res.pt_rot(ik, :, :), res.ctrl(ik, :, :, :), res.nopt(ik, :, :));
fprintf('[%s] k = 5: self %.3f  nonPT %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f\n', ...
    name, m.self, m.nopt, m.pt, m.ctrl, m.ratio, m.headroom);
end

function res = open_results(file, sig, overwrite)
res = [];
if exist(file, 'file') && ~overwrite
    L = load(file, 'res');
    res = L.res;
    if ~isequal(res.sig, sig)
        error(['%s was made with a different configuration. Set cfg.overwrite = true ' ...
               'to replace it, or set cfg.tag to write elsewhere.'], file);
    end
    fprintf('Resuming %s (%d of %d units done)\n', file, nnz(res.done), numel(res.done));
end
end

function save_atomic(file, res)
tmp = [file(1:end - 4) '_tmp.mat'];
save(tmp, 'res', '-v7');
movefile(tmp, file, 'f');
end
