%{
Filename: run_gabor_cue_transfer_incremental.m
Author:   Zitong Wang
Date:     2026-10-15

Description:
    Cue-transfer decoding of the Gabor-filter cell models as a function of the
    number of model cells, in the full cell space, for both models (default):
    simple (192 cells) and complex (48 cells). One run feeds:
      - Fig 15a: the 48-cell point of both models (simple down-sampled to the
        complex-cell count, as published);
      - Fig 15b: the whole models (simple 192, complex 48);
      - rebuttal point 9: self-decoding, PT, corrected control and headroom vs
        the number of cells, on the 10-cell grid of the neural sweep, with the
        participation ratio at every size (raw and noise-corrected).

    Grid per model = 10:10:N_total, N_total, and 48 (cfg.anchors): simple
    10, 20, ..., 190, 192 and 48; complex 10, 20, 30, 40, 48. Per size,
    cfg.n_draws random cell subsets (1 for the whole model). Within a draw the
    same cells are used for all six ordered rendering pairs, each decoded by
    neuron_decoding_rev/pt_transfer_core.m (10-fold trial partition,
    rotation-only and full PT, corrected control with cfg.n_perm permutations
    per fold, non-PT on the same cells), so model cells and neurons go through
    exactly the same decoder.

    Trials: the saved Poisson draws of the published analyses (load_gabor_model.m).

Run:
    From MATLAB, run this script (any current folder). By default it runs both
    models. To run one, set an override first, e.g.
        cfg_override = struct('models', {{'complex'}});
        run_gabor_cue_transfer_incremental
    Any cfg field below can be overridden this way (unknown fields are refused).
    Re-running the same configuration resumes: finished sizes are skipped. Clear
    cfg_override (clear cfg_override) before running another driver.

Output: results/decoding_outputs/gabor_cue_transfer_incremental/<model>.mat,
    one variable `res` per model, with the same fields as the neuron driver
    (neuron_decoding_rev/run_cue_transfer_incremental.m):
      res.model, res.population ('Gabor simple'), res.N_total, res.n_silent
      res.fig15a_N             the Fig 15a size (48)
      res.N_list [nN x 1], res.n_draws [nN x 1], res.done [nN x 1]
      res.pair_labels {6 x 1}  'ac->ec' = source (transformed, tested) ac,
                               target (decoder trained) ec
      res.pair_group  {6 x 1}  unordered pair: 'EC-AC', 'EC-EX', 'AC-EX'
      res.unit_idx    {nN x n_draws} cells used in each draw
      [nN x 6 pairs x n_draws x 10 folds]: self, nopt, pt_rot, pt_full, d_rot, d_full
      [nN x 6 x n_draws x 10 x n_perm]:    ctrl, ctrl_id
      res.pr_raw, res.pr_cv [nN x 3 cues x n_draws]
      res.runtime_s [nN x 1]
      res.summary  table, one row per size x {pooled, EC-AC, EC-EX, AC-EX}
    Unused draw slots (whole model) stay NaN.
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
cfg.N_step    = 10;                 % grid step (as the neural sweep and the rebuttal's Gabor sweep)
cfg.anchors   = 48;                 % extra grid points: the Fig 15a size (complex-cell count)
cfg.N_list    = [];                 % [] = automatic grid; or an explicit list (capped at the model size)
cfg.n_draws   = 20;                 % cell subsets per size (1 for the whole model)
cfg.n_perm    = 10;                 % control permutations per fold
cfg.pr_splits = 20;                 % half-splits of the noise-corrected PR
cfg.pairs     = {'ac','ec'; 'ec','ac'; 'ec','ex'; 'ex','ec'; 'ac','ex'; 'ex','ac'};   % {source, target}
cfg.base_seed = 31000;
cfg.data_root = '';                 % '' = Gabor_decoding/.../sampled_trials (load_gabor_model default)
cfg.out_root  = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag       = 'gabor_cue_transfer_incremental';
cfg.overwrite = false;              % false = resume an existing file with the same configuration

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

nP = size(cfg.pairs, 1);
t_all = tic;

%% Run
for im = 1:numel(cfg.models)
    model    = lower(cfg.models{im});
    model_id = gabor_model_id(model);

    pop  = load_gabor_model(model, struct('data_root', cfg.data_root));
    Ntot = pop.N;
    F = pop.n_rep;
    M = cfg.n_perm;
    src = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 1));
    tgt = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 2));

    % cell grid and draws per size
    if isempty(cfg.N_list)
        anchors = cfg.anchors(cfg.anchors <= Ntot);
        N_list = unique([cfg.N_step:cfg.N_step:Ntot, Ntot, anchors]);
    else
        N_list = unique(min(cfg.N_list(:)', Ntot));
    end
    nN = numel(N_list);
    D_list = repmat(cfg.n_draws, nN, 1);
    D_list(N_list == Ntot) = 1;
    D = max(D_list);

    % results file: resume it or start a new one
    file = fullfile(out_dir, sprintf('%s.mat', model));
    [~, data_name] = fileparts(pop.file);    % file name only, so another machine can resume the file
    sig  = struct('N_list', N_list, 'n_draws', cfg.n_draws, 'n_perm', M, 'pr_splits', cfg.pr_splits, ...
        'pairs', {cfg.pairs}, 'base_seed', cfg.base_seed, 'data_file', data_name);
    res = open_results(file, sig, cfg.overwrite);
    if isempty(res)
        res = struct();
        res.info        = 'Gabor-model cue-transfer decoding vs number of cells, full space. Field definitions: header of run_gabor_cue_transfer_incremental.m.';
        res.population  = pop.name;  res.model = model;
        res.N_total     = Ntot;
        res.n_silent    = pop.n_silent;
        res.fig15a_N    = 48;
        res.sig         = sig;
        res.cfg         = cfg;
        res.created     = char(datetime('now'));
        res.N_list      = N_list(:);
        res.n_draws     = D_list;
        res.pair_labels = strcat(cfg.pairs(:, 1), '->', cfg.pairs(:, 2));
        res.pair_group  = pair_group_names(cfg.pairs);
        res.cue_names   = pop.cue_names;
        res.done        = false(nN, 1);
        res.unit_idx    = cell(nN, D);
        for fn = {'self', 'nopt', 'pt_rot', 'pt_full', 'd_rot', 'd_full'}
            res.(fn{1}) = nan(nN, nP, D, F);
        end
        res.ctrl        = nan(nN, nP, D, F, M);
        res.ctrl_id     = nan(nN, nP, D, F, M);
        res.pr_raw      = nan(nN, 3, D);
        res.pr_cv       = nan(nN, 3, D);
        res.runtime_s   = nan(nN, 1);
        res.summary     = table();
    end
    fprintf('\n===== %s (%d cells; constant within AC/EC/EX: %s): N = %s | %d draws, %d folds, %d permutations =====\n', ...
        pop.name, Ntot, mat2str(pop.n_silent), mat2str(N_list), cfg.n_draws, F, M);
    calls_left = sum(D_list(~res.done)) * nP;    % for the time-left estimate of this model
    calls_done = 0;
    t_pop = tic;

    for iN = 1:nN
        if res.done(iN), continue; end
        N  = N_list(iN);
        Dn = D_list(iN);
        t0 = tic;
        fprintf('[%s] point %d/%d: N = %d, %d draws x %d pairs = %d calls\n', pop.name, iN, nN, N, Dn, nP, Dn * nP);

        % cell subsets (shared by the six pairs of a draw) and their PR
        idx = cell(1, Dn);
        for d = 1:Dn
            if N == Ntot
                idx{d} = 1:Ntot;
            else
                rsn = RandStream('twister', 'Seed', seed_of(cfg.base_seed, model_id, N, d, 0));
                idx{d} = sort(randperm(rsn, Ntot, N));
            end
            for c = 1:3
                pr = population_pr(zscore(pop.data{c}(:, idx{d})), pop.labels, ...
                    struct('n_split', cfg.pr_splits, 'seed', seed_of(cfg.base_seed, model_id, N, d, 10 + c)));
                res.pr_raw(iN, c, d) = pr.raw;
                res.pr_cv(iN, c, d)  = pr.cv;
            end
        end

        % decode every (draw, pair) in parallel
        nI = Dn * nP;
        seeds = zeros(nI, 1);
        for ii = 1:nI
            [d, p] = ind2sub([Dn, nP], ii);
            seeds(ii) = seed_of(cfg.base_seed, model_id, N, d, p);
        end
        data = pop.data;  labels = pop.labels;
        outs = cell(nI, 1);
        q = progress_queue(sprintf('%s N = %d', pop.name, N), nI);
        parfor ii = 1:nI
            [d, p] = ind2sub([Dn, nP], ii);
            o = struct('n_perm', M, 'seed', seeds(ii), 'same_neurons', true);
            outs{ii} = pt_transfer_core(data{src(p)}(:, idx{d}), data{tgt(p)}(:, idx{d}), labels, o);
            send(q, 1);
        end

        for ii = 1:nI
            [d, p] = ind2sub([Dn, nP], ii);
            o = outs{ii};
            res.self(iN, p, d, :)    = o.self;
            res.nopt(iN, p, d, :)    = o.nopt;
            res.pt_rot(iN, p, d, :)  = o.pt_rot;
            res.pt_full(iN, p, d, :) = o.pt_full;
            res.d_rot(iN, p, d, :)   = o.d_rot;
            res.d_full(iN, p, d, :)  = o.d_full;
            res.ctrl(iN, p, d, :, :)    = reshape(o.ctrl,    [1 1 1 F M]);
            res.ctrl_id(iN, p, d, :, :) = reshape(o.ctrl_id, [1 1 1 F M]);
        end
        res.unit_idx(iN, 1:Dn) = idx;
        res.done(iN)      = true;
        res.runtime_s(iN) = toc(t0);
        res.updated       = char(datetime('now'));
        res.summary       = build_summary(res);
        save_atomic(file, res);

        m = transfer_metrics(res.self(iN, :, :, :), res.pt_rot(iN, :, :, :), res.ctrl(iN, :, :, :, :), res.nopt(iN, :, :, :));
        calls_done = calls_done + nI;
        fprintf('[%s] N = %3d (%2d draws): self %.3f  nonPT %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f  PR %.1f | %.0f s | %s: ~%.1f h left\n', ...
            pop.name, N, Dn, m.self, m.nopt, m.pt, m.ctrl, m.ratio, m.headroom, ...
            mean(res.pr_raw(iN, :, 1:Dn), 'all'), res.runtime_s(iN), ...
            pop.name, toc(t_pop) / calls_done * (calls_left - calls_done) / 3600);
    end
end
fprintf('\nDone in %.1f min. Results in %s\n', toc(t_all) / 60, out_dir);

%% ---------------------------------------------------------------------------
function id = gabor_model_id(model)
% fixed index of a model, so seeds do not depend on which models a run includes
id = find(strcmp({'simple', 'complex'}, model));
assert(~isempty(id), 'Unknown Gabor model %s', model);
end

function s = seed_of(base, model_id, N, d, item)
% unique seed per (model, size, draw, item); item < 100, d < 100, N < 1000
s = base + model_id * 1e7 + N * 1e4 + d * 100 + item;
end

function g = pair_group_names(pairs)
% unordered rendering pair of each ordered pair, named as in the paper
names = containers.Map({'ac|ec', 'ec|ex', 'ac|ex'}, {'EC-AC', 'EC-EX', 'AC-EX'});
g = cell(size(pairs, 1), 1);
for i = 1:size(pairs, 1)
    g{i} = names(strjoin(sort(pairs(i, :)), '|'));
end
end

function T = build_summary(res)
groups = [{'pooled'}; unique(res.pair_group, 'stable')];
rows = {};
for iN = find(res.done)'
    for g = 1:numel(groups)
        if g == 1, P = 1:numel(res.pair_group); else, P = find(strcmp(res.pair_group, groups{g}))'; end
        m = transfer_metrics(res.self(iN, P, :, :), res.pt_rot(iN, P, :, :), res.ctrl(iN, P, :, :, :), res.nopt(iN, P, :, :));
        rows(end + 1, :) = {res.N_list(iN), res.n_draws(iN), groups{g}, m.self, m.nopt, m.pt, ...
            mean(res.pt_full(iN, P, :, :), 'all', 'omitnan'), m.ctrl, ...
            mean(res.ctrl_id(iN, P, :, :, :), 'all', 'omitnan'), m.ratio, m.diff, m.headroom, ...
            mean(res.pr_raw(iN, :, :), 'all', 'omitnan'), mean(res.pr_cv(iN, :, :), 'all', 'omitnan')}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'N', 'n_draws', 'pair', 'self', 'nopt', 'pt_rot', 'pt_full', ...
    'ctrl', 'ctrl_id', 'pt_over_ctrl', 'pt_minus_ctrl', 'headroom', 'pr_raw', 'pr_cv'});
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
    fprintf('Resuming %s (%d of %d grid points done)\n', file, nnz(res.done), numel(res.done));
end
end

function save_atomic(file, res)
% write to a temporary file first, so an interrupted save never corrupts the results
tmp = [file(1:end - 4) '_tmp.mat'];
save(tmp, 'res', '-v7');
movefile(tmp, file, 'f');
end
