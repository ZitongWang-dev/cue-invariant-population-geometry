%{
Filename: run_cue_transfer_incremental.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Cue-transfer decoding as a function of population size, in the full neuron
    space, for all four recorded populations (default). One run feeds:
      - Fig 3b: one point of each curve (V1 at its full size; V2 subsampled to
        the V1 count of the same monkey: FR 48, KO 109);
      - Fig 4c: self-decoding, non-PT, PT, corrected control and headroom vs N;
      - participation ratio vs N (raw and noise-corrected).

    N grid per population = 10:10:N_total, N_total, and the anchors of that
    monkey (FR: 24 48, KO: 54 109; 24/54 = the within-area halves size, so cue
    transfer has a size-matched reference for the halves analysis; 48/109 = the
    V1 size). Per N, cfg.n_draws random neuron subsets (1 at full size, where
    every subset is the whole population). Within a draw, the same neurons are
    used for all six ordered cue pairs, each decoded by pt_transfer_core.m
    (10-fold trial partition, rotation-only and full PT, corrected control with
    cfg.n_perm permutations per fold, non-PT on the same neurons).

Run:
    From MATLAB, run this script (any current folder). By default it runs all
    four populations. To split the work across machines, set an override in the
    workspace first, e.g.
        cfg_override = struct('populations', {{'KO','V2'; 'FR','V1'}});
        run_cue_transfer_incremental
    Any cfg field below can be overridden this way (unknown fields are refused).
    Re-running the same configuration resumes: finished N are skipped. Clear
    cfg_override (clear cfg_override) before running another driver.

Output: results/decoding_outputs/cue_transfer_incremental/<monkey>_<area>.mat,
    one variable `res` per population. Accuracies are kept per fold, so any
    pooling can be redone later:
      res.N_list [nN x 1], res.n_draws [nN x 1], res.done [nN x 1]
      res.pair_labels {6 x 1}  'ac->ec' = source (transformed, tested) ac,
                               target (decoder trained) ec
      res.pair_group  {6 x 1}  unordered pair: 'EC-AC', 'EC-EX', 'AC-EX'
      res.neuron_idx  {nN x n_draws} neurons used in each draw
      [nN x 6 pairs x n_draws x 10 folds]:
        self      target self-decoding (fold model on the target's held-out trials)
        nopt      non-PT: source held-out trials into the target model, no transform
        pt_rot    PT, rotation only (the headline, as in the paper)
        pt_full   PT, full similarity transform b*X*T + c
        d_rot, d_full   in-sample Procrustes residuals of the true fit
      [nN x 6 x n_draws x 10 x n_perm]:
        ctrl      corrected control, scored against the target labels (the floor)
        ctrl_id   the same predictions scored against the permuted identities (~1/50)
      res.pr_raw, res.pr_cv [nN x 3 cues x n_draws]   participation ratio per cue
      res.runtime_s [nN x 1]
      res.summary  table, one row per N x {pooled, EC-AC, EC-EX, AC-EX}; ratios
                   are ratios of pooled means (see transfer_metrics.m)
    Unused draw slots (full-size N) stay NaN.
%}

%% Initialization
clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
addpath(script_dir);

%% Configuration
cfg = struct();
cfg.populations = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};   % default: all four
cfg.timewindow  = [330 630];        % ms of the stored trace (30-330 ms after stimulus onset)
cfg.N_step      = 10;               % grid step of the neuron sweep
cfg.anchors     = struct('FR', [24 48], 'KO', [54 109]);   % extra grid points per monkey
cfg.N_list      = [];               % [] = automatic grid; or an explicit list (capped at N_total)
cfg.n_draws     = 20;               % neuron subsets per N (1 at full size)
cfg.n_perm      = 10;               % control permutations per fold
cfg.pr_splits   = 20;               % half-splits of the noise-corrected PR
cfg.pairs       = {'ac','ec'; 'ec','ac'; 'ec','ex'; 'ex','ec'; 'ac','ex'; 'ex','ac'};   % {source, target}
cfg.base_seed   = 1000;
cfg.out_root    = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag         = 'cue_transfer_incremental';
cfg.overwrite   = false;            % false = resume an existing file with the same configuration

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

all_pops = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};   % fixed order: seeds do not depend on cfg.populations
nP = size(cfg.pairs, 1);
t_all = tic;

%% Run
for ip = 1:size(cfg.populations, 1)
    monkey = cfg.populations{ip, 1};
    area   = cfg.populations{ip, 2};
    pop_id = find(strcmp(all_pops(:, 1), monkey) & strcmp(all_pops(:, 2), area));
    assert(~isempty(pop_id), 'Unknown population %s %s', monkey, area);

    pop = load_population(monkey, area, struct('timewindow', cfg.timewindow));
    Ntot = pop.N;
    F = pop.n_rep;
    M = cfg.n_perm;
    src = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 1));
    tgt = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 2));

    % N grid and draws per N
    if isempty(cfg.N_list)
        anchors = cfg.anchors.(monkey);
        N_list = unique([cfg.N_step:cfg.N_step:Ntot, Ntot, anchors(anchors <= Ntot)]);
    else
        N_list = unique(min(cfg.N_list(:)', Ntot));
    end
    nN = numel(N_list);
    D_list = repmat(cfg.n_draws, nN, 1);
    D_list(N_list == Ntot) = 1;
    D = max(D_list);

    % results file: resume it or start a new one
    file = fullfile(out_dir, sprintf('%s_%s.mat', monkey, area));
    sig  = struct('timewindow', cfg.timewindow, 'N_list', N_list, 'n_draws', cfg.n_draws, ...
        'n_perm', M, 'pr_splits', cfg.pr_splits, 'pairs', {cfg.pairs}, 'base_seed', cfg.base_seed);
    res = open_results(file, sig, cfg.overwrite);
    if isempty(res)
        res = struct();
        res.info        = 'Cue-transfer decoding vs N, full space. Field definitions: header of run_cue_transfer_incremental.m.';
        res.population  = pop.name;  res.monkey = monkey;  res.area = area;  res.N_total = Ntot;
        res.sig         = sig;
        res.cfg         = cfg;
        res.created     = char(datetime('now'));
        res.N_list      = N_list(:);
        res.n_draws     = D_list;
        res.pair_labels = strcat(cfg.pairs(:, 1), '->', cfg.pairs(:, 2));
        res.pair_group  = pair_group_names(cfg.pairs);
        res.cue_names   = pop.cue_names;
        res.done        = false(nN, 1);
        res.neuron_idx  = cell(nN, D);
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
    fprintf('\n===== %s: N = %s | %d draws, %d folds, %d permutations =====\n', ...
        pop.name, mat2str(N_list), cfg.n_draws, F, M);
    calls_left = sum(D_list(~res.done)) * nP;    % for the time-left estimate of this population
    calls_done = 0;
    t_pop = tic;

    for iN = 1:nN
        if res.done(iN), continue; end
        N  = N_list(iN);
        Dn = D_list(iN);
        t0 = tic;
        fprintf('[%s] point %d/%d: N = %d, %d draws x %d pairs = %d calls\n', pop.name, iN, nN, N, Dn, nP, Dn * nP);

        % neuron subsets (shared by the six pairs of a draw) and their PR
        idx = cell(1, Dn);
        for d = 1:Dn
            if N == Ntot
                idx{d} = 1:Ntot;
            else
                rsn = RandStream('twister', 'Seed', seed_of(cfg.base_seed, pop_id, N, d, 0));
                idx{d} = sort(randperm(rsn, Ntot, N));
            end
            for c = 1:3
                pr = population_pr(zscore(pop.data{c}(:, idx{d})), pop.labels, ...
                    struct('n_split', cfg.pr_splits, 'seed', seed_of(cfg.base_seed, pop_id, N, d, 10 + c)));
                res.pr_raw(iN, c, d) = pr.raw;
                res.pr_cv(iN, c, d)  = pr.cv;
            end
        end

        % decode every (draw, pair) in parallel
        nI = Dn * nP;
        seeds = zeros(nI, 1);
        for ii = 1:nI
            [d, p] = ind2sub([Dn, nP], ii);
            seeds(ii) = seed_of(cfg.base_seed, pop_id, N, d, p);
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
        res.neuron_idx(iN, 1:Dn) = idx;
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
function s = seed_of(base, pop_id, N, d, item)
% unique seed per (population, N, draw, item); item < 100, d < 100, N < 1000
s = base + pop_id * 1e7 + N * 1e4 + d * 100 + item;
end

function g = pair_group_names(pairs)
% unordered rendering pair of each ordered pair, named as in the paper
names = containers.Map({'ac|ec', 'ec|ex', 'ac|ex'}, {'EC-AC', 'EC-EX', 'AC-EX'});
g = cell(size(pairs, 1), 1);
for i = 1:size(pairs, 1)
    key = strjoin(sort(pairs(i, :)), '|');
    g{i} = names(key);
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
