%{
Filename: run_halves_incremental.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Within-area ("half-and-half") transfer decoding as a function of half size,
    in the full neuron space, for all four recorded populations (default). Two
    disjoint random halves of one area's neurons respond to the SAME rendering;
    half A (source) is mapped onto half B (target) and decoded by a model
    trained on half B. One run feeds:
      - Fig 4a: one point of each curve (half size 24 for FR, 54 for KO, the
        size used in the paper for both V1 and V2);
      - Fig 13: self-decoding, PT, corrected control and headroom vs half size;
      - participation ratio of each half.

    Half-size grid per population = 5, 10:10:h_max, h_max = floor(N_total/2),
    and the monkey's anchor (FR 24, KO 54). Per size, cfg.n_draws random
    splits; within a draw the same split serves all three renderings.

    Both halves share their trial rows (same pseudo-population trials), so
    same-session neurons in different halves carry the same trial's noise.
    pt_transfer_core.m holds out the same rows on both sides, so the target half
    never trains on the trial the source half is tested on. The non-PT column
    here maps different neurons onto each other and is only a sanity check
    (it should sit at chance).

Run:
    By default it runs all four populations. To split the work across machines:
        cfg_override = struct('populations', {{'KO','V2'; 'FR','V1'}});
        run_halves_incremental
    Re-running the same configuration resumes (finished sizes are skipped).
    Clear cfg_override before running another driver.

Output: results/decoding_outputs/halves_incremental/<monkey>_<area>.mat, variable `res`:
      res.h_list [nH x 1], res.done [nH x 1], res.renderings {'ac','ec','ex'}
      res.half_A, res.half_B {nH x n_draws}   neuron indices of the two halves
      [nH x 3 renderings x n_draws x 10 folds]:
        self (target half), nopt (sanity, different neurons), pt_rot, pt_full, d_rot, d_full
      [nH x 3 x n_draws x 10 x n_perm]: ctrl, ctrl_id
      res.pr_raw, res.pr_cv [nH x 3 renderings x n_draws x 2 halves (A, B)]
      res.runtime_s [nH x 1]
      res.summary  table, one row per size x {pooled, AC, EC, EX}
%}

%% Initialization
clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
addpath(script_dir);

%% Configuration
cfg = struct();
cfg.populations = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};   % default: all four
cfg.timewindow  = [330 630];
cfg.h_step      = 10;               % grid step of the half size (5 is always included)
cfg.anchors     = struct('FR', 24, 'KO', 54);   % half size used in the paper (Fig 4a)
cfg.h_list      = [];               % [] = automatic grid; or an explicit list (capped at floor(N/2))
cfg.n_draws     = 20;               % random half-splits per size
cfg.n_perm      = 10;
cfg.pr_splits   = 20;
cfg.base_seed   = 3000;
cfg.out_root    = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag         = 'halves_incremental';
cfg.overwrite   = false;

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

all_pops = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};
t_all = tic;

%% Run
for ip = 1:size(cfg.populations, 1)
    monkey = cfg.populations{ip, 1};
    area   = cfg.populations{ip, 2};
    pop_id = find(strcmp(all_pops(:, 1), monkey) & strcmp(all_pops(:, 2), area));
    assert(~isempty(pop_id), 'Unknown population %s %s', monkey, area);

    pop  = load_population(monkey, area, struct('timewindow', cfg.timewindow));
    Ntot = pop.N;
    F = pop.n_rep;
    M = cfg.n_perm;
    nR = numel(pop.cue_names);
    h_max = floor(Ntot / 2);

    if isempty(cfg.h_list)
        anchor = cfg.anchors.(monkey);
        h_list = unique([5, cfg.h_step:cfg.h_step:h_max, h_max, anchor(anchor <= h_max)]);
    else
        h_list = unique(min(cfg.h_list(:)', h_max));
    end
    nH = numel(h_list);
    D  = cfg.n_draws;

    file = fullfile(out_dir, sprintf('%s_%s.mat', monkey, area));
    sig  = struct('timewindow', cfg.timewindow, 'h_list', h_list, 'n_draws', D, ...
        'n_perm', M, 'pr_splits', cfg.pr_splits, 'base_seed', cfg.base_seed);
    res = open_results(file, sig, cfg.overwrite);
    if isempty(res)
        res = struct();
        res.info       = 'Within-area halves decoding vs half size, full space. Field definitions: header of run_halves_incremental.m.';
        res.population = pop.name;  res.monkey = monkey;  res.area = area;  res.N_total = Ntot;
        res.sig        = sig;
        res.cfg        = cfg;
        res.created    = char(datetime('now'));
        res.h_list     = h_list(:);
        res.renderings = pop.cue_names;
        res.done       = false(nH, 1);
        res.half_A     = cell(nH, D);
        res.half_B     = cell(nH, D);
        for fn = {'self', 'nopt', 'pt_rot', 'pt_full', 'd_rot', 'd_full'}
            res.(fn{1}) = nan(nH, nR, D, F);
        end
        res.ctrl       = nan(nH, nR, D, F, M);
        res.ctrl_id    = nan(nH, nR, D, F, M);
        res.pr_raw     = nan(nH, nR, D, 2);
        res.pr_cv      = nan(nH, nR, D, 2);
        res.runtime_s  = nan(nH, 1);
        res.summary    = table();
    end
    fprintf('\n===== %s halves: h = %s | %d draws, %d folds, %d permutations =====\n', ...
        pop.name, mat2str(h_list), D, F, M);
    calls_left = nnz(~res.done) * D * nR;     % for the time-left estimate of this population
    calls_done = 0;
    t_pop = tic;

    for iH = 1:nH
        if res.done(iH), continue; end
        h  = h_list(iH);
        t0 = tic;
        fprintf('[%s] point %d/%d: half size %d, %d draws x %d renderings = %d calls\n', pop.name, iH, nH, h, D, nR, D * nR);

        % random half-splits (one per draw, shared by the three renderings) and their PR
        A = cell(1, D);  B = cell(1, D);
        for d = 1:D
            rsn  = RandStream('twister', 'Seed', seed_of(cfg.base_seed, pop_id, h, d, 0));
            perm = randperm(rsn, Ntot);
            A{d} = sort(perm(1:h));
            B{d} = sort(perm(h + 1:2 * h));
            for r = 1:nR
                for side = 1:2
                    if side == 1, cols = A{d}; else, cols = B{d}; end
                    pr = population_pr(zscore(pop.data{r}(:, cols)), pop.labels, ...
                        struct('n_split', cfg.pr_splits, 'seed', seed_of(cfg.base_seed, pop_id, h, d, 10 + 2 * r + side)));
                    res.pr_raw(iH, r, d, side) = pr.raw;
                    res.pr_cv(iH, r, d, side)  = pr.cv;
                end
            end
        end

        nI = D * nR;
        seeds = zeros(nI, 1);
        for ii = 1:nI
            [d, r] = ind2sub([D, nR], ii);
            seeds(ii) = seed_of(cfg.base_seed, pop_id, h, d, r);
        end
        data = pop.data;  labels = pop.labels;
        outs = cell(nI, 1);
        q = progress_queue(sprintf('%s half = %d', pop.name, h), nI);
        parfor ii = 1:nI
            [d, r] = ind2sub([D, nR], ii);
            o = struct('n_perm', M, 'seed', seeds(ii), 'same_neurons', false);
            outs{ii} = pt_transfer_core(data{r}(:, A{d}), data{r}(:, B{d}), labels, o);
            send(q, 1);
        end

        for ii = 1:nI
            [d, r] = ind2sub([D, nR], ii);
            o = outs{ii};
            res.self(iH, r, d, :)    = o.self;
            res.nopt(iH, r, d, :)    = o.nopt;
            res.pt_rot(iH, r, d, :)  = o.pt_rot;
            res.pt_full(iH, r, d, :) = o.pt_full;
            res.d_rot(iH, r, d, :)   = o.d_rot;
            res.d_full(iH, r, d, :)  = o.d_full;
            res.ctrl(iH, r, d, :, :)    = reshape(o.ctrl,    [1 1 1 F M]);
            res.ctrl_id(iH, r, d, :, :) = reshape(o.ctrl_id, [1 1 1 F M]);
        end
        res.half_A(iH, :) = A;
        res.half_B(iH, :) = B;
        res.done(iH)      = true;
        res.runtime_s(iH) = toc(t0);
        res.updated       = char(datetime('now'));
        res.summary       = build_summary(res);
        save_atomic(file, res);

        m = transfer_metrics(res.self(iH, :, :, :), res.pt_rot(iH, :, :, :), res.ctrl(iH, :, :, :, :), res.nopt(iH, :, :, :));
        calls_done = calls_done + nI;
        fprintf('[%s] half = %3d: self %.3f  nonPT(sanity) %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f  PR %.1f | %.0f s | %s: ~%.1f h left\n', ...
            pop.name, h, m.self, m.nopt, m.pt, m.ctrl, m.ratio, m.headroom, mean(res.pr_raw(iH, :, :, :), 'all'), res.runtime_s(iH), ...
            pop.name, toc(t_pop) / calls_done * (calls_left - calls_done) / 3600);
    end
end
fprintf('\nDone in %.1f min. Results in %s\n', toc(t_all) / 60, out_dir);

%% ---------------------------------------------------------------------------
function s = seed_of(base, pop_id, h, d, item)
s = base + pop_id * 1e7 + h * 1e4 + d * 100 + item;
end

function T = build_summary(res)
groups = [{'pooled'}; upper(res.renderings(:))];
rows = {};
for iH = find(res.done)'
    for g = 1:numel(groups)
        if g == 1, R = 1:numel(res.renderings); else, R = g - 1; end
        m = transfer_metrics(res.self(iH, R, :, :), res.pt_rot(iH, R, :, :), res.ctrl(iH, R, :, :, :), res.nopt(iH, R, :, :));
        rows(end + 1, :) = {res.h_list(iH), groups{g}, m.self, m.nopt, m.pt, ...
            mean(res.pt_full(iH, R, :, :), 'all', 'omitnan'), m.ctrl, ...
            mean(res.ctrl_id(iH, R, :, :, :), 'all', 'omitnan'), m.ratio, m.diff, m.headroom, ...
            mean(res.pr_raw(iH, R, :, :), 'all', 'omitnan'), mean(res.pr_cv(iH, R, :, :), 'all', 'omitnan')}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'half_size', 'rendering', 'self', 'nopt_sanity', 'pt_rot', 'pt_full', ...
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
tmp = [file(1:end - 4) '_tmp.mat'];
save(tmp, 'res', '-v7');
movefile(tmp, file, 'f');
end
