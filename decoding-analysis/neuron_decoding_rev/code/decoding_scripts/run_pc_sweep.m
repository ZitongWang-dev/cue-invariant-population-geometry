%{
Filename: run_pc_sweep.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Transfer decoding as a function of the number of retained principal
    components k, on whole populations, for the three comparisons of the main
    figures (default: all three):
      'cue'     cue transfer, same neurons, six ordered rendering pairs per area
                (Fig 3b at k = 5, k sweep in Extended Data);
      'halves'  two disjoint halves of one area (floor(N/2) neurons each),
                same rendering, cfg.n_draws_halves random splits (Fig 4a at k = 5);
      'across'  all V1 neurons vs all V2 neurons, both directions, 3 x 3
                rendering combinations (Fig 4b at k = 5).

    PC decoding is done as in neuronal_decoding/.../procrustes_decoding_basic_incremental_pc.m:
    each side is projected onto the top-k principal directions of ITS OWN
    centered class-mean matrix, and the Procrustes rotation is a k x k map with
    k(k-1)/2 free parameters. Changes from that script, all inside
    pt_transfer_core.m: a 10-fold trial partition; the target basis, means and
    SVM come from training trials only, so self-decoding no longer leaks the
    test trials into the basis; cfg.n_perm control permutations per fold;
    non-PT (cue transfer only) projects the source trials onto the TARGET basis.

    k grid: cfg.k_base (default 1:10, 15:5:45, 49), with values above the
    comparison's maximum rank removed and that maximum added, e.g. FR V1 cue
    transfer ends at 48 and FR V1 halves (24 neurons) end at 24. The fold
    partition and the control permutations depend on the item and draw but not
    on k, so the k curves share their random numbers.

Run:
    By default it runs all three comparisons on all populations / both monkeys.
    To split the work across machines, e.g.
        cfg_override = struct('comparisons', {{'halves'}}, 'populations', {{'KO','V2'; 'FR','V1'}});
        run_pc_sweep
    'populations' applies to 'cue' and 'halves', 'monkeys' to 'across'.
    Re-running the same configuration resumes. Clear cfg_override before
    running another driver.

Output (variable `res` in each file):
    results/decoding_outputs/pc_sweep_cue/<monkey>_<area>.mat
        res.k_list [nK x 1]; res.pair_labels, res.pair_group {6 x 1}
        [nK x 6 pairs x 10 folds]: self, nopt, pt_rot, pt_full, d_rot, d_full
        [nK x 6 x 10 x n_perm]: ctrl, ctrl_id
        res.pr_raw, res.pr_cv [1 x 3 cues]; res.cumvar, res.cv_cumvar {1 x 3}
            (variance captured by the top k PCs of cue c = res.cumvar{c}(k))
    results/decoding_outputs/pc_sweep_halves/<monkey>_<area>.mat
        res.half_size; res.half_A, res.half_B {1 x n_draws}; res.done [n_draws x 1]
        [nK x 3 renderings x n_draws x 10]: self, pt_rot, pt_full, d_rot, d_full (nopt = NaN)
        [nK x 3 x n_draws x 10 x n_perm]: ctrl, ctrl_id
        res.pr_raw, res.pr_cv [3 renderings x n_draws x 2 halves]
    results/decoding_outputs/pc_sweep_across/<monkey>.mat
        res.dir_labels; res.combo_source_cue, res.combo_target_cue {9 x 1}; res.done [2 x 1]
        [2 directions x 9 combos x nK x 10]: self, pt_rot, pt_full, d_rot, d_full (nopt = NaN)
        [2 x 9 x nK x 10 x n_perm]: ctrl, ctrl_id
        res.pr_raw, res.pr_cv [2 areas x 3 cues]
    Each file also holds res.summary: one row per k (pooled, and per pair /
    rendering / direction) with ratios of pooled means (transfer_metrics.m).
%}

%% Initialization
clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
addpath(script_dir);

%% Configuration
cfg = struct();
cfg.comparisons    = {'cue', 'halves', 'across'};                    % default: all three
cfg.populations    = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};   % for 'cue' and 'halves'
cfg.monkeys        = {'FR', 'KO'};                                   % for 'across'
cfg.timewindow     = [330 630];
cfg.k_base         = [1:10, 15:5:45, 49];
cfg.n_draws_halves = 20;            % random half-splits ('halves' only)
cfg.n_perm         = 10;
cfg.pr_splits      = 20;
cfg.pairs          = {'ac','ec'; 'ec','ac'; 'ec','ex'; 'ex','ec'; 'ac','ex'; 'ex','ac'};   % {source, target}
cfg.base_seed      = 5000;
cfg.out_root       = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag_prefix     = 'pc_sweep';    % folders <tag_prefix>_cue / _halves / _across
cfg.overwrite      = false;

if exist('cfg_override', 'var') && isstruct(cfg_override)
    fo  = fieldnames(cfg_override);
    bad = setdiff(fo, fieldnames(cfg));
    assert(isempty(bad), 'cfg_override has unknown field(s): %s', strjoin(bad', ', '));
    for i = 1:numel(fo), cfg.(fo{i}) = cfg_override.(fo{i}); end
    fprintf('cfg_override applied: %s\n', strjoin(fo', ', '));
end
bad = setdiff(cfg.comparisons, {'cue', 'halves', 'across'});
assert(isempty(bad), 'Unknown comparison(s): %s', strjoin(bad, ', '));

if isempty(gcp('nocreate')), parpool('Processes'); end
t_all = tic;

%% Run
for ic = 1:numel(cfg.comparisons)
    out_dir = fullfile(cfg.out_root, sprintf('%s_%s', cfg.tag_prefix, cfg.comparisons{ic}));
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    switch cfg.comparisons{ic}
        case 'cue'
            for ip = 1:size(cfg.populations, 1)
                run_cue(cfg, cfg.populations{ip, 1}, cfg.populations{ip, 2}, out_dir);
            end
        case 'halves'
            for ip = 1:size(cfg.populations, 1)
                run_halves(cfg, cfg.populations{ip, 1}, cfg.populations{ip, 2}, out_dir);
            end
        case 'across'
            for im = 1:numel(cfg.monkeys)
                run_across(cfg, cfg.monkeys{im}, out_dir);
            end
    end
end
fprintf('\nDone in %.1f min.\n', toc(t_all) / 60);

%% ===========================================================================
function run_cue(cfg, monkey, area, out_dir)
% Cue transfer on the whole population of one area, every k, six ordered pairs.
pop_id = population_id(monkey, area);
pop = load_population(monkey, area, struct('timewindow', cfg.timewindow));
F = pop.n_rep;  M = cfg.n_perm;  nP = size(cfg.pairs, 1);
k_list = k_grid(cfg.k_base, min([pop.n_stim - 1, pop.N]));
nK = numel(k_list);
src = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 1));
tgt = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 2));

file = fullfile(out_dir, sprintf('%s_%s.mat', monkey, area));
sig  = struct('timewindow', cfg.timewindow, 'k_list', k_list, 'n_perm', M, ...
    'pr_splits', cfg.pr_splits, 'pairs', {cfg.pairs}, 'base_seed', cfg.base_seed);
res = open_results(file, sig, cfg.overwrite);
if ~isempty(res) && all(res.done)
    fprintf('%s cue PC sweep already done, skipped.\n', pop.name);
    return;
end
res = struct();
res.info        = 'PC sweep, cue transfer, whole population. Field definitions: header of run_pc_sweep.m.';
res.population  = pop.name;  res.monkey = monkey;  res.area = area;  res.N = pop.N;
res.sig         = sig;  res.cfg = cfg;  res.created = char(datetime('now'));
res.k_list      = k_list(:);
res.pair_labels = strcat(cfg.pairs(:, 1), '->', cfg.pairs(:, 2));
res.pair_group  = pair_group_names(cfg.pairs);
res.cue_names   = pop.cue_names;
res.done        = false;
for c = 1:3
    pr = population_pr(zscore(pop.data{c}), pop.labels, struct('n_split', cfg.pr_splits, 'seed', cfg.base_seed + c));
    res.pr_raw(1, c) = pr.raw;  res.pr_cv(1, c) = pr.cv;
    res.cumvar{1, c} = pr.cumvar;  res.cv_cumvar{1, c} = pr.cv_cumvar;
end
fprintf('\n===== %s cue PC sweep: N = %d, k = %s =====\n', pop.name, pop.N, mat2str(k_list));

t0 = tic;
nI = nK * nP;
data = pop.data;  labels = pop.labels;
seeds = zeros(nI, 1);
for ii = 1:nI
    [~, p] = ind2sub([nK, nP], ii);
    seeds(ii) = cfg.base_seed + pop_id * 1e7 + p;        % independent of k
end
outs = cell(nI, 1);
fprintf('[%s cue] %d k values x %d pairs = %d calls\n', pop.name, nK, nP, nI);
q = progress_queue(sprintf('%s cue PC sweep', pop.name), nI);
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
print_k5(res, pop.name, 'cue');
end

%% ===========================================================================
function run_halves(cfg, monkey, area, out_dir)
% Halves of one area (floor(N/2) each), every k, three renderings, n_draws splits.
pop_id = population_id(monkey, area);
pop = load_population(monkey, area, struct('timewindow', cfg.timewindow));
F = pop.n_rep;  M = cfg.n_perm;  D = cfg.n_draws_halves;  nR = numel(pop.cue_names);
h = floor(pop.N / 2);
k_list = k_grid(cfg.k_base, min([pop.n_stim - 1, h]));
nK = numel(k_list);

file = fullfile(out_dir, sprintf('%s_%s.mat', monkey, area));
sig  = struct('timewindow', cfg.timewindow, 'k_list', k_list, 'n_draws', D, 'n_perm', M, ...
    'pr_splits', cfg.pr_splits, 'base_seed', cfg.base_seed);
res = open_results(file, sig, cfg.overwrite);
if isempty(res)
    res = struct();
    res.info       = 'PC sweep, within-area halves. Field definitions: header of run_pc_sweep.m.';
    res.population = pop.name;  res.monkey = monkey;  res.area = area;  res.N = pop.N;
    res.half_size  = h;
    res.sig        = sig;  res.cfg = cfg;  res.created = char(datetime('now'));
    res.k_list     = k_list(:);
    res.renderings = pop.cue_names;
    res.done       = false(D, 1);
    res.half_A     = cell(1, D);  res.half_B = cell(1, D);
    for fn = {'self', 'nopt', 'pt_rot', 'pt_full', 'd_rot', 'd_full'}
        res.(fn{1}) = nan(nK, nR, D, F);
    end
    res.ctrl    = nan(nK, nR, D, F, M);
    res.ctrl_id = nan(nK, nR, D, F, M);
    res.pr_raw  = nan(nR, D, 2);  res.pr_cv = nan(nR, D, 2);
    res.runtime_s = nan(D, 1);
    res.summary = table();
end
fprintf('\n===== %s halves PC sweep: %d + %d neurons, k = %s, %d draws =====\n', pop.name, h, h, mat2str(k_list), D);

data = pop.data;  labels = pop.labels;
draws_left = nnz(~res.done);  draws_done = 0;  t_pop = tic;   % for the time-left estimate
for d = 1:D
    if res.done(d), continue; end
    t0 = tic;
    fprintf('[%s halves] draw %d/%d: %d k values x %d renderings = %d calls\n', pop.name, d, D, nK, nR, nK * nR);
    rsn  = RandStream('twister', 'Seed', cfg.base_seed + pop_id * 1e7 + d * 100);
    perm = randperm(rsn, pop.N);
    A = sort(perm(1:h));  B = sort(perm(h + 1:2 * h));
    for r = 1:nR
        for side = 1:2
            if side == 1, cols = A; else, cols = B; end
            pr = population_pr(zscore(data{r}(:, cols)), labels, ...
                struct('n_split', cfg.pr_splits, 'seed', cfg.base_seed + pop_id * 1e7 + d * 100 + 10 + 2 * r + side));
            res.pr_raw(r, d, side) = pr.raw;  res.pr_cv(r, d, side) = pr.cv;
        end
    end
    nI = nK * nR;
    seeds = zeros(nI, 1);
    for ii = 1:nI
        [~, r] = ind2sub([nK, nR], ii);
        seeds(ii) = cfg.base_seed + pop_id * 1e7 + d * 100 + r;   % independent of k
    end
    outs = cell(nI, 1);
    q = progress_queue(sprintf('%s halves draw %d', pop.name, d), nI);
    parfor ii = 1:nI
        [ik, r] = ind2sub([nK, nR], ii);
        o = struct('k', k_list(ik), 'n_perm', M, 'seed', seeds(ii), 'same_neurons', false);
        outs{ii} = pt_transfer_core(data{r}(:, A), data{r}(:, B), labels, o);
        send(q, 1);
    end
    for ii = 1:nI
        [ik, r] = ind2sub([nK, nR], ii);
        o = outs{ii};
        res.self(ik, r, d, :)    = o.self;
        res.nopt(ik, r, d, :)    = o.nopt;
        res.pt_rot(ik, r, d, :)  = o.pt_rot;
        res.pt_full(ik, r, d, :) = o.pt_full;
        res.d_rot(ik, r, d, :)   = o.d_rot;
        res.d_full(ik, r, d, :)  = o.d_full;
        res.ctrl(ik, r, d, :, :)    = reshape(o.ctrl,    [1 1 1 F M]);
        res.ctrl_id(ik, r, d, :, :) = reshape(o.ctrl_id, [1 1 1 F M]);
    end
    res.half_A{d} = A;  res.half_B{d} = B;
    res.done(d) = true;  res.runtime_s(d) = toc(t0);  res.updated = char(datetime('now'));
    groups = [{'pooled'}; upper(res.renderings(:))];
    sets = [{1:nR}; num2cell((1:nR)')];
    res.summary = sweep_summary(res, k_list, groups, sets);
    save_atomic(file, res);
    draws_done = draws_done + 1;
    fprintf('[%s halves] draw %2d/%d done | %.0f s | %s halves: ~%.1f h left\n', pop.name, d, D, res.runtime_s(d), ...
        pop.name, toc(t_pop) / draws_done * (draws_left - draws_done) / 3600);
end
if all(res.done), print_k5(res, pop.name, 'halves'); end
end

%% ===========================================================================
function run_across(cfg, monkey, out_dir)
% All V1 neurons vs all V2 neurons, every k, both directions, 9 rendering combinations.
mk_id = find(strcmp({'FR', 'KO'}, monkey));
assert(~isempty(mk_id), 'Unknown monkey %s', monkey);
pops = {load_population(monkey, 'V1', struct('timewindow', cfg.timewindow)), ...
        load_population(monkey, 'V2', struct('timewindow', cfg.timewindow))};
assert(isequal(pops{1}.labels, pops{2}.labels), 'V1 and V2 trial layouts differ.');
F = pops{1}.n_rep;  M = cfg.n_perm;
cues = pops{1}.cue_names;
[cs, ct] = ndgrid(1:3, 1:3);  cs = cs(:);  ct = ct(:);  nC = numel(cs);
k_list = k_grid(cfg.k_base, min([pops{1}.n_stim - 1, pops{1}.N, pops{2}.N]));
nK = numel(k_list);
dir_labels = {'V1-to-V2 (train V1, test V2)'; 'V2-to-V1 (train V2, test V1)'};
dir_target = [1; 2];  dir_source = [2; 1];

file = fullfile(out_dir, sprintf('%s.mat', monkey));
sig  = struct('timewindow', cfg.timewindow, 'k_list', k_list, 'n_perm', M, ...
    'pr_splits', cfg.pr_splits, 'base_seed', cfg.base_seed);
res = open_results(file, sig, cfg.overwrite);
if isempty(res)
    res = struct();
    res.info       = 'PC sweep, across areas, whole populations. Field definitions: header of run_pc_sweep.m.';
    res.monkey     = monkey;  res.N_area = [pops{1}.N, pops{2}.N];
    res.sig        = sig;  res.cfg = cfg;  res.created = char(datetime('now'));
    res.k_list     = k_list(:);
    res.dir_labels = dir_labels;
    res.dir_source_area = {'V2'; 'V1'};  res.dir_target_area = {'V1'; 'V2'};
    res.combo_source_cue = cues(cs)';  res.combo_target_cue = cues(ct)';
    res.cue_names  = cues;
    res.done       = false(2, 1);
    for fn = {'self', 'nopt', 'pt_rot', 'pt_full', 'd_rot', 'd_full'}
        res.(fn{1}) = nan(2, nC, nK, F);
    end
    res.ctrl    = nan(2, nC, nK, F, M);
    res.ctrl_id = nan(2, nC, nK, F, M);
    res.runtime_s = nan(2, 1);
    res.summary = table();
    for a = 1:2
        for c = 1:3
            pr = population_pr(zscore(pops{a}.data{c}), pops{a}.labels, ...
                struct('n_split', cfg.pr_splits, 'seed', cfg.base_seed + mk_id * 1e7 + 10 + 3 * a + c));
            res.pr_raw(a, c) = pr.raw;  res.pr_cv(a, c) = pr.cv;
        end
    end
end
fprintf('\n===== %s across-area PC sweep: V1 %d vs V2 %d neurons, k = %s =====\n', ...
    monkey, pops{1}.N, pops{2}.N, mat2str(k_list));

labels = pops{1}.labels;
for di = 1:2
    if res.done(di), continue; end
    t0 = tic;
    Xsrc = pops{dir_source(di)}.data;  Xtgt = pops{dir_target(di)}.data;
    nI = nK * nC;
    seeds = zeros(nI, 1);
    for ii = 1:nI
        [~, c] = ind2sub([nK, nC], ii);
        seeds(ii) = cfg.base_seed + mk_id * 1e7 + di * 1e5 + c;   % independent of k
    end
    outs = cell(nI, 1);
    fprintf('[%s across] direction %d/2: %s, %d k values x %d rendering combinations = %d calls\n', ...
        monkey, di, dir_labels{di}, nK, nC, nI);
    q = progress_queue(sprintf('%s across %s', monkey, dir_labels{di}), nI);
    parfor ii = 1:nI
        [ik, c] = ind2sub([nK, nC], ii);
        o = struct('k', k_list(ik), 'n_perm', M, 'seed', seeds(ii), 'same_neurons', false);
        outs{ii} = pt_transfer_core(Xsrc{cs(c)}, Xtgt{ct(c)}, labels, o);
        send(q, 1);
    end
    for ii = 1:nI
        [ik, c] = ind2sub([nK, nC], ii);
        o = outs{ii};
        res.self(di, c, ik, :)    = o.self;
        res.nopt(di, c, ik, :)    = o.nopt;
        res.pt_rot(di, c, ik, :)  = o.pt_rot;
        res.pt_full(di, c, ik, :) = o.pt_full;
        res.d_rot(di, c, ik, :)   = o.d_rot;
        res.d_full(di, c, ik, :)  = o.d_full;
        res.ctrl(di, c, ik, :, :)    = reshape(o.ctrl,    [1 1 1 F M]);
        res.ctrl_id(di, c, ik, :, :) = reshape(o.ctrl_id, [1 1 1 F M]);
    end
    res.done(di) = true;  res.runtime_s(di) = toc(t0);  res.updated = char(datetime('now'));
    res.summary = across_summary(res);
    save_atomic(file, res);
    m = transfer_metrics(res.self(di, :, k_list == min(5, max(k_list)), :), res.pt_rot(di, :, k_list == min(5, max(k_list)), :), ...
        res.ctrl(di, :, k_list == min(5, max(k_list)), :, :));
    fprintf('[%s %s] k = 5: self %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f | %.0f s\n', ...
        monkey, dir_labels{di}, m.self, m.pt, m.ctrl, m.ratio, m.headroom, res.runtime_s(di));
end
end

%% ===========================================================================
function k_list = k_grid(k_base, k_cap)
% values of k_base below the maximum rank, plus the maximum itself
k_list = unique([k_base(k_base < k_cap), k_cap]);
end

function id = population_id(monkey, area)
all_pops = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};
id = find(strcmp(all_pops(:, 1), monkey) & strcmp(all_pops(:, 2), area));
assert(~isempty(id), 'Unknown population %s %s', monkey, area);
end

function g = pair_group_names(pairs)
names = containers.Map({'ac|ec', 'ec|ex', 'ac|ex'}, {'EC-AC', 'EC-EX', 'AC-EX'});
g = cell(size(pairs, 1), 1);
for i = 1:size(pairs, 1)
    g{i} = names(strjoin(sort(pairs(i, :)), '|'));
end
end

function res = unpack(res, outs, dims, F, M)
% cue sweep: outs{ii} -> res.<field>(ik, p, :) with [ik, p] = ind2sub(dims, ii)
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
% cue / halves sweeps: k is the first dimension, the pair or rendering the second
rows = {};
for ik = 1:numel(k_list)
    for g = 1:numel(groups)
        S = sets{g};
        m = transfer_metrics(slice2(res.self, ik, S), slice2(res.pt_rot, ik, S), slice2(res.ctrl, ik, S), slice2(res.nopt, ik, S));
        rows(end + 1, :) = {k_list(ik), groups{g}, m.self, m.nopt, m.pt, mean(slice2(res.pt_full, ik, S), 'all', 'omitnan'), ...
            m.ctrl, mean(slice2(res.ctrl_id, ik, S), 'all', 'omitnan'), m.ratio, m.diff, m.headroom}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'k', 'set', 'self', 'nopt', 'pt_rot', 'pt_full', ...
    'ctrl', 'ctrl_id', 'pt_over_ctrl', 'pt_minus_ctrl', 'headroom'});
end

function x = slice2(A, ik, S)
% A(ik, S, :, ...) for an array of any number of trailing dimensions
idx = repmat({':'}, 1, ndims(A));
idx{1} = ik;  idx{2} = S;
x = A(idx{:});
end

function T = across_summary(res)
rows = {};
nC = numel(res.combo_source_cue);
same = find(strcmp(res.combo_source_cue, res.combo_target_cue))';
for di = find(res.done)'
    for ik = 1:numel(res.k_list)
        sets = {1:nC, same};  names = {'all 9 combos', 'same rendering'};
        for g = 1:2
            C = sets{g};
            m = transfer_metrics(res.self(di, C, ik, :), res.pt_rot(di, C, ik, :), res.ctrl(di, C, ik, :, :));
            rows(end + 1, :) = {res.dir_labels{di}, res.k_list(ik), names{g}, m.self, m.pt, ...
                mean(res.pt_full(di, C, ik, :), 'all', 'omitnan'), m.ctrl, ...
                mean(res.ctrl_id(di, C, ik, :, :), 'all', 'omitnan'), m.ratio, m.diff, m.headroom}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'direction', 'k', 'set', 'self', 'pt_rot', 'pt_full', ...
    'ctrl', 'ctrl_id', 'pt_over_ctrl', 'pt_minus_ctrl', 'headroom'});
end

function print_k5(res, name, kind)
ik = find(res.k_list == 5);
if isempty(ik), return; end
m = transfer_metrics(slice2(res.self, ik, ':'), slice2(res.pt_rot, ik, ':'), slice2(res.ctrl, ik, ':'), slice2(res.nopt, ik, ':'));
fprintf('[%s %s] k = 5: self %.3f  nonPT %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f\n', ...
    name, kind, m.self, m.nopt, m.pt, m.ctrl, m.ratio, m.headroom);
end

function res = open_results(file, sig, overwrite)
res = [];
if exist(file, 'file') && ~overwrite
    L = load(file, 'res');
    res = L.res;
    if ~isequal(res.sig, sig)
        error(['%s was made with a different configuration. Set cfg.overwrite = true ' ...
               'to replace it, or set cfg.tag_prefix to write elsewhere.'], file);
    end
    fprintf('Resuming %s (%d of %d units done)\n', file, nnz(res.done), numel(res.done));
end
end

function save_atomic(file, res)
tmp = [file(1:end - 4) '_tmp.mat'];
save(tmp, 'res', '-v7');
movefile(tmp, file, 'f');
end
