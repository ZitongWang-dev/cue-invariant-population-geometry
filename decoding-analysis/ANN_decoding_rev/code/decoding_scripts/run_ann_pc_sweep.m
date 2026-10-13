%{
Filename: run_ann_pc_sweep.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    CNN transfer decoding as a function of the number of retained principal
    components k, on whole layers, for two comparisons (default: both):
      'cue'        cue transfer within one pooling layer, same units, six ordered
                   rendering pairs (all eight layers by default);
      'hierarchy'  one layer vs another layer of the same network, both
                   directions, 3 x 3 rendering combinations (default AlexNet
                   pool2 <-> pool5, the CNN analogue of across-area decoding).

    The PC projection is the one of the neuron PC sweep
    (neuron_decoding_rev/run_pc_sweep.m), done inside the shared core
    pt_transfer_core.m: each side is projected onto the top-k principal
    directions of ITS OWN centered class-mean matrix (training trials only), and
    the rotation is a k x k map with k(k-1)/2 free parameters whatever the layer
    size, so CNN layers and neurons can be compared at the same rotation freedom
    (k = 5 in the main comparison). Non-PT (cue transfer only) projects the
    source trials onto the TARGET basis.

    k grid: cfg.k_base (default 1:10, 15:5:45, 49), with values above the
    comparison's maximum rank removed and that maximum added (every layer has at
    least 64 units, so the maximum is 49 = n_stim - 1). The fold partition and
    the control permutations depend on the item but not on k, so the k curves
    share their random numbers.

Run:
    By default it runs both comparisons. To split the work across machines, e.g.
        cfg_override = struct('comparisons', {{'cue'}}, 'layers', {{'vgg',4; 'vgg',5}});
        run_ann_pc_sweep
    'layers' applies to 'cue', 'layer_pairs' to 'hierarchy'. Re-running the same
    configuration resumes. Clear cfg_override before running another driver.

Output (variable `res` in each file):
    results/decoding_outputs/ann_pc_sweep_cue/<net>_pool<L>.mat
        res.net, res.layer, res.population, res.N, res.n_silent
        res.k_list [nK x 1]; res.pair_labels, res.pair_group {6 x 1}
        [nK x 6 pairs x 10 folds]: self, nopt, pt_rot, pt_full, d_rot, d_full
        [nK x 6 x 10 x n_perm]: ctrl, ctrl_id
        res.pr_raw, res.pr_cv [1 x 3 cues]; res.cumvar, res.cv_cumvar {1 x 3}
            (variance captured by the top k PCs of cue c = res.cumvar{c}(k))
    results/decoding_outputs/ann_pc_sweep_hierarchy/<net>_pool<lo>_pool<hi>.mat
        res.net, res.layers [lo hi], res.layer_names, res.N_layer [1 x 2]
        res.dir_labels; res.combo_source_cue, res.combo_target_cue {9 x 1}; res.done [2 x 1]
        [2 directions x 9 combos x nK x 10]: self, pt_rot, pt_full, d_rot, d_full (nopt = NaN)
        [2 x 9 x nK x 10 x n_perm]: ctrl, ctrl_id
        res.pr_raw, res.pr_cv [2 layers x 3 cues]
    Each file also holds res.summary: one row per k (pooled, and per pair /
    direction) with ratios of pooled means (transfer_metrics.m).
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
cfg.comparisons = {'cue', 'hierarchy'};                                 % default: both
cfg.layers      = {'alexnet',1; 'alexnet',2; 'alexnet',5; 'vgg',1; 'vgg',2; 'vgg',3; 'vgg',4; 'vgg',5};   % for 'cue'
cfg.layer_pairs = {'alexnet', 2, 5};                                    % for 'hierarchy': {net, lower, higher}
cfg.k_base      = [1:10, 15:5:45, 49];
cfg.n_perm      = 10;
cfg.pr_splits   = 20;
cfg.pairs       = {'ac','ec'; 'ec','ac'; 'ec','ex'; 'ex','ec'; 'ac','ex'; 'ex','ac'};   % {source, target}
cfg.base_seed   = 27000;
cfg.data_root   = '';               % '' = ANN_decoding/ANN_response_data (load_ann_layer default)
cfg.out_root    = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag_prefix  = 'ann_pc_sweep';   % folders <tag_prefix>_cue / _hierarchy
cfg.overwrite   = false;

if exist('cfg_override', 'var') && isstruct(cfg_override)
    fo  = fieldnames(cfg_override);
    bad = setdiff(fo, fieldnames(cfg));
    assert(isempty(bad), 'cfg_override has unknown field(s): %s', strjoin(bad', ', '));
    for i = 1:numel(fo), cfg.(fo{i}) = cfg_override.(fo{i}); end
    fprintf('cfg_override applied: %s\n', strjoin(fo', ', '));
end
bad = setdiff(cfg.comparisons, {'cue', 'hierarchy'});
assert(isempty(bad), 'Unknown comparison(s): %s', strjoin(bad, ', '));

if isempty(gcp('nocreate')), parpool('Processes'); end
t_all = tic;

%% Run
for ic = 1:numel(cfg.comparisons)
    out_dir = fullfile(cfg.out_root, sprintf('%s_%s', cfg.tag_prefix, cfg.comparisons{ic}));
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    switch cfg.comparisons{ic}
        case 'cue'
            for il = 1:size(cfg.layers, 1)
                run_cue(cfg, lower(cfg.layers{il, 1}), cfg.layers{il, 2}, out_dir);
            end
        case 'hierarchy'
            for ipair = 1:size(cfg.layer_pairs, 1)
                run_hierarchy(cfg, lower(cfg.layer_pairs{ipair, 1}), cfg.layer_pairs{ipair, 2}, cfg.layer_pairs{ipair, 3}, out_dir);
            end
    end
end
fprintf('\nDone in %.1f min.\n', toc(t_all) / 60);

%% ===========================================================================
function run_cue(cfg, net, layer, out_dir)
% Cue transfer on the whole layer, every k, six ordered pairs.
lay_id = ann_layer_id(net, layer);
pop = load_ann_layer(net, layer, struct('data_root', cfg.data_root));
F = pop.n_rep;  M = cfg.n_perm;  nP = size(cfg.pairs, 1);
k_list = k_grid(cfg.k_base, min([pop.n_stim - 1, pop.N]));
nK = numel(k_list);
src = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 1));
tgt = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 2));

file = fullfile(out_dir, sprintf('%s_pool%d.mat', net, layer));
[~, data_name] = fileparts(pop.file);
sig  = struct('k_list', k_list, 'n_perm', M, 'pr_splits', cfg.pr_splits, 'pairs', {cfg.pairs}, ...
    'base_seed', cfg.base_seed, 'data_file', data_name);
res = open_results(file, sig, cfg.overwrite);
if ~isempty(res) && all(res.done)
    fprintf('%s cue PC sweep already done, skipped.\n', pop.name);
    return;
end
res = struct();
res.info        = 'CNN PC sweep, cue transfer, whole layer. Field definitions: header of run_ann_pc_sweep.m.';
res.population  = pop.name;  res.net = net;  res.layer = layer;  res.N = pop.N;
res.n_silent    = pop.n_silent;
res.sig         = sig;  res.cfg = cfg;  res.created = char(datetime('now'));
res.k_list      = k_list(:);
res.pair_labels = strcat(cfg.pairs(:, 1), '->', cfg.pairs(:, 2));
res.pair_group  = pair_group_names(cfg.pairs);
res.cue_names   = pop.cue_names;
res.done        = false;
for c = 1:3
    pr = population_pr(zscore(pop.data{c}), pop.labels, ...
        struct('n_split', cfg.pr_splits, 'seed', cfg.base_seed + lay_id * 1e7 + 10 + c));
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
    seeds(ii) = cfg.base_seed + lay_id * 1e7 + p;        % independent of k
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
print_k5(res, pop.name);
end

%% ===========================================================================
function run_hierarchy(cfg, net, lo, hi, out_dir)
% Whole lower layer vs whole higher layer, every k, both directions, 9 rendering combinations.
pair_id = ann_layer_id(net, lo) * 10 + ann_layer_id(net, hi);
pops = {load_ann_layer(net, lo, struct('data_root', cfg.data_root)), ...
        load_ann_layer(net, hi, struct('data_root', cfg.data_root))};
assert(isequal(pops{1}.labels, pops{2}.labels), 'The two layers have different trial layouts.');
F = pops{1}.n_rep;  M = cfg.n_perm;
cues = pops{1}.cue_names;
[cs, ct] = ndgrid(1:3, 1:3);  cs = cs(:);  ct = ct(:);  nC = numel(cs);
k_list = k_grid(cfg.k_base, min([pops{1}.n_stim - 1, pops{1}.N, pops{2}.N]));
nK = numel(k_list);
names = {sprintf('pool%d', lo), sprintf('pool%d', hi)};
dir_labels = {sprintf('%s-to-%s (train %s, test %s)', names{1}, names{2}, names{1}, names{2}); ...
              sprintf('%s-to-%s (train %s, test %s)', names{2}, names{1}, names{2}, names{1})};
dir_target = [1; 2];  dir_source = [2; 1];

file = fullfile(out_dir, sprintf('%s_pool%d_pool%d.mat', net, lo, hi));
sig  = struct('k_list', k_list, 'n_perm', M, 'pr_splits', cfg.pr_splits, 'base_seed', cfg.base_seed);
res = open_results(file, sig, cfg.overwrite);
if isempty(res)
    res = struct();
    res.info        = 'CNN PC sweep, layer-to-layer, whole layers. Field definitions: header of run_ann_pc_sweep.m.';
    res.net         = net;  res.layers = [lo hi];
    res.layer_names = {pops{1}.name; pops{2}.name};
    res.N_layer     = [pops{1}.N, pops{2}.N];
    res.n_silent    = [pops{1}.n_silent; pops{2}.n_silent];
    res.sig         = sig;  res.cfg = cfg;  res.created = char(datetime('now'));
    res.k_list      = k_list(:);
    res.dir_labels  = dir_labels;
    res.dir_source_layer = names([2 1])';  res.dir_target_layer = names([1 2])';
    res.combo_source_cue = cues(cs)';  res.combo_target_cue = cues(ct)';
    res.cue_names   = cues;
    res.done        = false(2, 1);
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
                struct('n_split', cfg.pr_splits, 'seed', cfg.base_seed + pair_id * 1e6 + 10 + 3 * a + c));
            res.pr_raw(a, c) = pr.raw;  res.pr_cv(a, c) = pr.cv;
        end
    end
end
fprintf('\n===== %s layer-to-layer PC sweep: %s %d vs %s %d units, k = %s =====\n', ...
    net, names{1}, pops{1}.N, names{2}, pops{2}.N, mat2str(k_list));

labels = pops{1}.labels;
for di = 1:2
    if res.done(di), continue; end
    t0 = tic;
    Xsrc = pops{dir_source(di)}.data;  Xtgt = pops{dir_target(di)}.data;
    nI = nK * nC;
    seeds = zeros(nI, 1);
    for ii = 1:nI
        [~, c] = ind2sub([nK, nC], ii);
        seeds(ii) = cfg.base_seed + pair_id * 1e6 + di * 1e4 + c;   % independent of k
    end
    outs = cell(nI, 1);
    fprintf('[%s hierarchy] direction %d/2: %s, %d k values x %d rendering combinations = %d calls\n', ...
        net, di, dir_labels{di}, nK, nC, nI);
    q = progress_queue(sprintf('%s %s', net, dir_labels{di}), nI);
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
    res.summary = hierarchy_summary(res);
    save_atomic(file, res);
    ik5 = find(k_list == min(5, max(k_list)));
    m = transfer_metrics(res.self(di, :, ik5, :), res.pt_rot(di, :, ik5, :), res.ctrl(di, :, ik5, :, :));
    fprintf('[%s %s] k = %d: self %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f | %.0f s\n', ...
        net, dir_labels{di}, k_list(ik5), m.self, m.pt, m.ctrl, m.ratio, m.headroom, res.runtime_s(di));
end
end

%% ===========================================================================
function id = ann_layer_id(net, layer)
% fixed index of a layer, so seeds do not depend on which layers a run includes
all_layers = {'alexnet',1; 'alexnet',2; 'alexnet',5; 'vgg',1; 'vgg',2; 'vgg',3; 'vgg',4; 'vgg',5};
id = find(strcmp(all_layers(:, 1), net) & [all_layers{:, 2}]' == layer);
assert(~isempty(id), 'Unknown layer %s pool%d', net, layer);
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
% cue sweep: k is the first dimension, the pair the second
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

function T = hierarchy_summary(res)
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

function print_k5(res, name)
ik = find(res.k_list == 5);
if isempty(ik), return; end
m = transfer_metrics(res.self(ik, :, :), res.pt_rot(ik, :, :), res.ctrl(ik, :, :, :), res.nopt(ik, :, :));
fprintf('[%s cue] k = 5: self %.3f  nonPT %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f\n', ...
    name, m.self, m.nopt, m.pt, m.ctrl, m.ratio, m.headroom);
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
