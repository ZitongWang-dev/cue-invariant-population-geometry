%{
Filename: run_across_area.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Across-area (V1 <-> V2) transfer decoding in the full neuron space, for both
    monkeys (default), both directions and all 3 x 3 rendering combinations.
    Feeds Fig 4b (per training rendering, averaged over the three test
    renderings) and Fig 12 (all nine combinations).

    Directions, named as in the paper:
      'V1-to-V2'  decoder trained on V1 (target = V1), tested on V2 (source = V2)
      'V2-to-V1'  decoder trained on V2 (target = V2), tested on V1 (source = V1)
    Both areas are matched to the smaller one (V1: FR 48, KO 109): the smaller
    area is used whole and the larger one is subsampled, cfg.n_draws times. The
    same subsets serve both directions and all nine rendering combinations of a
    draw. V1 and V2 were recorded in separate sessions, but pt_transfer_core.m
    holds out the same rows on both sides anyway. The non-PT column maps
    different neurons onto each other and is only a sanity check (~chance).

Run:
    By default it runs both monkeys. To split the work across machines:
        cfg_override = struct('monkeys', {{'KO'}});
        run_across_area
    Re-running the same configuration resumes (a finished direction is skipped).
    Clear cfg_override before running another driver.

Output: results/decoding_outputs/across_area/<monkey>.mat, variable `res`:
      res.dir_labels {2 x 1}, res.dir_source_area, res.dir_target_area
      res.combo_source_cue, res.combo_target_cue {9 x 1}   (source = tested, target = trained)
      res.idx_V1, res.idx_V2 {1 x n_draws}   neurons used per draw
      [2 directions x 9 combos x n_draws x 10 folds]:
        self (target area), nopt (sanity, different neurons), pt_rot, pt_full, d_rot, d_full
      [2 x 9 x n_draws x 10 x n_perm]: ctrl, ctrl_id
      res.pr_raw, res.pr_cv [2 areas (V1, V2) x 3 cues x n_draws]
      res.done [2 x 1], res.runtime_s [2 x 1]
      res.summary  table: per direction, pooled / same rendering / per training
                   rendering (Fig 4b bars) / each of the nine combinations (Fig 12)
%}

%% Initialization
clearvars -except cfg_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
addpath(script_dir);

%% Configuration
cfg = struct();
cfg.monkeys    = {'FR', 'KO'};      % default: both
cfg.timewindow = [330 630];
cfg.n_draws    = 20;                % subsets of the larger area
cfg.n_perm     = 10;
cfg.pr_splits  = 20;
cfg.base_seed  = 4000;
cfg.out_root   = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag        = 'across_area';
cfg.overwrite  = false;

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

all_monkeys = {'FR', 'KO'};
dir_labels  = {'V1-to-V2 (train V1, test V2)'; 'V2-to-V1 (train V2, test V1)'};
dir_target  = [1; 2];               % area index (1 = V1, 2 = V2) the decoder is trained on
dir_source  = [2; 1];               % area index whose trials are transformed and tested
t_all = tic;

%% Run
for im = 1:numel(cfg.monkeys)
    monkey = cfg.monkeys{im};
    mk_id  = find(strcmp(all_monkeys, monkey));
    assert(~isempty(mk_id), 'Unknown monkey %s', monkey);

    pops = {load_population(monkey, 'V1', struct('timewindow', cfg.timewindow)), ...
            load_population(monkey, 'V2', struct('timewindow', cfg.timewindow))};
    assert(isequal(pops{1}.labels, pops{2}.labels), 'V1 and V2 trial layouts differ.');
    Nar = [pops{1}.N, pops{2}.N];
    n   = min(Nar);
    F   = pops{1}.n_rep;
    M   = cfg.n_perm;
    D   = cfg.n_draws;
    cues = pops{1}.cue_names;
    [cs, ct] = ndgrid(1:3, 1:3);    % combo c: source cue cs(c), target cue ct(c)
    cs = cs(:);  ct = ct(:);
    nC = numel(cs);

    file = fullfile(out_dir, sprintf('%s.mat', monkey));
    sig  = struct('timewindow', cfg.timewindow, 'n_draws', D, 'n_perm', M, ...
        'pr_splits', cfg.pr_splits, 'base_seed', cfg.base_seed, 'N_match', n);
    res = open_results(file, sig, cfg.overwrite);
    if isempty(res)
        res = struct();
        res.info       = 'Across-area decoding, full space. Field definitions: header of run_across_area.m.';
        res.monkey     = monkey;
        res.N_area     = Nar;          % [V1 V2] recorded neurons
        res.N_match    = n;
        res.sig        = sig;
        res.cfg        = cfg;
        res.created    = char(datetime('now'));
        res.dir_labels = dir_labels;
        res.dir_source_area = {'V2'; 'V1'};
        res.dir_target_area = {'V1'; 'V2'};
        res.combo_source_cue = cues(cs)';
        res.combo_target_cue = cues(ct)';
        res.cue_names  = cues;
        res.done       = false(2, 1);
        for fn = {'self', 'nopt', 'pt_rot', 'pt_full', 'd_rot', 'd_full'}
            res.(fn{1}) = nan(2, nC, D, F);
        end
        res.ctrl       = nan(2, nC, D, F, M);
        res.ctrl_id    = nan(2, nC, D, F, M);
        res.runtime_s  = nan(2, 1);
        res.summary    = table();

        % neuron subsets per draw (shared by both directions) and their PR
        res.idx_V1 = cell(1, D);
        res.idx_V2 = cell(1, D);
        res.pr_raw = nan(2, 3, D);
        res.pr_cv  = nan(2, 3, D);
        for d = 1:D
            rsn = RandStream('twister', 'Seed', seed_of(cfg.base_seed, mk_id, 0, d, 0));
            idx = cell(1, 2);
            for a = 1:2
                if Nar(a) == n, idx{a} = 1:n; else, idx{a} = sort(randperm(rsn, Nar(a), n)); end
                for c = 1:3
                    pr = population_pr(zscore(pops{a}.data{c}(:, idx{a})), pops{a}.labels, ...
                        struct('n_split', cfg.pr_splits, 'seed', seed_of(cfg.base_seed, mk_id, 0, d, 10 + 3 * a + c)));
                    res.pr_raw(a, c, d) = pr.raw;
                    res.pr_cv(a, c, d)  = pr.cv;
                end
            end
            res.idx_V1{d} = idx{1};
            res.idx_V2{d} = idx{2};
        end
    end
    fprintf('\n===== %s across areas: V1 %d / V2 %d neurons, matched to %d | %d draws, %d folds, %d permutations =====\n', ...
        monkey, Nar(1), Nar(2), n, D, F, M);

    data = {pops{1}.data, pops{2}.data};
    labels = pops{1}.labels;
    idxs = {res.idx_V1, res.idx_V2};
    for di = 1:2
        if res.done(di), continue; end
        t0 = tic;
        sa = dir_source(di);  ta = dir_target(di);
        nI = D * nC;
        seeds = zeros(nI, 1);
        for ii = 1:nI
            [d, c] = ind2sub([D, nC], ii);
            seeds(ii) = seed_of(cfg.base_seed, mk_id, di, d, c);
        end
        Xsrc = data{sa};  Xtgt = data{ta};
        Isrc = idxs{sa};  Itgt = idxs{ta};
        outs = cell(nI, 1);
        fprintf('[%s] direction %d/2: %s, %d draws x %d rendering combinations = %d calls\n', monkey, di, dir_labels{di}, D, nC, nI);
        q = progress_queue(sprintf('%s %s', monkey, dir_labels{di}), nI);
        parfor ii = 1:nI
            [d, c] = ind2sub([D, nC], ii);
            o = struct('n_perm', M, 'seed', seeds(ii), 'same_neurons', false);
            outs{ii} = pt_transfer_core(Xsrc{cs(c)}(:, Isrc{d}), Xtgt{ct(c)}(:, Itgt{d}), labels, o);
            send(q, 1);
        end

        for ii = 1:nI
            [d, c] = ind2sub([D, nC], ii);
            o = outs{ii};
            res.self(di, c, d, :)    = o.self;
            res.nopt(di, c, d, :)    = o.nopt;
            res.pt_rot(di, c, d, :)  = o.pt_rot;
            res.pt_full(di, c, d, :) = o.pt_full;
            res.d_rot(di, c, d, :)   = o.d_rot;
            res.d_full(di, c, d, :)  = o.d_full;
            res.ctrl(di, c, d, :, :)    = reshape(o.ctrl,    [1 1 1 F M]);
            res.ctrl_id(di, c, d, :, :) = reshape(o.ctrl_id, [1 1 1 F M]);
        end
        res.done(di)      = true;
        res.runtime_s(di) = toc(t0);
        res.updated       = char(datetime('now'));
        res.summary       = build_summary(res);
        save_atomic(file, res);

        m = transfer_metrics(res.self(di, :, :, :), res.pt_rot(di, :, :, :), res.ctrl(di, :, :, :, :), res.nopt(di, :, :, :));
        fprintf('[%s %s] all 9 combos: self %.3f  nonPT(sanity) %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.2f | %.0f s\n', ...
            monkey, dir_labels{di}, m.self, m.nopt, m.pt, m.ctrl, m.ratio, m.headroom, res.runtime_s(di));
    end
end
fprintf('\nDone in %.1f min. Results in %s\n', toc(t_all) / 60, out_dir);

%% ---------------------------------------------------------------------------
function s = seed_of(base, mk_id, di, d, item)
s = base + mk_id * 1e7 + di * 1e5 + d * 100 + item;
end

function T = build_summary(res)
rows = {};
nC = numel(res.combo_source_cue);
for di = find(res.done)'
    sets  = {1:nC, find(strcmp(res.combo_source_cue, res.combo_target_cue))'};
    names = {'all 9 combos', 'same rendering'};
    for c = 1:numel(res.cue_names)        % Fig 4b: one training rendering, all three test renderings
        sets{end + 1}  = find(strcmp(res.combo_target_cue, res.cue_names{c}))'; %#ok<AGROW>
        names{end + 1} = sprintf('train %s -> test AC|EC|EX', upper(res.cue_names{c})); %#ok<AGROW>
    end
    for c = 1:nC                          % Fig 12: every combination
        sets{end + 1}  = c; %#ok<AGROW>
        names{end + 1} = sprintf('train %s -> test %s', upper(res.combo_target_cue{c}), upper(res.combo_source_cue{c})); %#ok<AGROW>
    end
    for g = 1:numel(sets)
        C = sets{g};
        m = transfer_metrics(res.self(di, C, :, :), res.pt_rot(di, C, :, :), res.ctrl(di, C, :, :, :), res.nopt(di, C, :, :));
        rows(end + 1, :) = {res.dir_labels{di}, names{g}, m.self, m.nopt, m.pt, ...
            mean(res.pt_full(di, C, :, :), 'all', 'omitnan'), m.ctrl, ...
            mean(res.ctrl_id(di, C, :, :, :), 'all', 'omitnan'), m.ratio, m.diff, m.headroom}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'direction', 'set', 'self', 'nopt_sanity', 'pt_rot', 'pt_full', ...
    'ctrl', 'ctrl_id', 'pt_over_ctrl', 'pt_minus_ctrl', 'headroom'});
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
    fprintf('Resuming %s (%d of %d directions done)\n', file, nnz(res.done), numel(res.done));
end
end

function save_atomic(file, res)
tmp = [file(1:end - 4) '_tmp.mat'];
save(tmp, 'res', '-v7');
movefile(tmp, file, 'f');
end
