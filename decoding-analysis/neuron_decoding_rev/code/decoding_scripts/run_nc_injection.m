%{
Filename: run_nc_injection.m
Author:   Zitong Wang
Date:     2026-10-24

Description:
    Reviewer 2, major point 1 (rebuttal Point 8): inject structured noise
    correlations (NC) of realistic magnitude (mean pairwise r = 0.05-0.15)
    into the recorded pseudo-populations and test whether Procrustes-transfer
    (PT) decoding, scored against the corrected control, is robust to them.

    Conditions, for each whole population and all six ordered cue pairs:
      original      the recorded pseudo-population (real within-session NC)
      affine        NC removed by the deterministic affine trial shift of
                    Point 7 (noise_correlation_analysis_affine.m): the clean start
      inject <r>    NC injected into the affine-shifted data, with mean pairwise
                    NC r = cfg.levels (default 0, 0.05, 0.10, 0.15)

    Injection. One correlation matrix per population and seed, built from a
    few shared fluctuations (a factor model, a valid correlation matrix by
    construction):
        R = L*L' + D,   L = [a*1, b*B],   B (N x k) drawn from N(0, 1),
        D = diag(1 - a^2 - b^2*|B_i|^2)   (private variance)
    Off-diagonal entries are then roughly Gaussian with mean a^2 and SD
    b^2*sqrt(k). a sets the level; b and B are the same at every level, so the
    levels are paired. For every rendering and stimulus, 10 "trials" are drawn
    from N(0, R), and each neuron's 10 real spike counts are reordered so that
    their ranks follow its column (Iman-Conover). Every neuron keeps exactly its
    own counts per stimulus, so means, variances, z-scoring and the
    trial-averaged geometry are unchanged; only how trials line up across
    neurons changes. The result does not depend on the starting trial order,
    so no NC of the original data carries over. The same R is used for all
    stimuli and renderings.

    Target spread (SD of the pairwise NC across pairs): measured from the
    same-session pairs of all four recorded populations, as the SD of each
    pair's mean NC minus the estimation noise of that mean (cfg.nc_sd = []),
    or given as a number. It is kept fixed across levels.

    Calibration. Ten trials and integer counts weaken the achieved correlation,
    so a^2 and b^2 are rescaled (cfg.calib_iter rounds on cfg.calib_seeds
    separate draws) until the achieved mean at the top level and the achieved
    spread, measured with the NC metric below, match the targets. The achieved
    NC of every condition and seed is stored.

    NC metric (as Point 7): Pearson correlation of two neurons' spike counts
    over the 10 repeats of one stimulus; per pair, the mean over the 50 stimuli
    x 3 renderings (NaN, i.e. a silent neuron, excluded). Noise PR (check of the
    noise dimensionality): participation ratio of the correlation matrix of the
    residuals (trial minus its stimulus mean), mean over renderings.

    Decoding: pt_transfer_core.m unchanged (z-score per side, 10-fold row
    partition, rotation-only PT, corrected control with cfg.n_perm permutations
    per fold). The partition and permutations depend on (population, seed, pair)
    only, so all conditions of a seed are decoded on the same folds.

Run:
    From MATLAB, run this script (any current folder). It runs all four
    populations by default: 6 conditions x 6 pairs x 10 seeds = 360 core calls
    per population, ~2.2 h for all four on an 8-worker machine. Quick check of
    a machine first (~5 min, separate output folder):
        cfg_override = struct('populations', {{'FR','V1'}}, 'n_seeds', 1, 'n_perm', 2, ...
            'pairs', {{'ac','ec'}}, 'tag', 'nc_injection_test');
        run_nc_injection
    Any cfg field below can be overridden this way (unknown fields are refused).
    Re-running the same configuration resumes at the first unfinished seed.
    Clear cfg_override (clear cfg_override) before running another driver.
    Figures: visualization_scripts/plot_nc_injection.m.

Output: results/decoding_outputs/<tag>/<monkey>_<area>.mat, one variable res:
    res.cond_names {nC x 1}, res.cond_target [nC x 1] (NaN for original, affine)
    res.pair_labels {nP x 1} 'ac->ec' = source ac (transformed, tested), target ec
    res.pair_group  {nP x 1} unordered pair: 'EC-AC', 'EC-EX', 'AC-EX'
    res.done [nS x 1]
    [nC x nP x nS x 10 folds]:          self, nopt, pt_rot, pt_full
    [nC x nP x nS x 10 x n_perm]:       ctrl (the floor), ctrl_id (~1/50)
    [nC x nS] achieved NC:              nc_mean, nc_sd_raw, nc_sd_true,
                                        nc_mean_same, nc_mean_cross, noise_pr
    res.nc_mean_by_cue [nC x nS x 3]
    res.pair_nc_seed1 [n_pairs x nC]    per-pair mean NC of seed 1 (distributions)
    res.is_same_pair [n_pairs x 1]      pair recorded in the same session
    res.session_id [N x 1]
    res.target       struct: mean levels, sd, k, the measurement behind sd
    res.calib        struct: a2 per level, b2, calibration history
    res.n_capped [nS x 1]               neurons whose shared variance was capped
    res.runtime_s [nS x 1]
    res.summary      table, one row per condition x {pooled, EC-AC, EC-EX, AC-EX}
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
cfg.levels      = [0 0.05 0.10 0.15];   % target mean pairwise NC of the injected conditions
cfg.nc_sd       = [];               % target SD across pairs; [] = measured from same-session pairs
cfg.k_factors   = 5;                % random shared fluctuations (the structure of the noise)
cfg.max_shared  = 0.95;             % cap of a^2 + b^2*|B_i|^2 (private variance >= 0.05)
cfg.n_seeds     = 10;               % injection draws (each with its own fold partition)
cfg.n_perm      = 10;               % control permutations per fold
cfg.calib_iter  = 3;                % calibration rounds
cfg.calib_seeds = 2;                % draws per calibration round (separate from the run seeds)
cfg.pairs       = {'ac','ec'; 'ec','ac'; 'ec','ex'; 'ex','ec'; 'ac','ex'; 'ex','ac'};   % {source, target}
cfg.base_seed   = 3000;
cfg.out_root    = fullfile(script_dir, '..', '..', 'results', 'decoding_outputs');
cfg.tag         = 'nc_injection';
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
nS = cfg.n_seeds;
M  = cfg.n_perm;
k  = cfg.k_factors;
cond_names  = [{'original'; 'affine'}; arrayfun(@(r) sprintf('inject %.2f', r), cfg.levels(:), 'UniformOutput', false)];
cond_target = [NaN; NaN; cfg.levels(:)];
nC = numel(cond_names);
nL = numel(cfg.levels);
t_all = tic;

%% Target spread: same-session pairs of all four populations (the same on every machine)
if isempty(cfg.nc_sd)
    [nc_sd, sd_info] = measure_target_sd(all_pops, cfg.timewindow);
    fprintf(['Target NC spread from %d same-session pairs (4 populations): mean %.3f, SD of pair means %.3f, ' ...
        'estimation noise %.3f -> SD %.3f\n'], sd_info.n_pairs, sd_info.mean, sd_info.sd_raw, sd_info.est_sd, nc_sd);
else
    nc_sd = cfg.nc_sd;
    sd_info = struct('source', 'cfg.nc_sd');
    fprintf('Target NC spread from cfg.nc_sd: %.3f\n', nc_sd);
end

%% Run
for ip = 1:size(cfg.populations, 1)
    monkey = cfg.populations{ip, 1};
    area   = cfg.populations{ip, 2};
    pop_id = find(strcmp(all_pops(:, 1), monkey) & strcmp(all_pops(:, 2), area));
    assert(~isempty(pop_id), 'Unknown population %s %s', monkey, area);

    pop = load_population(monkey, area, struct('timewindow', cfg.timewindow));
    N = pop.N;
    F = pop.n_rep;
    labels = pop.labels;
    src = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 1));
    tgt = cellfun(@(c) find(strcmp(pop.cue_names, c)), cfg.pairs(:, 2));
    sess = session_ids(pop);
    Xa = cellfun(@(X) affine_shift(X, labels, sess), pop.data, 'UniformOutput', false);

    % results file: resume it or start a new one
    file = fullfile(out_dir, sprintf('%s_%s.mat', monkey, area));
    sig  = struct('timewindow', cfg.timewindow, 'levels', cfg.levels, 'nc_sd', nc_sd, 'k_factors', k, ...
        'max_shared', cfg.max_shared, 'n_seeds', nS, 'n_perm', M, 'calib_iter', cfg.calib_iter, ...
        'calib_seeds', cfg.calib_seeds, 'pairs', {cfg.pairs}, 'base_seed', cfg.base_seed);
    res = open_results(file, sig, cfg.overwrite);
    if isempty(res)
        fprintf('\n===== %s: calibrating the injection =====\n', pop.name);
        calib = calibrate_injection(Xa, labels, sess, cfg, nc_sd, pop_id);
        [~, is_same] = pair_nc(pop.data, labels, sess);

        res = struct();
        res.info        = 'NC injection into the recorded pseudo-population. Field definitions: header of run_nc_injection.m.';
        res.population  = pop.name;  res.monkey = monkey;  res.area = area;  res.N_total = N;
        res.sig         = sig;
        res.cfg         = cfg;
        res.created     = char(datetime('now'));
        res.cond_names  = cond_names;
        res.cond_target = cond_target;
        res.pair_labels = strcat(cfg.pairs(:, 1), '->', cfg.pairs(:, 2));
        res.pair_group  = pair_group_names(cfg.pairs);
        res.cue_names   = pop.cue_names;
        res.session_id  = sess;
        res.is_same_pair = is_same;
        res.target      = struct('levels', cfg.levels, 'sd', nc_sd, 'k', k, 'sd_measurement', sd_info);
        res.calib       = calib;
        res.done        = false(nS, 1);
        for fn = {'self', 'nopt', 'pt_rot', 'pt_full'}
            res.(fn{1}) = nan(nC, nP, nS, F);
        end
        res.ctrl        = nan(nC, nP, nS, F, M);
        res.ctrl_id     = nan(nC, nP, nS, F, M);
        for fn = {'nc_mean', 'nc_sd_raw', 'nc_sd_true', 'nc_mean_same', 'nc_mean_cross', 'noise_pr'}
            res.(fn{1}) = nan(nC, nS);
        end
        res.nc_mean_by_cue = nan(nC, nS, 3);
        res.pair_nc_seed1  = nan(numel(is_same), nC);
        res.n_capped    = nan(nS, 1);
        res.runtime_s   = nan(nS, 1);
        res.summary     = table();
        save_atomic(file, res);
    end
    a2 = res.calib.a2;  b2 = res.calib.b2;
    fprintf('\n===== %s: %d neurons, %d sessions | %d conditions x %d pairs x %d seeds | a^2 = %s, b^2 = %.4f =====\n', ...
        pop.name, N, max(sess), nC, nP, nS, mat2str(round(a2, 4)), b2);
    seeds_left = nnz(~res.done);
    seeds_run  = 0;
    t_pop = tic;

    for is = 1:nS
        if res.done(is), continue; end
        t0 = tic;

        % data of every condition for this seed (same draws at every level)
        dr = make_draws(seed_of(cfg.base_seed, pop_id, is, 50), N, k, pop.n_stim, F, numel(pop.data));
        [B, n_cap] = cap_loadings(dr.B, max(a2), b2, cfg.max_shared);
        Xc = cell(nC, numel(pop.data));
        Xc(1, :) = pop.data;
        Xc(2, :) = Xa;
        for l = 1:nL
            Xc(2 + l, :) = inject_nc(Xa, labels, a2(l), b2, B, dr);
        end

        % achieved noise correlation and noise dimensionality (manipulation check)
        for ic = 1:nC
            st = nc_summary(Xc(ic, :), labels, sess, true);
            res.nc_mean(ic, is)       = st.mean_all;
            res.nc_sd_raw(ic, is)     = st.sd_raw;
            res.nc_sd_true(ic, is)    = st.sd_true;
            res.nc_mean_same(ic, is)  = st.mean_same;
            res.nc_mean_cross(ic, is) = st.mean_cross;
            res.noise_pr(ic, is)      = st.noise_pr;
            res.nc_mean_by_cue(ic, is, :) = st.mean_by_cue;
            if is == 1, res.pair_nc_seed1(:, ic) = st.pair_mean; end
        end

        % decode every (condition, pair) in parallel; folds depend on (seed, pair) only
        nI = nC * nP;
        seeds = arrayfun(@(p) seed_of(cfg.base_seed, pop_id, is, p), 1:nP);
        outs = cell(nI, 1);
        fprintf('[%s] seed %d/%d: %d conditions x %d pairs = %d calls\n', pop.name, is, nS, nC, nP, nI);
        q = progress_queue(sprintf('%s seed %d', pop.name, is), nI);
        parfor ii = 1:nI
            [ic, p] = ind2sub([nC, nP], ii);
            o = struct('n_perm', M, 'seed', seeds(p), 'same_neurons', true);
            outs{ii} = pt_transfer_core(Xc{ic, src(p)}, Xc{ic, tgt(p)}, labels, o);
            send(q, 1);
        end

        for ii = 1:nI
            [ic, p] = ind2sub([nC, nP], ii);
            o = outs{ii};
            res.self(ic, p, is, :)    = o.self;
            res.nopt(ic, p, is, :)    = o.nopt;
            res.pt_rot(ic, p, is, :)  = o.pt_rot;
            res.pt_full(ic, p, is, :) = o.pt_full;
            res.ctrl(ic, p, is, :, :)    = reshape(o.ctrl,    [1 1 1 F M]);
            res.ctrl_id(ic, p, is, :, :) = reshape(o.ctrl_id, [1 1 1 F M]);
        end
        res.n_capped(is)  = n_cap;
        res.done(is)      = true;
        res.runtime_s(is) = toc(t0);
        res.updated       = char(datetime('now'));
        res.summary       = build_summary(res);
        save_atomic(file, res);

        seeds_run = seeds_run + 1;
        fprintf('[%s] seed %d done in %.0f s (%d neurons capped) | ~%.1f h left for %s\n', pop.name, is, ...
            res.runtime_s(is), n_cap, toc(t_pop) / seeds_run * (seeds_left - seeds_run) / 3600, pop.name);
        for ic = 1:nC
            m = transfer_metrics(res.self(ic, :, is, :), res.pt_rot(ic, :, is, :), res.ctrl(ic, :, is, :, :));
            fprintf('    %-12s NC %6.3f (SD %.3f, same %.3f, cross %.3f) noise PR %5.1f | self %.3f  PT %.3f  ctrl %.3f  PT/ctrl %.2f  headroom %.3f\n', ...
                cond_names{ic}, res.nc_mean(ic, is), res.nc_sd_true(ic, is), res.nc_mean_same(ic, is), ...
                res.nc_mean_cross(ic, is), res.noise_pr(ic, is), m.self, m.pt, m.ctrl, m.ratio, m.headroom);
        end
    end
end
fprintf('\nDone in %.1f min. Results in %s\n', toc(t_all) / 60, out_dir);

%% ---------------------------------------------------------------------------
function s = seed_of(base, pop_id, is, item)
% unique seed per (population, seed index, item); item < 100, is < 1000
s = base + pop_id * 1e7 + is * 1e3 + item;
end

function sess = session_ids(pop)
% Recording session of every neuron, in pseudo-population order (as
% build_affine_shift_config.m: sessions sorted by file name, good channels only).
sdir  = fullfile(fileparts(pop.file), 'ac');
files = dir(fullfile(sdir, '*.mat'));
names = sort({files.name});
names = names(~strcmpi(names, sprintf('%s_%s_ac.mat', pop.monkey, pop.area)));
sess = [];
for i = 1:numel(names)
    L = load(fullfile(sdir, names{i}), 'bad_channel');
    n_good = numel(setdiff(1:24, L.bad_channel));
    sess = [sess; repmat(i, n_good, 1)]; %#ok<AGROW>
end
assert(numel(sess) == pop.N, 'session_ids: %d neurons in the session files, %d in %s.', numel(sess), pop.N, pop.name);
end

function Y = affine_shift(X, labels, sess)
% Deterministic affine trial shift of Point 7 (noise_correlation_analysis_affine.m):
% neuron j of a session gets permutation j of t -> (a*t + b) mod 10, so neurons
% recorded together never keep their trial alignment.
n_rep = sum(labels == 1);
assert(n_rep == 10, 'affine_shift: the 24 patterns are defined for 10 repeats.');
pat = zeros(24, n_rep);
ab = [ones(10, 1), (0:9)'; 9 * ones(10, 1), (0:9)'; 3 * ones(4, 1), (0:3)'];
for i = 1:24
    pat(i, :) = mod(ab(i, 1) * (0:n_rep - 1) + ab(i, 2), n_rep) + 1;
end
pid = zeros(numel(sess), 1);
for s = unique(sess)'
    idx = find(sess == s);
    assert(numel(idx) <= 24, 'affine_shift: session %d has more than 24 neurons.', s);
    pid(idx) = 1:numel(idx);
end
Y = X;
for st = 1:max(labels)
    rows = find(labels == st);
    for n = 1:size(X, 2)
        Y(rows, n) = X(rows(pat(pid(n), :)), n);
    end
end
end

function dr = make_draws(seed, N, k, n_stim, n_rep, n_cue)
% Gaussian draws of one seed: loadings B, shared factors F, private noise E
rs = RandStream('twister', 'Seed', seed);
dr.B = randn(rs, N, k);
dr.F = randn(rs, n_rep, 1 + k, n_stim, n_cue);
dr.E = randn(rs, n_rep, N, n_stim, n_cue);
end

function [B, n_cap] = cap_loadings(B, a2max, b2, max_shared)
% Shrink the random loadings of a neuron whose shared variance a^2 + b^2*|B_i|^2
% would exceed max_shared at the top level (keeps the private variance > 0).
room = max_shared - a2max;
assert(room > 0, 'cap_loadings: a^2 = %.3f leaves no room below max_shared = %.2f.', a2max, max_shared);
nb = b2 * sum(B .^ 2, 2);
sc = ones(size(nb));
over = nb > room;
sc(over) = sqrt(room ./ nb(over));
B = B .* sc;
n_cap = nnz(over);
end

function Y = inject_nc(X, labels, a2, b2, B, dr)
% Reorder each neuron's trials within every (rendering, stimulus) block so that
% their ranks follow 10 draws from N(0, R), R = L*L' + D (Iman-Conover).
N = size(B, 1);
L = [sqrt(a2) * ones(N, 1), sqrt(b2) * B];
priv = 1 - sum(L .^ 2, 2);
assert(all(priv > 0), 'inject_nc: shared variance >= 1 for some neuron.');
Y = X;
for c = 1:numel(X)
    for s = 1:max(labels)
        rows = find(labels == s);
        Z = dr.F(:, :, s, c) * L' + dr.E(:, :, s, c) .* sqrt(priv');
        Xb = X{c}(rows, :);
        V  = sort(Xb, 1);
        [~, iz] = sort(Z, 1);
        Yb = zeros(size(Xb));
        Yb(sub2ind(size(Xb), iz, repmat(1:N, numel(rows), 1))) = V;   % k-th smallest count -> trial of the k-th smallest draw
        Y{c}(rows, :) = Yb;
    end
end
end

function [P, is_same] = pair_nc(X, labels, sess)
% Pairwise NC per (pair, rendering x stimulus): Pearson over the repeats of one
% stimulus; NaN when a neuron is silent for that stimulus.
N = size(X{1}, 2);
mask = triu(true(N), 1);
n_stim = max(labels);
P = nan(nnz(mask), n_stim * numel(X));
col = 0;
for c = 1:numel(X)
    for s = 1:n_stim
        col = col + 1;
        C = corr(X{c}(labels == s, :));
        P(:, col) = C(mask);
    end
end
[ii, jj] = find(mask);            % same (column-major) order as C(mask)
is_same = sess(ii) == sess(jj);
end

function st = nc_summary(X, labels, sess, with_pr)
% Achieved NC (Point 7 metric) and, optionally, the noise participation ratio
[P, is_same] = pair_nc(X, labels, sess);
n  = sum(~isnan(P), 2);
m  = mean(P, 2, 'omitnan');
v  = var(P, 0, 2, 'omitnan') ./ n;      % estimation noise of each pair's mean
ok = n > 1;
st.mean_all   = mean(m(ok));
st.sd_raw     = std(m(ok));
st.sd_true    = sqrt(max(var(m(ok)) - mean(v(ok)), 0));
st.mean_same  = mean(m(ok & is_same));
st.mean_cross = mean(m(ok & ~is_same));
n_stim = max(labels);
st.mean_by_cue = zeros(1, numel(X));
for c = 1:numel(X)
    st.mean_by_cue(c) = mean(mean(P(:, (c - 1) * n_stim + (1:n_stim)), 2, 'omitnan'), 'omitnan');
end
st.pair_mean = m;
st.noise_pr = NaN;
if with_pr, st.noise_pr = noise_pr(X, labels); end
end

function pr = noise_pr(X, labels)
% Participation ratio of the residual (noise) correlation matrix, mean over renderings
pr = zeros(1, numel(X));
for c = 1:numel(X)
    R = X{c};
    for s = 1:max(labels)
        rows = labels == s;
        R(rows, :) = R(rows, :) - mean(R(rows, :), 1);
    end
    sd = std(R, 0, 1);
    R  = R(:, sd > 0) ./ sd(sd > 0);
    ev = eig(cov(R));
    pr(c) = sum(ev)^2 / sum(ev .^ 2);
end
pr = mean(pr);
end

function [sd, info] = measure_target_sd(pops, timewindow)
% SD across same-session pairs of each pair's mean NC, minus the estimation
% noise of those means, pooled over the given populations (original data).
m_all = [];  v_all = [];
for i = 1:size(pops, 1)
    pop  = load_population(pops{i, 1}, pops{i, 2}, struct('timewindow', timewindow));
    sess = session_ids(pop);
    [P, is_same] = pair_nc(pop.data, pop.labels, sess);
    n  = sum(~isnan(P), 2);
    ok = is_same & n > 1;
    m_all = [m_all; mean(P(ok, :), 2, 'omitnan')]; %#ok<AGROW>
    v_all = [v_all; var(P(ok, :), 0, 2, 'omitnan') ./ n(ok)]; %#ok<AGROW>
end
sd = sqrt(max(var(m_all) - mean(v_all), 0));
info = struct('source', 'same-session pairs, 4 populations, original data', 'n_pairs', numel(m_all), ...
    'mean', mean(m_all), 'sd_raw', std(m_all), 'est_sd', sqrt(mean(v_all)), 'sd', sd);
end

function calib = calibrate_injection(Xa, labels, sess, cfg, nc_sd, pop_id)
% Rescale a^2 (top level) and b^2 until the achieved mean and spread match the
% targets; the other levels use the same factor on a^2.
k = cfg.k_factors;
top = max(cfg.levels);
[N] = size(Xa{1}, 2);
n_stim = max(labels);  n_rep = sum(labels == 1);
gm = 1;  gs = 1;
hist = nan(cfg.calib_iter, 4);
for it = 1:cfg.calib_iter
    a2 = top * gm;  b2 = nc_sd * gs / sqrt(k);
    mu = zeros(1, cfg.calib_seeds);  sd = mu;
    for cs = 1:cfg.calib_seeds
        dr = make_draws(seed_of(cfg.base_seed, pop_id, 900 + cs, 50), N, k, n_stim, n_rep, numel(Xa));
        B  = cap_loadings(dr.B, a2, b2, cfg.max_shared);
        st = nc_summary(inject_nc(Xa, labels, a2, b2, B, dr), labels, sess, false);
        mu(cs) = st.mean_all;  sd(cs) = st.sd_true;
    end
    hist(it, :) = [mean(mu), mean(sd), gm, gs];
    fprintf('  calibration round %d: achieved mean %.3f (target %.3f), SD %.3f (target %.3f)\n', ...
        it, mean(mu), top, mean(sd), nc_sd);
    if top > 0,   gm = gm * top / mean(mu); end
    if nc_sd > 0, gs = gs * nc_sd / mean(sd); end
end
calib = struct('a2', cfg.levels(:)' * gm, 'b2', nc_sd * gs / sqrt(k), 'gain_mean', gm, 'gain_sd', gs, ...
    'history', array2table(hist, 'VariableNames', {'achieved_mean', 'achieved_sd', 'gain_mean', 'gain_sd'}));
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
% one row per condition x {pooled, each unordered pair}, over the finished seeds
groups = [{'pooled'}; unique(res.pair_group, 'stable')];
done = find(res.done);
rows = {};
for ic = 1:numel(res.cond_names)
    for g = 1:numel(groups)
        if g == 1, P = 1:numel(res.pair_group); else, P = find(strcmp(res.pair_group, groups{g}))'; end
        m = transfer_metrics(res.self(ic, P, done, :), res.pt_rot(ic, P, done, :), ...
            res.ctrl(ic, P, done, :, :), res.nopt(ic, P, done, :));
        rows(end + 1, :) = {res.cond_names{ic}, res.cond_target(ic), groups{g}, numel(done), ...
            mean(res.nc_mean(ic, done)), mean(res.nc_sd_true(ic, done)), mean(res.noise_pr(ic, done)), ...
            m.self, m.nopt, m.pt, m.ctrl, mean(res.ctrl_id(ic, P, done, :, :), 'all', 'omitnan'), ...
            m.ratio, m.diff, m.headroom}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'condition', 'target_nc', 'pair', 'n_seeds', 'nc_mean', 'nc_sd', ...
    'noise_pr', 'self', 'nopt', 'pt_rot', 'ctrl', 'ctrl_id', 'pt_over_ctrl', 'pt_minus_ctrl', 'headroom'});
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
    fprintf('Resuming %s (%d of %d seeds done)\n', file, nnz(res.done), numel(res.done));
end
end

function save_atomic(file, res)
% write to a temporary file first, so an interrupted save never corrupts the results
tmp = [file(1:end - 4) '_tmp.mat'];
save(tmp, 'res', '-v7');
movefile(tmp, file, 'f');
end
