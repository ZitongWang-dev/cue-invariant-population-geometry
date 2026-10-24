%{
Filename: plot_nc_injection.m
Author:   Zitong Wang
Date:     2026-10-24

Description:
    Figures of run_nc_injection.m (rebuttal Point 8: noise correlations
    injected into the recorded pseudo-populations). Conditions: original,
    affine (NC removed, the clean start) and NC injected at mean pairwise
    r = 0, 0.05, 0.10, 0.15. Populations: V1 blue, V2 orange; KO solid, FR
    dashed. Every metric is a ratio of pooled means (transfer_metrics.m);
    error bars are percentile-bootstrap 95% CIs over injection seeds (each seed
    is one injection draw with its own fold partition). Changes from the
    affine condition are paired: the same resampled seeds for both conditions.

Run:
    plot_nc_injection                                   results of the default tag
    nc_tag = 'nc_injection_test'; plot_nc_injection     another tag (e.g. the quick test)

Output (results/figures/09_nc_injection/; 09_<tag>/ for another tag):
    nc_injection_metrics    self-decoding, PT and control, headroom, PT/control per condition
    nc_injection_change     paired change from the affine condition vs achieved mean NC
    nc_injection_check      achieved NC (all / same- / cross-session pairs), its spread,
                            noise participation ratio, per-pair NC distributions
    nc_injection.csv        every plotted number
%}

clearvars -except nc_tag; clc; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'helpers'));
S = figure_settings();
if ~exist('nc_tag', 'var'), nc_tag = 'nc_injection'; end
out_dir = '09_nc_injection';
if ~strcmp(nc_tag, 'nc_injection'), out_dir = ['09_' nc_tag]; end
fprintf('plot_nc_injection (%s)\n', nc_tag);
mstyle = struct('KO', '-', 'FR', '--');

%% Load the finished populations
R = struct();
pops = {};
for ip = 1:numel(S.pops)
    f = fullfile(S.results_root, nc_tag, [S.pops{ip} '.mat']);
    if ~exist(f, 'file'), continue; end
    L = load(f, 'res');
    if ~any(L.res.done), continue; end
    R.(S.pops{ip}) = L.res;
    pops{end + 1} = S.pops{ip}; %#ok<SAGROW>
end
assert(~isempty(pops), 'No finished results in %s', fullfile(S.results_root, nc_tag));
r1 = R.(pops{1});
nC = numel(r1.cond_names);
inj = find(~isnan(r1.cond_target))';
i_orig = find(strcmp(r1.cond_names, 'original'));
i_aff  = find(strcmp(r1.cond_names, 'affine'));
xlab = r1.cond_names;
xlab{i_aff} = 'affine (NC removed)';

%% Bootstrap over seeds
B = struct();
seed = 30000;
for ip = 1:numel(pops)
    seed = seed + 1;
    B.(pops{ip}) = seed_bootstrap(R.(pops{ip}), i_aff, S, seed);
end
offs = ((1:numel(pops)) - (numel(pops) + 1) / 2) * 0.09;

%% Figure 1: metrics per condition
spec = {'self', 'Self-decoding', NaN; 'pt', 'PT (solid) and control (dotted)', NaN; ...
        'headroom', 'Headroom', 0; 'ratio', 'PT / control', 1};
fig = figure('Visible', 'off', 'Position', [100 100 1500 430], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:size(spec, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(pops));
    for ip = 1:numel(pops)
        [c, ls] = pop_style(pops{ip}, S, mstyle);
        b = B.(pops{ip});
        x = (1:nC) + offs(ip);
        h(ip) = draw_conditions(ax, x, b, spec{s, 1}, c, ls, inj, i_orig, i_aff);
        if strcmp(spec{s, 1}, 'pt')
            draw_conditions(ax, x, b, 'ctrl', c, ':', inj, i_orig, i_aff);
        end
    end
    if ~isnan(spec{s, 3}), yline(ax, spec{s, 3}, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off'); end
    set(ax, 'XTick', 1:nC, 'XTickLabel', xlab, 'XTickLabelRotation', 35);
    xlim(ax, [0.5, nC + 0.5]);  box(ax, 'on');  title(ax, spec{s, 2});
    if s == 1, legend(ax, h, strrep(pops, '_', ' '), 'Location', 'best'); end
end
title(tl, sprintf(['Injected noise correlations, whole populations, pooled over 6 cue pairs ' ...
    '(diamond: original, square: NC removed; bars: 95%% CI over %d seeds)'], nnz(r1.done)), 'FontSize', 11);
save_figure(fig, S, out_dir, 'nc_injection_metrics');

%% Figure 2: paired change from the affine (NC removed) condition vs achieved mean NC
spec2 = {'self', '\Delta self-decoding'; 'pt', '\Delta PT'; 'ctrl', '\Delta control'; 'headroom', '\Delta headroom'};
fig = figure('Visible', 'off', 'Position', [100 100 1500 400], 'Color', 'w');
tl  = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for s = 1:size(spec2, 1)
    ax = nexttile(tl);  hold(ax, 'on');
    h = gobjects(1, numel(pops));
    for ip = 1:numel(pops)
        [c, ls] = pop_style(pops{ip}, S, mstyle);
        b = B.(pops{ip});
        x = mean(R.(pops{ip}).nc_mean(:, R.(pops{ip}).done), 2)';
        y = b.d_est.(spec2{s, 1});  lo = b.d_lo.(spec2{s, 1});  hi = b.d_hi.(spec2{s, 1});
        errorbar(ax, x([i_aff, inj]), y([i_aff, inj]), y([i_aff, inj]) - lo([i_aff, inj]), ...
            hi([i_aff, inj]) - y([i_aff, inj]), 'LineStyle', 'none', 'Color', c, 'CapSize', 3, 'HandleVisibility', 'off');
        errorbar(ax, x(i_orig), y(i_orig), y(i_orig) - lo(i_orig), hi(i_orig) - y(i_orig), ...
            'LineStyle', 'none', 'Color', c, 'CapSize', 3, 'HandleVisibility', 'off');
        h(ip) = plot(ax, x(inj), y(inj), ls, 'Color', c, 'Marker', 'o', 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
        plot(ax, x(i_orig), y(i_orig), 'd', 'Color', c, 'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    yline(ax, 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
    box(ax, 'on');  xlabel(ax, 'Achieved mean pairwise NC');  title(ax, spec2{s, 2});
    if s == 1, legend(ax, h, strrep(pops, '_', ' '), 'Location', 'best'); end
end
title(tl, ['Change from the NC-removed (affine) data, paired over seeds; circles: injected levels, ' ...
    'diamond: original data'], 'FontSize', 11);
save_figure(fig, S, out_dir, 'nc_injection_change');

%% Figure 3: manipulation check
fig = figure('Visible', 'off', 'Position', [100 100 1300 820], 'Color', 'w');
tl  = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% (A) achieved mean NC: all pairs, same-session, cross-session
ax = nexttile(tl);  hold(ax, 'on');
h = gobjects(1, numel(pops) + 3);
for ip = 1:numel(pops)
    [c, ls] = pop_style(pops{ip}, S, mstyle);
    r = R.(pops{ip});  d = r.done;  x = (1:nC) + offs(ip);
    h(ip) = plot(ax, x, mean(r.nc_mean(:, d), 2), ls, 'Color', c, 'Marker', 'o', 'MarkerFaceColor', c, 'MarkerSize', 4);
    plot(ax, x, mean(r.nc_mean_same(:, d), 2), '^', 'Color', c, 'MarkerSize', 5, 'HandleVisibility', 'off');
    plot(ax, x, mean(r.nc_mean_cross(:, d), 2), 'v', 'Color', c, 'MarkerSize', 5, 'HandleVisibility', 'off');
end
plot(ax, inj, r1.cond_target(inj), '_', 'Color', 'k', 'MarkerSize', 14, 'LineWidth', 1.5, 'HandleVisibility', 'off');
h(end - 2) = plot(ax, NaN, NaN, 'ko', 'MarkerFaceColor', 'k');
h(end - 1) = plot(ax, NaN, NaN, 'k^');
h(end)     = plot(ax, NaN, NaN, 'kv');
yline(ax, 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
set(ax, 'XTick', 1:nC, 'XTickLabel', xlab, 'XTickLabelRotation', 35);  xlim(ax, [0.5, nC + 0.5]);  box(ax, 'on');
ylabel(ax, 'Mean pairwise NC');  title(ax, 'Achieved NC (black bars: targets)');
legend(ax, h, [strrep(pops, '_', ' '), {'all pairs', 'same-session', 'cross-session'}], 'Location', 'northwest', 'NumColumns', 2);
% (B) spread across pairs (estimation noise removed)
ax = nexttile(tl);  hold(ax, 'on');
for ip = 1:numel(pops)
    [c, ls] = pop_style(pops{ip}, S, mstyle);
    r = R.(pops{ip});  d = r.done;  x = (1:nC) + offs(ip);
    plot(ax, x, mean(r.nc_sd_true(:, d), 2), ls, 'Color', c, 'Marker', 'o', 'MarkerFaceColor', c, 'MarkerSize', 4);
end
yline(ax, r1.target.sd, '--k', sprintf('target %.3f', r1.target.sd), 'LabelHorizontalAlignment', 'left');
set(ax, 'XTick', 1:nC, 'XTickLabel', xlab, 'XTickLabelRotation', 35);  xlim(ax, [0.5, nC + 0.5]);  box(ax, 'on');
ylabel(ax, 'SD of pairwise NC across pairs');  title(ax, 'Spread of NC (estimation noise removed)');
% (C) noise dimensionality
ax = nexttile(tl);  hold(ax, 'on');
for ip = 1:numel(pops)
    [c, ls] = pop_style(pops{ip}, S, mstyle);
    r = R.(pops{ip});  d = r.done;  x = (1:nC) + offs(ip);
    plot(ax, x, mean(r.noise_pr(:, d), 2), ls, 'Color', c, 'Marker', 'o', 'MarkerFaceColor', c, 'MarkerSize', 4);
end
set(ax, 'XTick', 1:nC, 'XTickLabel', xlab, 'XTickLabelRotation', 35);  xlim(ax, [0.5, nC + 0.5]);  box(ax, 'on');
ylabel(ax, 'Participation ratio');  title(ax, 'Dimensionality of the trial-to-trial noise');
% (D) per-pair NC distributions of one population (seed 1)
ax = nexttile(tl);  hold(ax, 'on');
r = r1;  edges = -0.4:0.02:0.6;
sets = {r.pair_nc_seed1(r.is_same_pair, i_orig), 'original, same-session pairs', [0.2 0.2 0.2], '-'; ...
        r.pair_nc_seed1(~r.is_same_pair, i_orig), 'original, cross-session pairs', [0.6 0.6 0.6], '-'; ...
        r.pair_nc_seed1(:, i_aff), 'affine (NC removed)', [0.6 0.6 0.6], ':'};
cl = lines(numel(inj));
for l = 1:numel(inj)
    sets(end + 1, :) = {r.pair_nc_seed1(:, inj(l)), r.cond_names{inj(l)}, cl(l, :), '-'}; %#ok<SAGROW>
end
for i = 1:size(sets, 1)
    histogram(ax, sets{i, 1}, edges, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', ...
        'EdgeColor', sets{i, 3}, 'LineStyle', sets{i, 4}, 'LineWidth', 1.3, 'DisplayName', sets{i, 2});
end
xline(ax, 0, '-', 'Color', [0.3 0.3 0.3], 'HandleVisibility', 'off');
box(ax, 'on');  xlabel(ax, 'Per-pair mean NC (50 stimuli x 3 renderings)');  ylabel(ax, 'Density');
title(ax, sprintf('%s, seed 1: per-pair NC', strrep(pops{1}, '_', ' ')));
legend(ax, 'Location', 'northeast');
title(tl, 'Manipulation check: the injection reaches the target NC and changes the noise structure', 'FontSize', 11);
save_figure(fig, S, out_dir, 'nc_injection_check');

%% Table of every plotted number
rows = {};
for ip = 1:numel(pops)
    r = R.(pops{ip});  b = B.(pops{ip});  d = r.done;
    for ic = 1:nC
        rows(end + 1, :) = {strrep(pops{ip}, '_', ' '), r.cond_names{ic}, r.cond_target(ic), nnz(d), ...
            mean(r.nc_mean(ic, d)), mean(r.nc_sd_true(ic, d)), mean(r.nc_mean_same(ic, d)), ...
            mean(r.nc_mean_cross(ic, d)), mean(r.noise_pr(ic, d)), ...
            b.est.self(ic), b.lo.self(ic), b.hi.self(ic), b.est.pt(ic), b.lo.pt(ic), b.hi.pt(ic), ...
            b.est.ctrl(ic), b.lo.ctrl(ic), b.hi.ctrl(ic), b.est.ratio(ic), b.lo.ratio(ic), b.hi.ratio(ic), ...
            b.est.headroom(ic), b.lo.headroom(ic), b.hi.headroom(ic), ...
            b.d_est.headroom(ic), b.d_lo.headroom(ic), b.d_hi.headroom(ic), ...
            b.d_est.self(ic), b.d_est.pt(ic), b.d_est.ctrl(ic)}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'population', 'condition', 'target_nc', 'n_seeds', 'nc_mean', 'nc_sd', ...
    'nc_same', 'nc_cross', 'noise_pr', 'self', 'self_lo', 'self_hi', 'pt', 'pt_lo', 'pt_hi', ...
    'ctrl', 'ctrl_lo', 'ctrl_hi', 'pt_over_ctrl', 'ratio_lo', 'ratio_hi', 'headroom', 'headroom_lo', 'headroom_hi', ...
    'd_headroom', 'd_headroom_lo', 'd_headroom_hi', 'd_self', 'd_pt', 'd_ctrl'});
save_table(T, S, out_dir, 'nc_injection');

%% ===========================================================================
function [c, ls] = pop_style(pop, S, mstyle)
% V1 blue / V2 orange; KO solid / FR dashed (as the other rev figures)
parts = strsplit(pop, '_');
c  = S.c.(parts{2});
ls = mstyle.(parts{1});
end

function h = draw_conditions(ax, x, b, f, c, ls, inj, i_orig, i_aff)
% injected levels as a connected line; original (diamond) and affine (square) apart
y = b.est.(f);  lo = b.lo.(f);  hi = b.hi.(f);
errorbar(ax, x, y, y - lo, hi - y, 'LineStyle', 'none', 'Color', c, 'CapSize', 3, 'HandleVisibility', 'off');
h = plot(ax, x(inj), y(inj), ls, 'Color', c, 'Marker', 'o', 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
plot(ax, x(i_aff), y(i_aff), 's', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 6, 'HandleVisibility', 'off');
plot(ax, x(i_orig), y(i_orig), 'd', 'Color', c, 'MarkerSize', 7, 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

function out = seed_bootstrap(r, ref, S, seed)
% Point estimates, percentile-bootstrap 95% CIs over seeds, and paired changes
% from condition ref (the same resampled seeds for both conditions).
done = find(r.done);
nSd  = numel(done);
nC   = numel(r.cond_names);
Cm   = mean(r.ctrl, 5, 'omitnan');
f    = {'self', 'pt', 'ctrl', 'headroom', 'ratio'};
E    = metrics_on(r, Cm, done, nC);
rs   = RandStream('twister', 'Seed', seed);
BB   = zeros(S.n_boot, nC, numel(f));
for b = 1:S.n_boot
    BB(b, :, :) = reshape(metrics_on(r, Cm, done(randi(rs, nSd, nSd, 1)), nC), [1 nC numel(f)]);
end
DB = BB - BB(:, ref, :);
for i = 1:numel(f)
    out.est.(f{i})   = E(:, i)';
    out.lo.(f{i})    = prctile(BB(:, :, i), 2.5, 1);
    out.hi.(f{i})    = prctile(BB(:, :, i), 97.5, 1);
    out.d_est.(f{i}) = E(:, i)' - E(ref, i);
    out.d_lo.(f{i})  = prctile(DB(:, :, i), 2.5, 1);
    out.d_hi.(f{i})  = prctile(DB(:, :, i), 97.5, 1);
end
end

function V = metrics_on(r, Cm, sel, nC)
% pooled metrics of every condition over the seeds sel (repeats allowed)
V = zeros(nC, 5);
for ic = 1:nC
    m = transfer_metrics(r.self(ic, :, sel, :), r.pt_rot(ic, :, sel, :), Cm(ic, :, sel, :));
    V(ic, :) = [m.self, m.pt, m.ctrl, m.headroom, m.ratio];
end
end
