%{
Filename: visualization_procrustes_decoding_basic_incremental_pr.m
Author:   Zitong Wang
Date:     2025-07

Description:
    Figures for the incremental neural-vs-Gabor participation-ratio (PR)
    comparison. Figure 1 (below) plots accuracy and PR against population
    size on a dual y-axis; Figure 2 (accuracy vs PR) will be appended here.

    FIGURE 1 -- accuracy and PR vs population size (dual y-axis):
      - LEFT  axis : decoding accuracy
                     * rotation-only Procrustes transfer (solid, marker 'o')
                     * within-condition self-decoding     (dashed, marker 's')
      - RIGHT axis : participation ratio (PR)              (dotted, marker '^')
    All populations are overlaid, colour-coded:
      FR V1, FR V2, KO V1, KO V2, and the Gabor model (black).

    Accuracy is the rotation-only transfer (column 6 of .acc) to match the
    main results; self-decoding is column 1 (genAcc). Each quantity is
    averaged across the six ordered rendering pairs (PR is pooled across the
    two conditions and the six pairs). Error bars are +/-1 SEM over the pooled
    resampling distribution (std / sqrt(n)) -- the same convention as
    visualization_procrustes_decoding_basic_incremental.

Path note:
    This script is expected to live in
        neuronal_decoding/code/visualization_scripts/
    Neural results tree:  neuronal_decoding/results/decoding_outputs/...   (../../results)
    Gabor  results tree:  Gabor_decoding/results/decoding_outputs/...      (../../../Gabor_decoding/results)
    Paths below are anchored to this file's own location (mfilename), so the
    current working directory does not matter.
%}

%% Configuration
clc; clear;

% Resolve this script's folder. mfilename('fullpath') only points at the saved
% file when the script is Run as a whole file; running cell-by-cell or pasting
% into the Command Window returns a temp editor path, so fall back to pwd
% (assumes the working directory is .../code/visualization_scripts/).
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir) || contains(script_dir, fullfile('AppData','Local','Temp'))
    script_dir = pwd;
end

% --- Result-tree base paths (anchored to this script's location) ---
neural_base = fullfile(script_dir, '..','..','results','decoding_outputs', ...
    'procrustes_decoding_basic_incremental_pr_results');
gabor_base  = fullfile(script_dir, '..','..','..','Gabor_decoding','results','decoding_outputs', ...
    'procrustes_Gabor_incremental_pr_results');

% --- Gabor model to overlay ---
gabor_filter_name = 'even and odd combined';

% --- Populations: {display name, directory holding its 6 pair files, colour} ---
populations = {
    'FR V1',            fullfile(neural_base,'FR','V1'),            [0.00 0.45 0.74];
    'FR V2',            fullfile(neural_base,'FR','V2'),            [0.30 0.75 0.93];
    'KO V1',            fullfile(neural_base,'KO','V1'),            [0.72 0.10 0.10];
    'KO V2',            fullfile(neural_base,'KO','V2'),            [0.95 0.50 0.50];
    'Gabor (even+odd)', fullfile(gabor_base, gabor_filter_name),    [0.00 0.00 0.00];
};

% --- Metric columns inside .acc ---
col_self = 1;   % within-condition self-decoding (genAcc)
col_rot  = 6;   % rotation-only Procrustes transfer

% --- The six ordered-pair files (variable name inside each == file stem) ---
pair_files = {'acec_results','ecex_results','acex_results', ...
              'ecac_results','exec_results','exac_results'};

% --- Display options ---
show_err = true;          % overlay +/-1 SEM error bars over the pooled resampling distribution
show_endpoint_N = true;   % (Figure 2) label each curve's largest-N endpoint with its N
n_boot = 2000;            % (Figures 3-4) bootstrap iterations for ratio CIs
fit_method = 'mm';    % (Figure 5) curve-fit form: 'satexp' = a*(1-exp(-N/tau)), 'mm' = a*N/(N+k)

% --- Output figure location ---
fig_dir = fullfile(script_dir, '..','..','results','figures','incremental_pr');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Load and pool each population
nPop = size(populations,1);
P = struct('name',{},'color',{},'N',{}, ...
           'self_m',{},'self_sem',{},'rot_m',{},'rot_sem',{},'pr_m',{},'pr_sem',{});

for p = 1:nPop
    pdir = populations{p,2};
    if ~exist(pdir,'dir')
        warning('Missing directory, skipping %s: %s', populations{p,1}, pdir);
        continue;
    end
    try
        [N, self_m, self_sem, rot_m, rot_sem, pr_m, pr_sem] = ...
            pool_population(pdir, pair_files, col_self, col_rot);
    catch ME
        warning('Could not load %s (%s); skipping.', populations{p,1}, ME.message);
        continue;
    end
    P(end+1) = struct('name',populations{p,1},'color',populations{p,3}, ...
        'N',N,'self_m',self_m,'self_sem',self_sem, ...
        'rot_m',rot_m,'rot_sem',rot_sem,'pr_m',pr_m,'pr_sem',pr_sem);
    fprintf('Loaded %-18s : %d sizes (N = %d ... %d)\n', ...
        populations{p,1}, numel(N), N(1), N(end));
end

if isempty(P)
    error('No populations loaded. Check neural_base / gabor_base paths.');
end

%% Plot: dual-y accuracy (left) + PR (right) vs number of units
figure('Color','w','Position',[100 100 950 620]);
ax = gca; hold(ax,'on');

pop_handles = gobjects(1,numel(P));   % colour proxies for the population legend

for p = 1:numel(P)
    c = P(p).color;
    x = P(p).N;

    % ----- LEFT axis: accuracy (self solid, PT dashed) -----
    yyaxis left
    if show_err
        errorbar(x, P(p).self_m, P(p).self_sem, '-',  'Color',c, 'Marker','s', ...
            'MarkerFaceColor',c, 'MarkerSize',4, 'LineWidth',1.5, 'CapSize',3);
        errorbar(x, P(p).rot_m,  P(p).rot_sem,  '--', 'Color',c, 'Marker','o', ...
            'MarkerSize',4, 'LineWidth',1.5, 'CapSize',3);
    else
        plot(x, P(p).self_m, '-',  'Color',c, 'Marker','s', ...
            'MarkerFaceColor',c, 'MarkerSize',4, 'LineWidth',1.5);
        plot(x, P(p).rot_m,  '--', 'Color',c, 'Marker','o', ...
            'MarkerSize',4, 'LineWidth',1.5);
    end
    pop_handles(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);  % colour proxy (population)

    % ----- RIGHT axis: PR -----
    yyaxis right
    if show_err
        errorbar(x, P(p).pr_m, P(p).pr_sem, ':', 'Color',c, 'Marker','^', ...
            'MarkerSize',4, 'LineWidth',1.3, 'CapSize',3);
    else
        plot(x, P(p).pr_m, ':', 'Color',c, 'Marker','^', ...
            'MarkerSize',4, 'LineWidth',1.3);
    end
end

% ----- Axis styling -----
yyaxis left
ylabel('Decoding accuracy');
ylim([0 1]);
ax.YAxis(1).Color = 'k';

yyaxis right
ylabel('Participation ratio');
yl = ylim; ylim([0 yl(2)]);
ax.YAxis(2).Color = 'k';

xlabel('Number of units (neurons / filters)');
grid on; box on;

% ----- Linestyle key: fake black lines naming the metric for each style -----
yyaxis left
h_self = plot(nan, nan, '-k',  'Marker','s', 'MarkerFaceColor','k', 'LineWidth',1.5);
h_pt   = plot(nan, nan, '--k', 'Marker','o', 'LineWidth',1.5);
h_pr   = plot(nan, nan, ':k',  'Marker','^', 'LineWidth',1.5);

legend([pop_handles, h_self, h_pt, h_pr], ...
       [{P.name}, {'self-decoding','PT-decoding','PR'}], ...
       'Location','best');
sgtitle('Decoding accuracy and participation ratio vs population size');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig1_accuracy_pr_vs_N.png'));

%% Figure 2: transfer (rotation) and self-decoding vs participation ratio
%  Reuses the pooled P struct from Figure 1. Each population is a connected
%  curve through the (PR, accuracy) plane, ordered by N (smallest N at the
%  lower-left). 2-D error bars are +/-1 SEM on both axes.
figure('Color','w','Position',[100 100 1180 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% Panel A: PR vs rotation-only PT transfer
axA = nexttile;
leg_handles = plot_acc_vs_pr(P, 'rot_m', 'rot_sem', 'o', show_err, show_endpoint_N);
xlabel('Participation ratio'); ylabel('Rotation-only PT transfer accuracy');
title('Cue-transfer (rotation) vs effective dimensionality');
grid on; box on;

% Panel B: PR vs self-decoding
nexttile;
plot_acc_vs_pr(P, 'self_m', 'self_sem', 's', show_err, show_endpoint_N);
xlabel('Participation ratio'); ylabel('Self-decoding accuracy');
title('Self-decoding vs effective dimensionality');
grid on; box on;

lg = legend(axA, leg_handles, {P.name});
lg.Layout.Tile = 'east';
sgtitle('Decoding accuracy vs participation ratio');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig2_accuracy_vs_pr.png'));

%% Load per-unit (neuron-repeat) data for the ratio figures
%  Self-contained reload of the six pair files per population so Figures 1-2
%  stay untouched. A "unit" is one neuron-resampling repeat within one pair
%  file; this is the resampling unit used for the bootstrap CIs below.
U = struct('name',{},'color',{},'N',{},'pr_m',{},'pr_sem',{}, ...
           'unit_pt',{},'unit_self',{},'unit_pr',{});
for p = 1:nPop
    pdir = populations{p,2};
    if ~exist(pdir,'dir'), continue; end
    try
        Up = load_units(pdir, pair_files, col_self, col_rot);
    catch ME
        warning('Could not load units for %s (%s); skipping.', populations{p,1}, ME.message);
        continue;
    end
    U(end+1) = struct('name',populations{p,1},'color',populations{p,3}, ...
        'N',Up.N,'pr_m',Up.pr_m,'pr_sem',Up.pr_sem, ...
        'unit_pt',Up.unit_pt,'unit_self',Up.unit_self,'unit_pr',Up.unit_pr);
end

rng(1);   % reproducible bootstrap CIs

%% Figure 3: normalized cue-transfer (PT / self) vs participation ratio
%  Per N: pooled mean PT (rotation transfer over all six directions) divided
%  by pooled mean self-decoding (all three conditions) -- the symmetric
%  PT/mean(self) ratio-of-means. y error = percentile bootstrap over
%  neuron-repeats; x = pooled PR mean (x error = PR SEM, Figure-2 convention).
figure('Color','w','Position',[100 100 720 560]); hold on;
leg_h = gobjects(1,numel(U));
for p = 1:numel(U)
    c = U(p).color;
    [r, lo, hi] = bootstrap_ratio(U(p).unit_pt, U(p).unit_self, n_boot);
    errorbar(U(p).pr_m, r, r-lo, hi-r, U(p).pr_sem, U(p).pr_sem, '-', ...
        'Color',c, 'Marker','o', 'MarkerFaceColor',c, 'MarkerSize',4, ...
        'LineWidth',1.4, 'CapSize',3);
    leg_h(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);
end
xlabel('Participation ratio');
ylabel('Normalized transfer   (PT / self-decoding)');
title('Cue-transfer efficiency vs effective dimensionality');
grid on; box on;
legend(leg_h, {U.name}, 'Location','best');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig3_normalized_transfer_vs_pr.png'));

%% Figure 4: transfer-per-dimension (PT / PR) vs neuron count
%  Matched-population-size view requested by the reviewer: at each N, pooled
%  mean PT divided by pooled mean PR (a transfer-per-effective-dimension
%  index). y error = percentile bootstrap over neuron-repeats.
figure('Color','w','Position',[100 100 720 560]); hold on;
leg_h = gobjects(1,numel(U));
for p = 1:numel(U)
    c = U(p).color;
    [r, lo, hi] = bootstrap_ratio(U(p).unit_pt, U(p).unit_pr, n_boot);
    errorbar(U(p).N, r, r-lo, hi-r, '-', ...
        'Color',c, 'Marker','o', 'MarkerFaceColor',c, 'MarkerSize',4, ...
        'LineWidth',1.4, 'CapSize',3);
    leg_h(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);
end
xlabel('Number of units (neurons / filters)');
ylabel('Transfer per effective dimension   (PT / PR)');
title('Transfer per dimension vs population size');
grid on; box on;
legend(leg_h, {U.name}, 'Location','best');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig4_transfer_per_dim_vs_N.png'));

%% Figure 5: where PR and PT saturate, viewed four ways
%  (a) PR vs N  -> horizontal asymptote PR_inf (PR(N) fit)
%  (b) PT vs N  -> horizontal asymptote PT_inf (PT(N) fit; usually unconverged)
%  (c) PT vs PR -> vertical line = PR_inf carried over from the PR(N) fit,
%      plus the limiting corner (PR_inf, PT_inf).
%  (d) PT vs PR -> vertical line read purely from the trajectory's shape:
%      extend the curve and stop where its tail turns parallel to the y-axis
%      (slope dPT/dPR exceeds vert_slope), regardless of whether PT is a
%      reachable accuracy. (c) and (d) are two readings of the same limiting
%      dimensionality -- (c) from the independent PR(N) fit, (d) from the
%      joint curve -- so comparing them is a self-consistency check.
nP = numel(P);
prinf_all = zeros(1,nP); ptinf_all = zeros(1,nP);
prrate_all = zeros(1,nP); ptrate_all = zeros(1,nP);
prfun_all = cell(1,nP);  ptfun_all = cell(1,nP);
fprintf('\nCurve fits (method = %s):\n', fit_method);
for p = 1:nP
    [prinf_all(p), prrate_all(p), prfun_all{p}] = fit_curve(P(p).N, P(p).pr_m,  fit_method);
    [ptinf_all(p), ptrate_all(p), ptfun_all{p}] = fit_curve(P(p).N, P(p).rot_m, fit_method);
    fprintf('  %-18s PR_inf = %5.2f   PT_inf = %5.2f\n', P(p).name, prinf_all(p), ptinf_all(p));
end

ycap = 1.05;        % display cap on the accuracy axes

figure('Color','w','Position',[100 100 1300 760]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax_a = nexttile(1);   % (a) PR vs N            [top-left]
ax_c = nexttile(2);   % (c) PT vs PR, PR_inf   [top-right]
ax_b = nexttile(3);   % (b) PT vs N            [bottom-left]
ax_d = nexttile(4);   % (d) PT vs PR, tail     [bottom-right]

% ---- (a) PR vs N : horizontal asymptote PR_inf ----
axes(ax_a); hold on;
for p = 1:nP
    c = P(p).color; Nv = P(p).N;
    Ng = linspace(min(Nv), max(Nv)*1.5, 300);
    plot(Nv, P(p).pr_m, 'o', 'Color',c, 'MarkerFaceColor',c, 'MarkerSize',4);
    plot(Ng, prfun_all{p}(Ng), '-', 'Color',c, 'LineWidth',1.4);
    yline(prinf_all(p), '--', sprintf('%.1f',prinf_all(p)), 'Color',c, 'FontSize',8, ...
        'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
end
ylabel('Participation ratio');
title('(a) PR vs N   (asymptote PR_\infty)'); grid on; box on;

% ---- (b) PT vs N : horizontal asymptote PT_inf ----
axes(ax_b); hold on;
for p = 1:nP
    c = P(p).color; Nv = P(p).N;
    Ng = linspace(min(Nv), max(Nv)*1.5, 300);
    plot(Nv, P(p).rot_m, 'o', 'Color',c, 'MarkerFaceColor',c, 'MarkerSize',4);
    plot(Ng, ptfun_all{p}(Ng), '-', 'Color',c, 'LineWidth',1.4);
    if ptinf_all(p) <= ycap
        yline(ptinf_all(p), '--', sprintf('%.2f',ptinf_all(p)), 'Color',c, 'FontSize',8, ...
            'LabelHorizontalAlignment','left', 'LabelVerticalAlignment','bottom');
    end
end
xlabel('Number of units (N)'); ylabel('Rotation-only PT transfer accuracy');
title('(b) PT vs N   (asymptote PT_\infty, often unconverged)'); grid on; box on;
ylim([0 ycap]);
linkaxes([ax_a ax_b],'x'); ax_a.XTickLabel = [];   % shared N axis, label only on (b)

% ---- (c) PT vs PR : vertical line = PR_inf (from the PR(N) fit) + corner ----
axes(ax_c); hold on;
leg_h = gobjects(1,nP);
for p = 1:nP
    c = P(p).color; Nv = P(p).N;
    prinf = prinf_all(p); ptinf = ptinf_all(p);
    prfun = prfun_all{p}; ptfun = ptfun_all{p};

    plot(P(p).pr_m, P(p).rot_m, 'o', 'Color',c, 'MarkerFaceColor',c, 'MarkerSize',4);

    % extend N until PR saturates (and a little past PT=1 if it crosses)
    N_pt1  = invert_fit(ptinf, ptrate_all(p), fit_method, 1.0);
    N_pr99 = invert_fit(prinf, prrate_all(p), fit_method, 0.999*prinf);
    N_max = N_pr99;
    if isfinite(N_pt1), N_max = max(N_max, N_pt1*1.10); end
    N_max = max(min(N_max, max(Nv)*1000), max(Nv)*1.5);
    Ng = linspace(min(Nv), N_max, 600);
    PRg = prfun(Ng); PTg = ptfun(Ng);
    plot_traj(PRg, PTg, c);

    xline(prinf, '--', sprintf('%.1f',prinf), 'Color',c, 'FontSize',8, ...
        'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','top');
    if ptinf <= ycap
        plot(prinf, ptinf, 'p', 'Color',c, 'MarkerFaceColor',c, 'MarkerSize',10);
    end
    leg_h(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);
end
xlabel('Participation ratio'); ylabel('Rotation-only PT transfer accuracy');
title('(c) PT vs PR   (vertical = PR_\infty from PR-vs-N fit)'); grid on; box on;
ylim([0 ycap]);

% ---- (d) PT vs PR : vertical asymptote from extrapolating PR(PT) to PT -> inf ----
%  Fit PR as a function of PT (axes swapped) and read its asymptote PR_v -- the
%  PR the joint curve approaches as PT grows without bound. This is the vertical
%  asymptote implied by the PR-vs-PT shape itself, independent of the PR(N) fit
%  used in (c). It is a long extrapolation from the observed PT range out to
%  PT -> inf, so PR_v is a soft, curve-intrinsic estimate.
axes(ax_d); hold on;
for p = 1:nP
    c = P(p).color;
    [prv, ~, prvfun] = fit_curve(P(p).rot_m, P(p).pr_m, fit_method);   % x = PT, y = PR

    plot(P(p).pr_m, P(p).rot_m, 'o', 'Color',c, 'MarkerFaceColor',c, 'MarkerSize',4);

    % sweep PT past 1 (ignoring reachability) so the curve bends toward PR_v
    PTg = linspace(min(P(p).rot_m), 1.5, 400);
    PRg = prvfun(PTg);
    plot_traj(PRg, PTg, c);   % solid PT<=1, dashed in the PT>1 extrapolation

    xline(prv, '--', sprintf('%.1f',prv), 'Color',c, 'FontSize',8, ...
        'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','top');
end
xlabel('Participation ratio'); ylabel('Rotation-only PT transfer accuracy');
title('(d) PT vs PR   (vertical = PR_v, extrapolating PR(PT) to PT\rightarrow\infty)'); grid on; box on;
ylim([0 ycap]);
linkaxes([ax_c ax_d],'y');   % shared accuracy axis; x autoscales per panel

lg = legend(ax_c, leg_h, {P.name}); lg.Layout.Tile = 'east';
sgtitle('Saturation of PR and PT, and two readings of the PR-vs-PT asymptote');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig5_saturation_four_views.png'));

%% ----------------------------------------------------------------------
function [N, self_m, self_sem, rot_m, rot_sem, pr_m, pr_sem] = ...
    pool_population(pdir, pair_files, col_self, col_rot)
% Load the six ordered-pair result files in pdir, pool across pairs, and
% return per-neuron-count mean and SEM (std/sqrt(n) over the pooled
% resampling distribution) for self-decoding accuracy, rotation-transfer
% accuracy, and PR. Accuracy pools all resampling rows across the six pairs;
% PR pools both conditions across the six pairs -- the same SEM convention
% as visualization_procrustes_decoding_basic_incremental.

nPairs = numel(pair_files);
res = cell(1,nPairs);
for k = 1:nPairs
    f = fullfile(pdir, [pair_files{k} '.mat']);
    S = load(f);
    res{k} = S.(pair_files{k});   % 1 x nSeq cell of structs
end

nSeq     = numel(res{1});
N        = zeros(1,nSeq);
self_m   = zeros(1,nSeq); self_sem = zeros(1,nSeq);
rot_m    = zeros(1,nSeq); rot_sem  = zeros(1,nSeq);
pr_m     = zeros(1,nSeq); pr_sem   = zeros(1,nSeq);

for s = 1:nSeq
    N(s) = res{1}{s}.neuron_num;

    self_dist = [];   % all resampling rows, pooled across pairs
    rot_dist  = [];
    pr_dist   = [];   % per-repeat PR, pooled across the two conditions and six pairs
    for k = 1:nPairs
        c = res{k}{s};
        self_dist = [self_dist; c.acc(:,col_self)];
        rot_dist  = [rot_dist;  c.acc(:,col_rot )];
        pr_dist   = [pr_dist;   c.pr_stim1(:); c.pr_stim2(:)];
    end

    self_m(s) = mean(self_dist); self_sem(s) = std(self_dist)/sqrt(numel(self_dist));
    rot_m(s)  = mean(rot_dist);  rot_sem(s)  = std(rot_dist)/sqrt(numel(rot_dist));
    pr_m(s)   = mean(pr_dist);   pr_sem(s)   = std(pr_dist)/sqrt(numel(pr_dist));
end

% Ensure ascending population size for clean line plotting
[N, order] = sort(N);
self_m = self_m(order); self_sem = self_sem(order);
rot_m  = rot_m(order);  rot_sem  = rot_sem(order);
pr_m   = pr_m(order);   pr_sem   = pr_sem(order);
end

function h = plot_acc_vs_pr(P, accfield, semfield, mk, show_err, label_end)
% Plot accuracy-vs-PR connected curves (one per population) with 2-D SEM
% error bars. accfield/semfield name the accuracy mean/SEM fields in P (e.g.
% 'rot_m'/'rot_sem' or 'self_m'/'self_sem'); mk is the marker. Points are
% ordered by N (smallest N at the lower-left). Returns one clean line handle
% per population for the legend.
hold on;
h = gobjects(1,numel(P));
for p = 1:numel(P)
    c  = P(p).color;
    x  = P(p).pr_m;        xe = P(p).pr_sem;
    y  = P(p).(accfield);  ye = P(p).(semfield);
    if show_err
        errorbar(x, y, ye, ye, xe, xe, '-', 'Color',c, 'Marker',mk, ...
            'MarkerFaceColor',c, 'MarkerSize',4, 'LineWidth',1.4, 'CapSize',3);
    else
        plot(x, y, '-', 'Color',c, 'Marker',mk, 'MarkerFaceColor',c, ...
            'MarkerSize',4, 'LineWidth',1.4);
    end
    h(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);  % clean legend proxy
    if label_end
        text(x(end), y(end), sprintf('  %d', P(p).N(end)), ...
            'Color',c, 'FontSize',8, 'VerticalAlignment','middle');
    end
end
end

function U = load_units(pdir, pair_files, col_self, col_rot)
% Reload the six ordered-pair files and return per-neuron-repeat ("unit")
% means needed for the ratio bootstraps, plus the pooled PR mean/SEM. A unit
% is one neuron-resampling repeat within one pair file (the resampling unit).
% Returns matrices [nUnits x nSeq] for PT, self, and PR per unit, sorted by N.
nPairs = numel(pair_files);
res = cell(1,nPairs);
for k = 1:nPairs
    S = load(fullfile(pdir, [pair_files{k} '.mat']));
    res{k} = S.(pair_files{k});
end

nSeq   = numel(res{1});
nrep   = numel(res{1}{1}.pr_stim1);   % neuron_sample_repeat
nUnits = nPairs * nrep;

N         = zeros(1,nSeq);
pr_m      = zeros(1,nSeq);
pr_sem    = zeros(1,nSeq);
unit_pt   = zeros(nUnits,nSeq);
unit_self = zeros(nUnits,nSeq);
unit_pr   = zeros(nUnits,nSeq);

for s = 1:nSeq
    N(s) = res{1}{s}.neuron_num;
    pr_pool = [];
    u = 0;
    for k = 1:nPairs
        c   = res{k}{s};
        ids = c.acc_repeat_id;
        upt   = accumarray(ids, c.acc(:,col_rot ), [], @mean);  % nrep x 1, per-unit mean PT
        uself = accumarray(ids, c.acc(:,col_self), [], @mean);  % nrep x 1, per-unit mean self
        upr   = (c.pr_stim1(:) + c.pr_stim2(:)) / 2;            % nrep x 1, per-unit mean PR
        rows = u + (1:nrep);
        unit_pt(rows,s)   = upt;
        unit_self(rows,s) = uself;
        unit_pr(rows,s)   = upr;
        u = u + nrep;
        pr_pool = [pr_pool; c.pr_stim1(:); c.pr_stim2(:)];
    end
    pr_m(s)   = mean(pr_pool);
    pr_sem(s) = std(pr_pool) / sqrt(numel(pr_pool));   % Figure-2 convention
end

% Sort by N ascending (match Figures 1-2)
[N, order] = sort(N);
pr_m      = pr_m(order);      pr_sem    = pr_sem(order);
unit_pt   = unit_pt(:,order); unit_self = unit_self(:,order);
unit_pr   = unit_pr(:,order);

U = struct('N',N, 'pr_m',pr_m, 'pr_sem',pr_sem, ...
           'unit_pt',unit_pt, 'unit_self',unit_self, 'unit_pr',unit_pr);
end

function [r, lo, hi] = bootstrap_ratio(num_units, den_units, nboot)
% Ratio-of-means with a percentile bootstrap over units (rows). num_units and
% den_units are [nUnits x nSeq] per-unit means, paired by row. For each column
% the point estimate is mean(num)/mean(den); the CI resamples the units (the
% same indices for numerator and denominator, preserving pairing).
[nUnits, nSeq] = size(num_units);
r  = zeros(1,nSeq);
lo = zeros(1,nSeq);
hi = zeros(1,nSeq);
for s = 1:nSeq
    nu = num_units(:,s);
    de = den_units(:,s);
    r(s) = mean(nu) / mean(de);
    bs = zeros(nboot,1);
    for b = 1:nboot
        idx = randi(nUnits, nUnits, 1);   % resample neuron-repeats with replacement
        bs(b) = mean(nu(idx)) / mean(de(idx));
    end
    lo(s) = prctile(bs, 2.5);
    hi(s) = prctile(bs, 97.5);
end
end

function [yinf, rate, fun] = fit_curve(x, y, method)
% Least-squares fit of a 2-parameter saturating curve via fminsearch (no
% toolbox dependency). Both forms share the asymptote yinf = y(N->inf):
%   'satexp' : y = yinf * (1 - exp(-x/rate))      (rate = tau, time-constant)
%   'mm'     : y = yinf * x / (x + rate)           (rate = k, half-saturation N)
% Returns the asymptote yinf, the rate parameter, and a handle to the fit.
x = x(:); y = y(:);
switch lower(method)
    case 'satexp'
        model = @(b,xx) b(1) * (1 - exp(-xx ./ abs(b(2))));
    case 'mm'
        model = @(b,xx) b(1) * xx ./ (xx + abs(b(2)));
    otherwise
        error('fit_curve:method', 'Unknown fit_method "%s" (use ''satexp'' or ''mm'').', method);
end
sse  = @(b) sum((y - model(b,x)).^2);
b0   = [max(y)*1.05, median(x)];
opts = optimset('Display','off','MaxFunEvals',5000,'MaxIter',5000);
b    = fminsearch(sse, b0, opts);
yinf = b(1);
rate = abs(b(2));
fun  = @(xx) model([yinf, rate], xx);
end

function N = invert_fit(yinf, rate, method, ytarget)
% N at which the fitted curve reaches ytarget; Inf if ytarget >= asymptote.
if ytarget >= yinf
    N = Inf; return;
end
switch lower(method)
    case 'satexp'
        N = -rate * log(1 - ytarget/yinf);
    case 'mm'
        N = rate * ytarget / (yinf - ytarget);
    otherwise
        N = Inf;
end
end

function plot_traj(PRg, PTg, c)
% Plot a PR-vs-PT trajectory: solid where PT is a valid accuracy (<=1),
% dashed in the extrapolated PT>1 region.
k = find(PTg > 1, 1);
if isempty(k)
    plot(PRg, PTg, '-', 'Color',c, 'LineWidth',1.6);
elseif k == 1
    plot(PRg, PTg, '--', 'Color',c, 'LineWidth',1.6);
else
    plot(PRg(1:k-1),   PTg(1:k-1),   '-',  'Color',c, 'LineWidth',1.6);
    plot(PRg(k-1:end), PTg(k-1:end), '--', 'Color',c, 'LineWidth',1.6);
end
end