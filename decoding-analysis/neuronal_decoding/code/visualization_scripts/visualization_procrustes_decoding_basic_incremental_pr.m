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

script_dir = fileparts(mfilename('fullpath'));

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
saveas(gcf, fullfile(fig_dir, 'fig1_accuracy_pr_vs_N.png'));

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
saveas(gcf, fullfile(fig_dir, 'fig2_accuracy_vs_pr.png'));

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