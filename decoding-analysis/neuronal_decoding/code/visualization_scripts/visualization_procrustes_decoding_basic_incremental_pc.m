%{
Filename: visualization_procrustes_decoding_basic_incremental_pc.m
Author:   Zitong Wang
Date:     2026-07

Description:
    Figures for the PC-dimension Procrustes transfer analysis
    (procrustes_decoding_basic_incremental_pc.m). The x-axis is the number of
    retained principal components, k. All quantities are pooled across the six
    ordered rendering pairs; all four populations overlaid, colour = population.

    FIGURE 1 -- decoding accuracy vs k:
      PANEL 1 (raw accuracy):
        - self-decoding        within-condition genAcc      (.acc col 1, solid,  key 's')
        - PT decoding          rotation-only transfer        (.acc col 6, dashed, key 'o')
        - PT control decoding  rotation-only shuffle null     (.acc col 4, dotted, key '^')
      PANEL 2 (normalized by self-decoding):
        - PT decoding / self-decoding        (dashed)
        - control decoding / self-decoding   (dotted), y = 1 = self-decoding ceiling.

    FIGURE 2 -- residual Procrustes shape distance vs k (.pdist):
      PANEL 1 (raw): d_before (col 1, solid, key 's') and d_rot (col 2, dashed,
                     key 'o'); the two share MATLAB's stim2-target denominator.
      PANEL 2 (after / before): d_rot ./ d_before, y = 1 = no improvement from
                     rotation (lower = rotation shrank the residual mismatch more).
      (col 3, d_full, is not plotted.)

    FIGURE 3 -- ratio summary vs k (two subplots, ratios as percentages):
      PANEL 1 (ratios to self): PT/self (dashed, key 'o') and control/self
                     (dotted, key '^'); both sit at or below 100%.
      PANEL 2 (PT / control, solid, key 's'): how far PT pulls above its own
                     shuffle null; runs above 100%. Split from panel 1 because it
                     lives on a different scale. Both panels share a 100% line.
      PT/self and control/self repeat the Figure 1 panel-2 ratios.

    PT column choice: the shuffle control is rotation-only (T alone, no scale, no
    translation), so its matched real quantity is the rotation-only transfer
    (col 6), NOT the full b*T+c transfer (col 3). col_pt is set to 6 for a
    like-for-like comparison; set col_pt = 3 (and pt_name) to plot the full
    transform instead.

    Conventions match visualization_procrustes_decoding_basic_incremental_pr.m:
      - colour = population, linestyle = quantity, with a black-proxy linestyle key;
      - error is +/-1 SEM over the pooled resampling (trial-split) distribution
        (std/sqrt(n)); normalized/ratio panels use ratio-of-means with a percentile
        bootstrap over the pooled resampling rows (same estimator as the pr script).

    Deliberate departures (documented):
      - x-axis is PC count (.pc_num), not neuron count; no PR fields are read
        (dimension is controlled, not measured);
      - the dense k-axis (up to 49 points) is drawn as mean lines with shaded
        bands rather than per-point error bars, to stay legible; markers appear
        only in the linestyle key, not on the data curves.

Path note:
    Expected to live in neuronal_decoding/code/visualization_scripts/.
    Neural results tree: ../../results/decoding_outputs/...  (anchored via mfilename).
%}

%% Configuration
clc; clear;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir) || contains(script_dir, fullfile('AppData','Local','Temp'))
    script_dir = pwd;
end

neural_base = fullfile(script_dir, '..','..','results','decoding_outputs', ...
    'procrustes_decoding_basic_incremental_pc_results');

% --- Populations: {display name, directory holding its 6 pair files, colour} ---
populations = {
    'FR V1', fullfile(neural_base,'FR','V1'), [0.00 0.45 0.74];
    'FR V2', fullfile(neural_base,'FR','V2'), [0.30 0.75 0.93];
    'KO V1', fullfile(neural_base,'KO','V1'), [0.72 0.10 0.10];
    'KO V2', fullfile(neural_base,'KO','V2'), [0.95 0.50 0.50];
};

% --- Metric columns inside .acc ---
col_self = 1;   % within-condition self-decoding (genAcc)
col_pt   = 6;   % PT decoding: rotation-only transfer (set 3 for full b*T+c)
col_ctrl = 4;   % PT control: rotation-only broken-correspondence shuffle
pt_name  = 'PT (rotation)';   % legend label for the PT curve

% --- Columns inside .pdist (Figure 2) ---
pcol_before = 1;   % d_before (no rotation)
pcol_after  = 2;   % d_rot    (rotation only)

% --- The six ordered-pair files (variable name inside each == file stem) ---
pair_files = {'acec_results','ecex_results','acex_results', ...
              'ecac_results','exec_results','exac_results'};

% --- Display options ---
show_band = true;    % shaded +/-1 SEM (panel 1) / bootstrap CI (panel 2) bands
n_boot    = 2000;    % bootstrap iterations for the normalized-panel CIs
chance    = 1/50;    % 50-class chance level (reference line, panel 1)

fig_dir = fullfile(script_dir, '..','..','results','figures','incremental_pc');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Load and pool each population
nPop = size(populations,1);
P = struct('name',{},'color',{},'k',{}, ...
           'self_m',{},'self_sem',{},'pt_m',{},'pt_sem',{},'ctrl_m',{},'ctrl_sem',{}, ...
           'rpt',{},'rpt_lo',{},'rpt_hi',{},'rct',{},'rct_lo',{},'rct_hi',{}, ...
           'rpc',{},'rpc_lo',{},'rpc_hi',{}, ...
           'bef_m',{},'bef_sem',{},'aft_m',{},'aft_sem',{}, ...
           'rd',{},'rd_lo',{},'rd_hi',{});

for p = 1:nPop
    pdir = populations{p,2};
    if ~exist(pdir,'dir')
        warning('Missing directory, skipping %s: %s', populations{p,1}, pdir);
        continue;
    end
    try
        [kx, self_mat, pt_mat, ctrl_mat, bef_mat, aft_mat] = ...
            pool_population(pdir, pair_files, col_self, col_pt, col_ctrl, pcol_before, pcol_after);
    catch ME
        warning('Could not load %s (%s); skipping.', populations{p,1}, ME.message);
        continue;
    end

    self_m = mean(self_mat,1); self_sem = std(self_mat,0,1)./sqrt(size(self_mat,1));
    pt_m   = mean(pt_mat,1);   pt_sem   = std(pt_mat,0,1)  ./sqrt(size(pt_mat,1));
    ctrl_m = mean(ctrl_mat,1); ctrl_sem = std(ctrl_mat,0,1)./sqrt(size(ctrl_mat,1));
    bef_m  = mean(bef_mat,1);  bef_sem  = std(bef_mat,0,1) ./sqrt(size(bef_mat,1));
    aft_m  = mean(aft_mat,1);  aft_sem  = std(aft_mat,0,1) ./sqrt(size(aft_mat,1));

    [rpt, rpt_lo, rpt_hi] = bootstrap_ratio(pt_mat,   self_mat, n_boot);
    [rct, rct_lo, rct_hi] = bootstrap_ratio(ctrl_mat, self_mat, n_boot);
    [rpc, rpc_lo, rpc_hi] = bootstrap_ratio(pt_mat,   ctrl_mat, n_boot);   % PT / control
    [rd,  rd_lo,  rd_hi ] = bootstrap_ratio(aft_mat,  bef_mat,  n_boot);   % after / before

    P(end+1) = struct('name',populations{p,1},'color',populations{p,3},'k',kx, ...
        'self_m',self_m,'self_sem',self_sem,'pt_m',pt_m,'pt_sem',pt_sem, ...
        'ctrl_m',ctrl_m,'ctrl_sem',ctrl_sem, ...
        'rpt',rpt,'rpt_lo',rpt_lo,'rpt_hi',rpt_hi, ...
        'rct',rct,'rct_lo',rct_lo,'rct_hi',rct_hi, ...
        'rpc',rpc,'rpc_lo',rpc_lo,'rpc_hi',rpc_hi, ...
        'bef_m',bef_m,'bef_sem',bef_sem,'aft_m',aft_m,'aft_sem',aft_sem, ...
        'rd',rd,'rd_lo',rd_lo,'rd_hi',rd_hi);   %#ok<SAGROW>
    fprintf('Loaded %-6s : k = 1..%d\n', populations{p,1}, kx(end));
end

if isempty(P)
    error('No populations loaded. Check neural_base path (has the pc compute script been run?).');
end

%% Figure: two panels vs number of PCs
figure('Color','w','Position',[100 100 1180 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ----- Panel 1: raw accuracies -----
axA = nexttile; hold(axA,'on');
pop_handles = gobjects(1,numel(P));
for p = 1:numel(P)
    c = P(p).color; x = P(p).k;
    plot_band(x, P(p).self_m, P(p).self_sem, c, '-',  show_band);
    plot_band(x, P(p).pt_m,   P(p).pt_sem,   c, '--', show_band);
    plot_band(x, P(p).ctrl_m, P(p).ctrl_sem, c, ':',  show_band);
    pop_handles(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);  % colour proxy
end
yline(chance, '-', 'chance', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, ...
      'LabelHorizontalAlignment','left', 'FontSize',8, 'HandleVisibility','off');
xlabel('Number of PCs (k)'); ylabel('Decoding accuracy');
ylim([0 1]); grid on; box on;
title('Self, PT, and control decoding');

% linestyle key (black proxies, markers only here)
h_self = plot(nan, nan, '-k',  'Marker','s', 'MarkerFaceColor','k', 'LineWidth',1.5);
h_pt   = plot(nan, nan, '--k', 'Marker','o', 'LineWidth',1.5);
h_ctrl = plot(nan, nan, ':k',  'Marker','^', 'LineWidth',1.5);

% ----- Panel 2: normalized by self-decoding -----
nexttile; hold on;
for p = 1:numel(P)
    c = P(p).color; x = P(p).k;
    plot_ci(x, P(p).rpt, P(p).rpt_lo, P(p).rpt_hi, c, '--', show_band);
    plot_ci(x, P(p).rct, P(p).rct_lo, P(p).rct_hi, c, ':',  show_band);
end
yline(1, '-', 'self ceiling', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, ...
      'LabelHorizontalAlignment','left', 'FontSize',8, 'HandleVisibility','off');
xlabel('Number of PCs (k)'); ylabel('Accuracy / self-decoding');
yl = ylim; ylim([0 max(1.05, yl(2))]); grid on; box on;
title('PT and control normalized by self');

% ----- Shared legend -----
lg = legend(axA, [pop_handles, h_self, h_pt, h_ctrl], ...
       [{P.name}, {'self-decoding', pt_name, 'control (shuffle)'}]);
lg.Layout.Tile = 'east';
sgtitle('Procrustes transfer decoding vs PC dimension');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig_pc_self_pt_control.png'));

%% Figure 2: residual Procrustes shape distance vs number of PCs
figure('Color','w','Position',[100 100 1180 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ----- Panel 1: d_before and d_rot -----
axD = nexttile; hold(axD,'on');
pop_handles2 = gobjects(1,numel(P));
for p = 1:numel(P)
    c = P(p).color; x = P(p).k;
    plot_band(x, P(p).bef_m, P(p).bef_sem, c, '-',  show_band);
    plot_band(x, P(p).aft_m, P(p).aft_sem, c, '--', show_band);
    pop_handles2(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);
end
xlabel('Number of PCs (k)'); ylabel('Residual Procrustes distance');
grid on; box on;
title('Shape distance before vs after rotation');

h_bef = plot(nan, nan, '-k',  'Marker','s', 'MarkerFaceColor','k', 'LineWidth',1.5);
h_aft = plot(nan, nan, '--k', 'Marker','o', 'LineWidth',1.5);

% ----- Panel 2: after / before -----
nexttile; hold on;
for p = 1:numel(P)
    c = P(p).color; x = P(p).k;
    plot_ci(x, P(p).rd, P(p).rd_lo, P(p).rd_hi, c, '-', show_band);
end
yline(1, '-', 'no improvement', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, ...
      'LabelHorizontalAlignment','left', 'FontSize',8, 'HandleVisibility','off');
xlabel('Number of PCs (k)'); ylabel('d_{rot} / d_{before}  (after / before)');
yl = ylim; ylim([0 max(1.05, yl(2))]); grid on; box on;
title('Rotation''s residual reduction');

lg2 = legend(axD, [pop_handles2, h_bef, h_aft], ...
       [{P.name}, {'before rotation (d_{before})','after rotation (d_{rot})'}]);
lg2.Layout.Tile = 'east';
sgtitle('Procrustes shape distance vs PC dimension');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig_pc_procrustes_distance.png'));

%% Figure 3: ratio summary vs number of PCs (two subplots, ratios as percentages)
figure('Color','w','Position',[100 100 1180 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ----- Panel 1: ratios to self-decoding -----
axR = nexttile; hold(axR,'on');
pop_handles3 = gobjects(1,numel(P));
for p = 1:numel(P)
    c = P(p).color; x = P(p).k;
    plot_ci(x, 100*P(p).rpt, 100*P(p).rpt_lo, 100*P(p).rpt_hi, c, '--', show_band);   % PT / self
    plot_ci(x, 100*P(p).rct, 100*P(p).rct_lo, 100*P(p).rct_hi, c, ':',  show_band);   % control / self
    pop_handles3(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);
end
yline(100, '-', '100%', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, ...
      'LabelHorizontalAlignment','left', 'FontSize',8, 'HandleVisibility','off');
xlabel('Number of PCs (k)'); ylabel('Ratio to self-decoding (%)');
ylim([0 50])
grid on; box on;
title('PT and control, normalized by self');

% linestyle key (markers only here)
h_r1 = plot(nan, nan, '--k', 'Marker','o', 'LineWidth',1.5);
h_r2 = plot(nan, nan, ':k',  'Marker','^', 'LineWidth',1.5);
h_r3 = plot(nan, nan, '-k',  'Marker','s', 'MarkerFaceColor','k', 'LineWidth',1.5);

% ----- Panel 2: PT / control -----
nexttile; hold on;
for p = 1:numel(P)
    c = P(p).color; x = P(p).k;
    plot_ci(x, 100*P(p).rpc, 100*P(p).rpc_lo, 100*P(p).rpc_hi, c, '-', show_band);   % PT / control
end
yline(100, '-', '100%', 'Color',[0.6 0.6 0.6], 'LineWidth',0.8, ...
      'LabelHorizontalAlignment','left', 'FontSize',8, 'HandleVisibility','off');
xlabel('Number of PCs (k)'); ylabel('PT / control (%)');
grid on; box on;
title('PT relative to its shuffle null');

% ----- Shared legend -----
lg3 = legend(axR, [pop_handles3, h_r1, h_r2, h_r3], ...
       [{P.name}, {'PT / self', 'control / self', 'PT / control'}]);
lg3.Layout.Tile = 'east';
sgtitle('Ratio summary vs PC dimension');

% linestyle key (markers only here)
h_r1 = plot(nan, nan, '--k', 'Marker','o', 'LineWidth',1.5);
h_r2 = plot(nan, nan, ':k',  'Marker','^', 'LineWidth',1.5);
h_r3 = plot(nan, nan, '-k',  'Marker','s', 'MarkerFaceColor','k', 'LineWidth',1.5);
legend([pop_handles3, h_r1, h_r2, h_r3], ...
       [{P.name}, {'PT / self', 'control / self', 'PT / control'}], ...
       'Location','eastoutside');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig_pc_ratio_summary.png'));

%% ------------------------------------------------------------------------
function [kx, self_mat, pt_mat, ctrl_mat, bef_mat, aft_mat] = pool_population(pdir, pair_files, col_self, col_pt, col_ctrl, pcol_before, pcol_after)
% Load the six ordered-pair files once and pool their per-split rows for both the
% accuracy metrics (.acc) and the residual Procrustes distances (.pdist). Returns
% the PC list kx (1 x nSeq) and [nRows x nSeq] matrices (nRows = 6 pairs x
% repeats), columns sorted by ascending k.
nPairs = numel(pair_files);
res = cell(1,nPairs);
for k = 1:nPairs
    S = load(fullfile(pdir, [pair_files{k} '.mat']));
    res{k} = S.(pair_files{k});   % 1 x nSeq cell of structs
end

nSeq  = numel(res{1});
nrep  = size(res{1}{1}.acc, 1);
nRows = nPairs * nrep;

kx       = zeros(1,nSeq);
self_mat = zeros(nRows,nSeq);
pt_mat   = zeros(nRows,nSeq);
ctrl_mat = zeros(nRows,nSeq);
bef_mat  = zeros(nRows,nSeq);
aft_mat  = zeros(nRows,nSeq);

for s = 1:nSeq
    kx(s) = res{1}{s}.pc_num;
    u = 0;
    for k = 1:nPairs
        c    = res{k}{s};
        rows = u + (1:size(c.acc,1));
        self_mat(rows,s) = c.acc(:,col_self);
        pt_mat(rows,s)   = c.acc(:,col_pt);
        ctrl_mat(rows,s) = c.acc(:,col_ctrl);
        bef_mat(rows,s)  = c.pdist(:,pcol_before);
        aft_mat(rows,s)  = c.pdist(:,pcol_after);
        u = u + size(c.acc,1);
    end
end

[kx, order] = sort(kx);
self_mat = self_mat(:,order);
pt_mat   = pt_mat(:,order);
ctrl_mat = ctrl_mat(:,order);
bef_mat  = bef_mat(:,order);
aft_mat  = aft_mat(:,order);
end

function [r, lo, hi] = bootstrap_ratio(num_mat, den_mat, nboot)
% Ratio-of-means with a percentile bootstrap over rows (pooled resampling rows),
% same estimator as the pr script. num_mat/den_mat are [nRows x nSeq], paired by
% row; the same resampled indices are used for numerator and denominator.
[nRows, nSeq] = size(num_mat);
r  = zeros(1,nSeq); lo = zeros(1,nSeq); hi = zeros(1,nSeq);
for s = 1:nSeq
    nu = num_mat(:,s); de = den_mat(:,s);
    r(s) = mean(nu) / mean(de);
    bs = zeros(nboot,1);
    for b = 1:nboot
        idx = randi(nRows, nRows, 1);
        bs(b) = mean(nu(idx)) / mean(de(idx));
    end
    lo(s) = prctile(bs, 2.5);
    hi(s) = prctile(bs, 97.5);
end
end

function plot_band(x, m, sem, c, ls, show_band)
% mean line + optional translucent +/-1 SEM band. Data lines carry no markers
% (the key supplies them) and are hidden from the legend.
if show_band
    fill([x fliplr(x)], [m-sem fliplr(m+sem)], c, ...
        'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
end
plot(x, m, ls, 'Color',c, 'LineWidth',1.8, 'HandleVisibility','off');
end

function plot_ci(x, r, lo, hi, c, ls, show_band)
% ratio line + optional translucent bootstrap-CI band.
if show_band
    fill([x fliplr(x)], [lo fliplr(hi)], c, ...
        'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
end
plot(x, r, ls, 'Color',c, 'LineWidth',1.8, 'HandleVisibility','off');
end
