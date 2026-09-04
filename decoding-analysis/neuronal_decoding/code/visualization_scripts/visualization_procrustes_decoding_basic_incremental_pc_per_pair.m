%{
Filename: visualization_procrustes_decoding_basic_incremental_pc_per_pair.m
Author:   Zitong Wang
Date:     2026-08

Description:
    Per-pair Procrustes residual figures for the PC-dimension transfer analysis
    (procrustes_decoding_basic_incremental_pc.m), written to answer Reviewer
    point 10. The companion visualization_procrustes_decoding_basic_incremental_pc.m
    POOLS the residual across all six ordered rendering pairs; this script instead
    UN-POOLS it into the three unordered cue pairs and compares them directly:

        EC-AC  =  acec_results + ecac_results
        EC-EX  =  ecex_results + exec_results
        AC-EX  =  acex_results + exac_results   (the pair the reviewer flags)

    Reviewer point 10 (the concern). The AC rendering is the first PC of the same
    patches that define EX, so AC-EX alignment could be inflated by shared stimulus
    construction rather than by genuine cue-invariant coding. If it were inflated,
    AC-EX should be the TIGHTEST pair: the SMALLEST Procrustes residual. This
    script plots the residual per pair so that prediction can be checked directly.

    Quantities (from .pdist, all sharing MATLAB's stim2-target denominator, so they
    are directly comparable across pairs; smaller = tighter alignment):
        d_rot   (.pdist col 2)  residual after rotation only  (b = 1, T only)
        d_full  (.pdist col 3)  MATLAB procrustes d (optimal b, T, c) -- the floor
    d_before (col 1) is the no-rotation baseline and is not plotted here.

    FIGURE 1 -- d_rot vs k, one tile per population (2x2), three coloured pair
                curves with +/-1 SEM bands. x-axis = number of retained PCs k.
    FIGURE 2 -- d_full vs k, same layout.
    FIGURE 3 -- collapsed summary: grouped bars per population (left panel d_rot,
                right panel d_full), bar height = mean over k, error bar = the
                typical per-k sampling SEM (mean across k of the per-k SEM, so the
                bar's uncertainty reflects sampling noise, NOT k-dependence). Makes
                the across-pair ordering legible at a glance.

    Caveat (inherited from the compute script): each cue is first projected onto its
    OWN top-k mean-PCA basis, so these are residual shape distances WITHIN that
    basis, not the full cross-cue rotation magnitude on the untruncated manifold.

Conventions:
    colour = cue pair (AC-EX drawn in red, the pair under suspicion); one tile per
    population; +/-1 SEM over the pooled resampling rows (2 ordered files x
    trial_sample_repeat rows per k). Path anchoring matches the pooled script.

Path note:
    Expected to live in neuronal_decoding/code/visualization_scripts/.
    Neural results tree: ../../results/decoding_outputs/...  (anchored via mfilename).
    Figures written to:  ../../results/figures/incremental_pc/.
%}

%% Configuration
clc; clear;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir) || contains(script_dir, fullfile('AppData','Local','Temp'))
    script_dir = pwd;
end

neural_base = fullfile(script_dir, '..','..','results','decoding_outputs', ...
    'procrustes_decoding_basic_incremental_pc_results');

% --- Populations: {display name, directory holding its 6 ordered-pair files} ---
populations = {
    'FR V1', fullfile(neural_base,'FR','V1');
    'FR V2', fullfile(neural_base,'FR','V2');
    'KO V1', fullfile(neural_base,'KO','V1');
    'KO V2', fullfile(neural_base,'KO','V2');
};

% --- Unordered cue pairs: {display name, {two ordered files}, colour} ---
%     Each unordered pair pools the two orderings so the residual is direction-free.
pair_groups = {
    'EC-AC', {'acec_results','ecac_results'}, [0.00 0.45 0.74];   % blue
    'EC-EX', {'ecex_results','exec_results'}, [0.20 0.62 0.30];   % green
    'AC-EX', {'acex_results','exac_results'}, [0.83 0.15 0.15];   % red (reviewer's pair)
};

% --- Columns inside .pdist ---
pcol_rot  = 2;   % d_rot   (rotation only, b = 1)
pcol_full = 3;   % d_full  (optimal b, T, c)

% --- Display options ---
show_band = true;    % shaded +/-1 SEM bands on the k-resolved curves

fig_dir = fullfile(script_dir, '..','..','results','figures','incremental_pc');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

nPop  = size(populations,1);
nPair = size(pair_groups,1);

%% Load: per population, per unordered pair -> per-k mean/SEM for d_rot and d_full
% D(p,u) holds everything needed to plot population p, pair u.
D = struct('k',{},'rot_m',{},'rot_sem',{},'full_m',{},'full_sem',{});
have_pop = false(1,nPop);

for p = 1:nPop
    pdir = populations{p,2};
    if ~exist(pdir,'dir')
        warning('Missing directory, skipping %s: %s', populations{p,1}, pdir);
        continue;
    end
    ok = true;
    for u = 1:nPair
        try
            [kx, Mrot]  = load_pair_metric(pdir, pair_groups{u,2}, pcol_rot);
            [~,  Mfull] = load_pair_metric(pdir, pair_groups{u,2}, pcol_full);
        catch ME
            warning('Could not load %s / %s (%s); skipping population.', ...
                populations{p,1}, pair_groups{u,1}, ME.message);
            ok = false; break;
        end
        D(p,u).k       = kx;
        D(p,u).rot_m   = mean(Mrot,1);
        D(p,u).rot_sem = std(Mrot,0,1)./sqrt(size(Mrot,1));
        D(p,u).full_m  = mean(Mfull,1);
        D(p,u).full_sem= std(Mfull,0,1)./sqrt(size(Mfull,1));
    end
    have_pop(p) = ok;
    if ok, fprintf('Loaded %-6s : k = 1..%d\n', populations{p,1}, D(p,1).k(end)); end
end

if ~any(have_pop)
    error('No populations loaded. Check neural_base (has the pc compute script run?).');
end

%% Console summary: per-pair ordering (the reviewer's exact test)
fprintf('\n=== Mean residual across k (smaller = tighter alignment) ===\n');
fprintf('%-6s | %-7s | %8s %8s %8s | loosest / tightest\n','pop','metric','EC-AC','EC-EX','AC-EX');
for p = 1:nPop
    if ~have_pop(p), continue; end
    for which = {'rot','full'}
        w = which{1};
        vals = zeros(1,nPair);
        for u = 1:nPair, vals(u) = mean(D(p,u).([w '_m'])); end
        [~,imax] = max(vals); [~,imin] = min(vals);
        fprintf('%-6s | d_%-5s | %8.4f %8.4f %8.4f | loosest=%s tightest=%s\n', ...
            populations{p,1}, w, vals(1), vals(2), vals(3), ...
            pair_groups{imax,1}, pair_groups{imin,1});
    end
end
fprintf(['Inflation hypothesis (point 10) predicts AC-EX = TIGHTEST (smallest). ' ...
         'Observe whether AC-EX is ever the tightest above.\n\n']);

%% Figure 1 -- d_rot per pair, one tile per population
plot_metric_grid(D, have_pop, populations, pair_groups, 'rot', show_band, ...
    'Residual after rotation (d_{rot})', ...
    'Per-pair Procrustes residual d_{rot} vs PC dimension');
saveas(gcf, fullfile(fig_dir, 'fig_pc_procrustes_distance_per_pair_drot.png'));

%% Figure 2 -- d_full per pair, one tile per population
plot_metric_grid(D, have_pop, populations, pair_groups, 'full', show_band, ...
    'Full Procrustes residual (d_{full})', ...
    'Per-pair Procrustes residual d_{full} vs PC dimension');
saveas(gcf, fullfile(fig_dir, 'fig_pc_procrustes_distance_per_pair_dfull.png'));

%% Figure 3 -- collapsed summary: grouped bars per population
% Bar height = mean over k; error bar = mean over k of the per-k SEM (typical
% sampling noise, decoupled from the systematic k-dependence).
rot_bar  = nan(nPop,nPair);  rot_err  = nan(nPop,nPair);
full_bar = nan(nPop,nPair);  full_err = nan(nPop,nPair);
for p = 1:nPop
    if ~have_pop(p), continue; end
    for u = 1:nPair
        rot_bar(p,u)  = mean(D(p,u).rot_m);   rot_err(p,u)  = mean(D(p,u).rot_sem);
        full_bar(p,u) = mean(D(p,u).full_m);  full_err(p,u) = mean(D(p,u).full_sem);
    end
end

figure('Color','w','Position',[100 100 1180 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
pair_colors = cell2mat(pair_groups(:,3));

axB1 = nexttile;
grouped_bars(rot_bar, rot_err, pair_colors, populations(:,1));
ylabel('d_{rot}  (mean over k)'); title('Rotation-only residual');
axB2 = nexttile;
grouped_bars(full_bar, full_err, pair_colors, populations(:,1));
ylabel('d_{full}  (mean over k)'); title('Full Procrustes residual');

% pair legend (proxies)
hbar = gobjects(1,nPair);
hold(axB1,'on');
for u = 1:nPair
    hbar(u) = bar(axB1, nan, nan, 'FaceColor', pair_groups{u,3});
end
lg = legend(axB1, hbar, pair_groups(:,1), 'Location','northwest');
title(lg,'cue pair');
sgtitle(['Per-pair residual summary  (AC-EX is never the tightest ' ...
         '-> no shared-construction inflation)']);
saveas(gcf, fullfile(fig_dir, 'fig_pc_procrustes_distance_per_pair_summary.png'));

fprintf('Figures written to %s\n', fig_dir);

%% ------------------------------------------------------------------------
function [kx, M] = load_pair_metric(pdir, ordered_files, pcol)
% Pool the two ordered files of one unordered pair. Returns the PC list kx
% (1 x nSeq, ascending) and M = [nRows x nSeq] of the requested .pdist column,
% nRows = (#ordered files) x trial_sample_repeat, columns sorted by ascending k.
nf = numel(ordered_files);
res = cell(1,nf);
for f = 1:nf
    S = load(fullfile(pdir, [ordered_files{f} '.mat']));
    res{f} = S.(ordered_files{f});   % 1 x nSeq cell of structs
end
nSeq = numel(res{1});
nrep = size(res{1}{1}.pdist, 1);
kx = zeros(1,nSeq);
M  = zeros(nf*nrep, nSeq);
for s = 1:nSeq
    kx(s) = res{1}{s}.pc_num;
    u = 0;
    for f = 1:nf
        rows = u + (1:nrep);
        M(rows,s) = res{f}{s}.pdist(:,pcol);
        u = u + nrep;
    end
end
[kx, order] = sort(kx);
M = M(:,order);
end

function plot_metric_grid(D, have_pop, populations, pair_groups, which, show_band, ylab, suptitle_str)
% 2x2 grid, one tile per population; three coloured pair curves + SEM bands.
nPop  = size(populations,1);
nPair = size(pair_groups,1);
figure('Color','w','Position',[100 100 1100 780]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax0 = gobjects(1,1);
for p = 1:nPop
    ax = nexttile; hold(ax,'on');
    if p == 1, ax0 = ax; end
    if ~have_pop(p)
        title(ax, sprintf('%s (missing)', populations{p,1})); axis(ax,'off'); continue;
    end
    for u = 1:nPair
        c = pair_groups{u,3};
        x = D(p,u).k;
        m = D(p,u).([which '_m']);  sem = D(p,u).([which '_sem']);
        if show_band
            fill([x fliplr(x)], [m-sem fliplr(m+sem)], c, ...
                'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off','Parent',ax);
        end
        plot(ax, x, m, '-', 'Color',c, 'LineWidth',1.8);
    end
    xlabel(ax,'Number of PCs (k)'); ylabel(ax, ylab);
    title(ax, populations{p,1}); grid(ax,'on'); box(ax,'on');
end
if isgraphics(ax0)
    lg = legend(ax0, pair_groups(:,1), 'Location','best');
    title(lg,'cue pair');
end
sgtitle(suptitle_str);
end

function grouped_bars(vals, errs, colors, xlabels)
% Grouped bar chart (rows = groups/populations, cols = pairs) with SEM error bars.
hb = bar(vals, 'grouped'); hold on;
nPair = size(vals,2);
for u = 1:nPair
    hb(u).FaceColor = colors(u,:);
    x = hb(u).XEndPoints;
    errorbar(x, vals(:,u), errs(:,u), 'k', 'linestyle','none', 'LineWidth',0.9, ...
        'HandleVisibility','off');
end
set(gca,'XTick',1:numel(xlabels),'XTickLabel',xlabels);
grid on; box on; ylim([0 max(vals(:)+errs(:),[],'omitnan')*1.12]);
end
