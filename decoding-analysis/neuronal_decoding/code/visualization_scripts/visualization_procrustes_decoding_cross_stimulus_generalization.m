%{
Filename: visualization_procrustes_decoding_cross_stimulus_generalization.m
Author:   Zitong Wang
Date:     2026-07-30

Description:
    First-look figures for the cross-stimulus (CCGP-style) generalization
    analysis produced by procrustes_decoding_cross_stimulus_generalization.m.

    One figure per monkey/area. Each figure has THREE panels, one per unordered
    rendering pair -- AC<->EC, EC<->EX, AC<->EX -- so the three pairs are shown
    individually and NOT averaged into a single six-pair curve. Within a panel
    the two ordered directions of that pair (e.g. ac->ec and ec->ac) are pooled;
    to display a single direction instead, drop one file stem from that pair's
    entry in the `pairs` table below (the framework handles 1 or 2 directions).

    Axes are deliberately simple: x = number of units, y = decoding accuracy.
    Curves (columns of .acc):
        self_decode  (1) within-cue ceiling
        pt_gen       (3) MAIN -- rotation fit on 40, applied to the held-out 10
        pt_ceiling   (4) in-sample PT ceiling (rotation fit on the 10)
        pt_floor     (5) self-consistent-shuffle null
        chance       (7)
    no_transform (2) and rand_rot_40 (6) are defined in the style table and can
    be shown by adding their column to `measures_to_plot`.

    Aggregation: for each neuron-count, per-PARTITION means are formed first
    (partition = the CCGP resampling unit) and pooled across the pair's
    directions; the plotted value is the mean over partitions and the error bar
    is +/-1 SEM over partitions. This is a deliberate departure from the sibling
    _pr script's pool-all-rows SEM: it makes partition the unit of resampling and
    avoids treating the (by construction) replicated self_decode rows -- constant
    across the trial loop within a (partition, neuron) cell -- as independent.
    (The point estimate is unchanged, since partitions carry equal row counts.)

Path note:
    Expected to live in neuronal_decoding/code/visualization_scripts/.
    Results tree: neuronal_decoding/results/decoding_outputs/...  (../../results),
    anchored to this file's own location via mfilename.
%}

%% Configuration
clc; clear;
close all
% Resolve this script's folder (fall back to pwd when run cell-by-cell).
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir) || contains(script_dir, fullfile('AppData','Local','Temp'))
    script_dir = pwd;
end

% --- Hold-out level to visualize (selects the HoldStim<n> results folder) ---
n_stim_hold = 2;   % must match the compute run you want to view

results_base = fullfile(script_dir, '..','..','results','decoding_outputs', ...
    'procrustes_decoding_cross_stimulus_generalization_results', ...
    sprintf('HoldStim%d', n_stim_hold));

fig_dir = fullfile(script_dir, '..','..','results','figures','cross_stimulus_generalization', ...
    sprintf('HoldStim%d', n_stim_hold));
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- Areas to plot: {monkey, area} per row. One figure each. ---
areas_to_plot = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};

% --- Unordered pairs: {panel title, {ordered-direction file stems}} ---
%     Two stems  -> the panel pools both directions.
%     One stem   -> single-direction panel (just delete the other stem).
pairs = {
    'AC \leftrightarrow EC', {'acec_results','ecac_results'};
    'EC \leftrightarrow EX', {'ecex_results','exec_results'};
    'AC \leftrightarrow EX', {'acex_results','exac_results'};
};

% --- Measure style table, indexed by .acc column (1..7) ---
% Labels track the hold-out level: fit set = 50 - n_stim_hold, hold-out = n_stim_hold.
n_fit = 50 - n_stim_hold;
m_name  = {'self-decoding', 'no-transform', ...
           sprintf('PT-gen (fit %d)', n_fit), ...
           sprintf('PT-ceiling (fit %d)', n_stim_hold), ...
           'PT-floor (shuffle)', ...
           sprintf('rand-rot (fit %d)', n_fit), ...
           'chance'};
m_color = {[0.50 0.50 0.50],[0.70 0.70 0.70],[0.00 0.45 0.74],[0.30 0.75 0.93], ...
           [0.85 0.33 0.10],[0.95 0.55 0.45],[0.00 0.00 0.00]};
m_style = {'--',':','-','--','-','--',':'};
m_mark  = {'s','none','o','^','v','none','none'};
m_lw    = [1.5 1.2 2.2 1.5 1.5 1.2 1.2];

% Which measures to draw (subset of columns 1..7). All seven shown here:
% 1 self-decoding, 2 no-transform, 3 PT-gen, 4 PT-ceiling, 5 PT-floor,
% 6 rand-rot (fit 40), 7 chance.
measures_to_plot = [1 2 3 4 5 6 7];

show_err = true;   % +/-1 SEM error bars over the partition distribution

%% Build one figure per area
for a = 1:size(areas_to_plot,1)
    monkey = areas_to_plot{a,1};
    vp     = areas_to_plot{a,2};
    pdir   = fullfile(results_base, monkey, vp);
    if ~exist(pdir,'dir')
        warning('Missing results dir, skipping %s %s: %s', monkey, vp, pdir);
        continue;
    end

    figure('Color','w','Position',[100 100 1250 430]);
    tiledlayout(1, size(pairs,1), 'TileSpacing','compact','Padding','compact');

    leg_h = gobjects(1,numel(measures_to_plot));
    leg_names = cell(1,numel(measures_to_plot));

    for pr = 1:size(pairs,1)
        nexttile; hold on;

        dir_structs = load_pair(pdir, pairs{pr,2});   % 1xD cell of result cells

        for ci = 1:numel(measures_to_plot)
            col = measures_to_plot(ci);
            [N, m, sem] = pair_partition_stats(dir_structs, col);
            if show_err
                h = errorbar(N, m, sem, m_style{col}, 'Color',m_color{col}, ...
                    'Marker',m_mark{col}, 'MarkerFaceColor',m_color{col}, ...
                    'MarkerSize',4, 'LineWidth',m_lw(col), 'CapSize',3);
            else
                h = plot(N, m, m_style{col}, 'Color',m_color{col}, ...
                    'Marker',m_mark{col}, 'MarkerFaceColor',m_color{col}, ...
                    'MarkerSize',4, 'LineWidth',m_lw(col));
            end
            if pr == 1
                leg_h(ci) = h;
                leg_names{ci} = m_name{col};
            end
        end

        ylim([0 1]); grid on; box on;
        xlabel('Number of units');
        if pr == 1, ylabel('Decoding accuracy'); end
        title(pairs{pr,1});
    end

    lg = legend(leg_h, leg_names);
    lg.Layout.Tile = 'east';
    sgtitle(sprintf('%s %s  --  cross-stimulus generalization (hold %d)', monkey, vp, n_stim_hold));

    % saveas(gcf, fullfile(fig_dir, sprintf('cross_stim_gen_HoldStim%d_%s_%s.png', n_stim_hold, monkey, vp)));
end

%% ---------- helpers ----------
function res = load_pair(pdir, files)
% files: 1xD cell of file/variable stems. Returns 1xD cell, each a 1xnSeq cell
% of per-neuron-count structs (the ordered directions of one unordered pair).
D = numel(files);
res = cell(1,D);
for d = 1:D
    S = load(fullfile(pdir, [files{d} '.mat']));
    res{d} = S.(files{d});
end
end

function [N, m, sem] = pair_partition_stats(dir_structs, col)
% For each neuron-count, pool per-PARTITION means of column `col` across the
% pair's directions, then return mean and +/-1 SEM over partitions. Sorted by N.
D    = numel(dir_structs);
nSeq = numel(dir_structs{1});
N   = zeros(1,nSeq);
m   = zeros(1,nSeq);
sem = zeros(1,nSeq);
for s = 1:nSeq
    N(s) = dir_structs{1}{s}.neuron_num;
    part_means = [];
    for d = 1:D
        c = dir_structs{d}{s};
        % per-partition mean of this measure (partition = CCGP resampling unit)
        pm = accumarray(c.partition_id, c.acc(:,col), [], @mean);
        part_means = [part_means; pm];
    end
    m(s)   = mean(part_means);
    sem(s) = std(part_means) / sqrt(numel(part_means));
end
[N, ord] = sort(N);
m   = m(ord);
sem = sem(ord);
end
