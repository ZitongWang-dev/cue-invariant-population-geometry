%{
Filename: visualization_procrustes_meancentroid_cross_stimulus_generalization.m
Author:   Zitong Wang
Date:     2026-08-06

Description:
    Figures for the MEAN-LEVEL nearest-centroid generalization analysis produced
    by procrustes_meancentroid_cross_stimulus_generalization.m.

    One figure per monkey/area. Each figure is a 2 x 3 grid:
        columns = the three unordered rendering pairs (AC<->EC, EC<->EX, AC<->EX),
                  each pooling its two ordered directions;
        rows    = the two reference sets --
                  ROW 1 GLOBAL  : transformed held-out source means decoded against
                                  all 50 target means (chance = 1/50);
                  ROW 2 WITHIN  : decoded against only the held-out target means
                                  (chance = 1/n_stim_hold).
    Within each panel, three curves vs number of units:
        gen           rotation fit on F (true correspondence) -- the main measure
        no-transform  raw source means, no rotation (the same-neuron baseline)
        rand          rotation fit on F with shuffled correspondence (floor)
    A dotted line marks analytic chance for that reference set.

    The key read is gen vs no-transform (does the rotation ADD cross-cue
    alignment beyond what the raw same-neuron responses already give?), with rand
    and the chance line marking the floor.

    Aggregation: per-PARTITION means are formed first (partition = the resampling
    unit) and pooled across the pair's two directions; plotted value is the mean
    over partitions, error bar +/-1 SEM over partitions. Same convention as the
    trial-level cross-stimulus viz.

Path note:
    Expected in neuronal_decoding/code/visualization_scripts/. Results tree
    anchored to this file via mfilename.
%}

%% Configuration
clc; clear;

% Resolve this script's folder (fall back to pwd when run cell-by-cell).
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir) || contains(script_dir, fullfile('AppData','Local','Temp'))
    script_dir = pwd;
end

% --- Hold-out level to visualize (selects the HoldStim<n> results folder) ---
n_stim_hold = 2;   % must match the compute run you want to view

results_base = fullfile(script_dir, '..','..','results','decoding_outputs', ...
    'procrustes_meancentroid_cross_stimulus_generalization_results', ...
    sprintf('HoldStim%d', n_stim_hold));

fig_dir = fullfile(script_dir, '..','..','results','figures', ...
    'meancentroid_cross_stimulus_generalization', sprintf('HoldStim%d', n_stim_hold));
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- Areas to plot: {monkey, area} per row. One figure each. ---
areas_to_plot = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};

% --- Unordered pairs (columns): {panel title, {ordered-direction file stems}} ---
%     Two stems -> pooled; one stem -> single-direction (delete the other).
pairs = {
    'AC \leftrightarrow EC', {'acec_results','ecac_results'};
    'EC \leftrightarrow EX', {'ecex_results','exec_results'};
    'AC \leftrightarrow EX', {'acex_results','exac_results'};
};

% --- Reference sets (rows): {row label, [gen nt rand] columns in .acc, chance} ---
ref_sets = {
    'Global', [1 2 3], 1/50;
    'Within', [4 5 6], 1/n_stim_hold;
};

% --- Measure styles, in the gen/no-transform/rand order used within each ref set ---
n_fit = 50 - n_stim_hold;
meas_name  = {sprintf('gen (fit %d)', n_fit), 'no-transform', 'rand (shuffle)'};
meas_color = {[0.00 0.45 0.74], [0.50 0.50 0.50], [0.85 0.33 0.10]};
meas_style = {'-', '--', '-'};
meas_mark  = {'o', 's', 'v'};
meas_lw    = [2.2 1.5 1.5];

%% Build one figure per area
for a = 1:size(areas_to_plot,1)
    monkey = areas_to_plot{a,1};
    vp     = areas_to_plot{a,2};
    pdir   = fullfile(results_base, monkey, vp);
    if ~exist(pdir,'dir')
        warning('Missing results dir, skipping %s %s: %s', monkey, vp, pdir);
        continue;
    end

    figure('Color','w','Position',[80 80 1250 720]);
    tiledlayout(size(ref_sets,1), size(pairs,1), 'TileSpacing','compact','Padding','compact');

    leg_h     = gobjects(1, numel(meas_name)+1);   % 3 measures + chance
    leg_names = [meas_name, {'chance'}];
    got_leg   = false;

    for r = 1:size(ref_sets,1)
        cols   = ref_sets{r,2};
        chance = ref_sets{r,3};
        for pc = 1:size(pairs,1)
            nexttile; hold on;

            dir_structs = load_pair(pdir, pairs{pc,2});
            for mi = 1:numel(cols)
                [N, m, sem] = pair_partition_stats(dir_structs, cols(mi));
                h = errorbar(N, m, sem, meas_style{mi}, 'Color',meas_color{mi}, ...
                    'Marker',meas_mark{mi}, 'MarkerFaceColor',meas_color{mi}, ...
                    'MarkerSize',4, 'LineWidth',meas_lw(mi), 'CapSize',3);
                if ~got_leg, leg_h(mi) = h; end
            end

            yline(chance, ':k', 'LineWidth',1.2);   % analytic chance for this ref set
            if ~got_leg
                leg_h(end) = plot(nan, nan, ':k', 'LineWidth',1.2);  % legend proxy
                got_leg = true;
            end

            ylim([0 1]); grid on; box on;
            if r == 1,                title(pairs{pc,1}); ylim([0 0.2]); end
            if r == size(ref_sets,1), xlabel('Number of units'); end
            if pc == 1,               ylabel(sprintf('%s  --  accuracy', ref_sets{r,1})); end
        end
    end

    lg = legend(leg_h, leg_names);
    lg.Layout.Tile = 'east';
    sgtitle(sprintf('%s %s  --  mean nearest-centroid generalization (hold %d)', monkey, vp, n_stim_hold));

    saveas(gcf, fullfile(fig_dir, sprintf('meancentroid_HoldStim%d_%s_%s.png', n_stim_hold, monkey, vp)));
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
        pm = accumarray(c.partition_id, c.acc(:,col), [], @mean);  % per-partition mean
        part_means = [part_means; pm];
    end
    m(s)   = mean(part_means);
    sem(s) = std(part_means) / sqrt(numel(part_means));
end
[N, ord] = sort(N);
m   = m(ord);
sem = sem(ord);
end
