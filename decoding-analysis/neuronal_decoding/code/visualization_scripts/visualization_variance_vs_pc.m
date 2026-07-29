%{
Filename: visualization_variance_vs_pc.m
Author:   Zitong Wang
Date:     2026-07

Description:
    Signal-variance spectrum vs number of principal components, k, for each cue.
    This quantity is NOT saved by procrustes_decoding_basic_incremental_pc.m
    (that script keeps only the PC directions, not the singular values), and it
    cannot be recovered from the results .mat files -- so this script recomputes
    it directly from the raw neural data, using the SAME construction the compute
    script's mean_basis uses:

        per cue:  z-score per unit over all 500 trials
               -> take_average(.,10)  =>  50 x N centered mean manifold
               -> svd  =>  lambda_i = s_i^2
               -> variance ratio      = lambda / sum(lambda)
               -> cumulative variance = cumsum(variance ratio)

    This is the fraction of the across-stimulus (signal) variance captured by the
    top-k PCs -- i.e. exactly the geometry whose top-k the decoding retains. It
    uses ALL trials (split-independent), so it is the clean reference spectrum;
    the decoding builds its stim1 basis from training trials only, so per-split
    spectra wobble slightly around this, but the shape is essentially identical.

      PANEL 1: cumulative variance explained vs k (with 80% / 90% references).
      PANEL 2: marginal (per-PC) variance ratio vs k -- the scree curve.

    All four populations overlaid, colour = population. The three cues (AC/EC/EX)
    are summarised per population as the mean curve with a min-max band across
    cues (set plot_cues_separately = true to draw them as individual thin lines).

Dependencies:
    multiclass_svmloader_PT.m (on the MATLAB path) and the raw data tree
    neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat.

Path note:
    Expected in neuronal_decoding/code/visualization_scripts/; raw data resolved
    as ../../neuronal_data (same anchor the compute script uses).
%}

%% Configuration
clc; clear;

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir) || contains(script_dir, fullfile('AppData','Local','Temp'))
    script_dir = pwd;
end

raw_base = fullfile(script_dir, '..','..','neuronal_data');

populations = {
    'FR V1', 'FR','V1', [0.00 0.45 0.74];
    'FR V2', 'FR','V2', [0.30 0.75 0.93];
    'KO V1', 'KO','V1', [0.72 0.10 0.10];
    'KO V2', 'KO','V2', [0.95 0.50 0.50];
};

cues = {'AC',1; 'EC',2; 'EX',3};    % name, index into three_stim_array
timewindow = [330 630];             % must match the compute script
kmax = 49;                          % cap (signal manifold rank <= 49)
plot_cues_separately = false;       % false: mean +/- min-max band across cues

fig_dir = fullfile(script_dir, '..','..','results','figures','incremental_pc');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Compute per-cue spectra for each population
nPop = size(populations,1);
V = struct('name',{},'color',{},'k',{}, ...
           'cum',{},'cum_lo',{},'cum_hi',{},'mar',{},'mar_lo',{},'mar_hi',{}, ...
           'cum_cues',{},'mar_cues',{});

for p = 1:nPop
    monkey = populations{p,2}; vp = populations{p,3};
    f = fullfile(raw_base, monkey, vp, sprintf('%s_%s_allstim.mat', monkey, vp));
    if ~exist(f,'file')
        warning('Missing raw data, skipping %s: %s', populations{p,1}, f);
        continue;
    end
    S = load(f, 'three_stim_array');
    spike_data = S.three_stim_array;

    nCues = size(cues,1);
    cum_cell = cell(1,nCues);
    mar_cell = cell(1,nCues);
    for ci = 1:nCues
        data = multiclass_svmloader_PT(spike_data{cues{ci,2}}, timewindow);  % 500 x N
        data = zscore(data);
        means = take_average(data, 10);           % 50 x N
        Xc = means - mean(means,1);
        s  = svd(Xc);                             % descending singular values
        lambda = s.^2;
        vr  = lambda / sum(lambda);               % marginal variance ratio
        cvr = cumsum(vr);                         % cumulative
        kk  = min(kmax, numel(vr));
        mar_cell{ci} = vr(1:kk).';
        cum_cell{ci} = cvr(1:kk).';
    end

    % common length across cues (same N -> same length; guard anyway)
    L = min(cellfun(@numel, cum_cell));
    cum_mat = cell2mat(cellfun(@(v) v(1:L), cum_cell.', 'UniformOutput', false));  % nCues x L
    mar_mat = cell2mat(cellfun(@(v) v(1:L), mar_cell.', 'UniformOutput', false));

    V(end+1) = struct('name',populations{p,1},'color',populations{p,4},'k',1:L, ...
        'cum',mean(cum_mat,1),'cum_lo',min(cum_mat,[],1),'cum_hi',max(cum_mat,[],1), ...
        'mar',mean(mar_mat,1),'mar_lo',min(mar_mat,[],1),'mar_hi',max(mar_mat,[],1), ...
        'cum_cues',{cum_mat},'mar_cues',{mar_mat});   %#ok<SAGROW>
    fprintf('Computed %-6s : N-derived spectrum over k = 1..%d\n', populations{p,1}, L);
end

if isempty(V)
    error('No populations computed. Check raw_base path and multiclass_svmloader_PT on path.');
end

%% Figure: cumulative (left) and marginal scree (right)
figure('Color','w','Position',[100 100 1180 480]);
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% ----- Panel 1: cumulative variance -----
axA = nexttile; hold(axA,'on');
pop_handles = gobjects(1,numel(V));
for p = 1:numel(V)
    c = V(p).color; x = V(p).k;
    if plot_cues_separately
        plot(x, V(p).cum_cues.', '-', 'Color',[c 0.35], 'LineWidth',0.8, 'HandleVisibility','off');
    else
        fill([x fliplr(x)], [V(p).cum_lo fliplr(V(p).cum_hi)], c, ...
            'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
    end
    plot(x, V(p).cum, '-', 'Color',c, 'LineWidth',1.8, 'HandleVisibility','off');
    pop_handles(p) = plot(nan, nan, '-', 'Color',c, 'LineWidth',2.5);
end
yline(0.90, ':', '90%', 'Color',[0.6 0.6 0.6], 'FontSize',8, 'HandleVisibility','off');
yline(0.80, ':', '80%', 'Color',[0.6 0.6 0.6], 'FontSize',8, 'HandleVisibility','off');
xlabel('Number of PCs (k)'); ylabel('Cumulative signal variance explained');
ylim([0 1]); grid on; box on;
title('Cumulative variance vs PC dimension');

% ----- Panel 2: marginal scree -----
nexttile; hold on;
for p = 1:numel(V)
    c = V(p).color; x = V(p).k;
    if plot_cues_separately
        plot(x, V(p).mar_cues.', '-', 'Color',[c 0.35], 'LineWidth',0.8, 'HandleVisibility','off');
    else
        fill([x fliplr(x)], [V(p).mar_lo fliplr(V(p).mar_hi)], c, ...
            'FaceAlpha',0.12, 'EdgeColor','none', 'HandleVisibility','off');
    end
    plot(x, V(p).mar, '-', 'Color',c, 'LineWidth',1.8, 'HandleVisibility','off');
end
xlabel('Number of PCs (k)'); ylabel('Variance ratio per PC (scree)');
set(gca,'YScale','log'); grid on; box on;
title('Per-PC variance (scree)');

lg = legend(axA, pop_handles, {V.name});
lg.Layout.Tile = 'east';
sgtitle('Signal-variance spectrum of the 50-condition mean manifold');

%% Save (uncomment for production)
% saveas(gcf, fullfile(fig_dir, 'fig_pc_signal_variance.png'));

%% ------------------------------------------------------------------------
function stim_data_trial_averged = take_average(stim_data, number_of_average)
% identical to the compute script's helper
[~, neuron] = size(stim_data);
stim_data_trial_averged = zeros(50, neuron);
for i = 1:50
    temp = stim_data(i*number_of_average-(number_of_average-1):i*number_of_average, :);
    stim_data_trial_averged(i,:) = mean(temp, 1);
end
end
