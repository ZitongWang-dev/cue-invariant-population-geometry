%{
% FILENAME: visualize_procrustes_alignment_pca.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   This script visualizes the geometric alignment of neural population
%   codes for different stimulus types (AC, EC, EX). It loads trial-averaged
%   firing rates, applies a rotation-only Procrustes transformation to align
%   the 'AC' and 'EC' neural response geometries to the 'EX' geometry, and
%   then uses PCA to find a shared low-dimensional space.
%
%   The final output is a 2x3 subplot figure showing the population data in
%   the first two principal components. The top row displays the original,
%   unaligned data, while the bottom row shows the data after alignment.
%   Stimulus image thumbnails are overlaid on the scatter plot for clarity.
%}

%% House-keeping
clc; clear; close all;

%% USER-MODIFIABLE PARAMETERS ----------------------------------------------
monkey      = 'KO';             % Monkey ID
vp          = 'V1';             % Visual pathway / area

timewindow  = [330 630];        % Time window for trial-averaging (ms)
points_idx  = 1:50;             % Indices of trials/images to plot
circle_size = 525;              % Marker size for highlighted points
patch_width = 1.5;              % Desired width of each thumbnail (in data-units)

%% PATHS --------------------------------------------------------------------
% Directory containing *.mat spiking files and image thumbnails
img_path     = fullfile(pwd, 'Stimuli_concepts_pngs');

% Build wildcard pattern (used with DIR—not LOAD!)
mat_pattern  = fullfile(pwd, monkey, vp, sprintf('%s*%s*allstim.mat', monkey, vp));

%% LOAD SPIKING DATA --------------------------------------------------------
files = dir(mat_pattern);
assert(~isempty(files), ['No matching .mat file found for pattern:\n' mat_pattern]);
matfile_path      = fullfile(files(1).folder, files(1).name);

% Expecting variable `three_stim_array` inside the MAT file
spike_data_struct = load(matfile_path, 'three_stim_array');
spike_data        = spike_data_struct.three_stim_array;

%% COMBINE TRIALS: compute mean firing per trial ---------------------------
combined_array = cell(size(spike_data));
for i = 1:numel(spike_data)
    spike_onestim     = spike_data{i};           
    means_per_neuron  = cellfun(@mean, spike_onestim, 'UniformOutput', false);
    combined_array{i} = means_per_neuron;
end

%% Load the spikes for the time window and normalize it.
data_trial_averaged = cell(size(combined_array));

for i = 1:numel(combined_array)
    [tmp, ~]               = spike_loader(combined_array{i}, timewindow, 1);
    tmp                    = tmp';                     
    tmp                    = zscore(tmp);              
    data_trial_averaged{i} = tmp;
end
%% PCA of 3 STIM TYPES + 2 TRANSFORMED TYPES -----------------------
ac = data_trial_averaged{1};  
ec = data_trial_averaged{2};  
ex = data_trial_averaged{3};  

% Procrustes alignment (AC→EX, EC→EX)
[~, ac_trans_full, transform4AC] = procrustes(ex, ac);
ac_trans             = ac * transform4AC.T;

[~, ec_trans_full, transform4EC] = procrustes(ex, ec);
ec_trans             = ec * transform4EC.T;

% PCA on concatenated data
all_data = [ac; ec; ex; ac_trans; ec_trans];
[~, score, ~] = pca(all_data);

% Slice scores back into groups (first two PCs)
ac_pca        = score(1:50,       1:2);
ec_pca        = score(51:100,     1:2);
ex_pca        = score(101:150,    1:2);
ac_trans_pca  = score(151:200,    1:2);
ec_trans_pca  = score(201:250,    1:2);

%% VISUALISATION ------------------------------------------------------------
figure('Color','w');

ax(1) = subplot(2,3,1);
    plot_images(ex_pca,       'ex', 's', patch_width, img_path, points_idx, circle_size);
    title('EX'); xlabel('PC1'); ylabel('PC2');
    xlim([-15 15]); ylim([-8 12]); colormap gray;

ax(2) = subplot(2,3,2);
    plot_images(ec_pca,       'ec', 'l', patch_width, img_path, points_idx, circle_size);
    title('EC'); xlabel('PC1'); ylabel('PC2');
    xlim([-15 15]); ylim([-8 12]); colormap gray;

ax(3) = subplot(2,3,3);
    plot_images(ac_pca,       'ac', 's', patch_width, img_path, points_idx, circle_size);
    title('AC'); xlabel('PC1'); ylabel('PC2');
    xlim([-15 15]); ylim([-8 12]); colormap gray;

ax(4) = subplot(2,3,4);
    plot_images(ex_pca,       'ex', 's', patch_width, img_path, points_idx, circle_size);
    title('EX'); xlabel('PC1'); ylabel('PC2');
    xlim([-15 20]); ylim([-8 12]); colormap gray;

ax(5) = subplot(2,3,5);
    plot_images(ec_trans_pca, 'ec', 'l', patch_width, img_path, points_idx, circle_size);
    title('Transformed EC'); xlabel('PC1'); ylabel('PC2');
    xlim([-15 20]); ylim([-8 12]); colormap gray;

ax(6) = subplot(2,3,6);
    plot_images(ac_trans_pca, 'ac', 's', patch_width, img_path, points_idx, circle_size);
    title('Transformed AC'); xlabel('PC1'); ylabel('PC2');
    xlim([-15 20]); ylim([-8 12]); colormap gray;

% Synchronise axis limits across all subplots
linkaxes(ax, 'xy');

sgtitle(sprintf('%s%s extrema of AC', monkey, vp));

%% -------------------------------------------------------------------------
function plot_images(coords, stim_type, img_size_code, patch_width, img_path, points_idx, circle_size)
% plot_images  Overlay image thumbnails at specified 2-D coordinates.
%   patch_width : desired width of each thumbnail in data-units (PC1/PC2 space)

    hold on;

    % Plot each thumbnail at its coordinate, preserving aspect ratio
    for ii = 1:numel(points_idx)
        idx      = points_idx(ii);
        fname    = sprintf('%s_%03d_%s0.png', stim_type, idx, img_size_code);
        im       = double(imread(fullfile(img_path, fname)));
        im = flipud(im);
        [h, w]   = size(im);

        % absolute width = patch_width; height scaled by aspect ratio
        scaleX = patch_width;
        scaleY = patch_width * (h / w);

        xx = linspace(-scaleX/2, scaleX/2, h) + coords(idx,1);
        yy = linspace(-scaleY/2, scaleY/2, w) + coords(idx,2);

        image('XData', xx, 'YData', yy, 'CData', im, 'CDataMapping', 'direct');
    end

    % Highlight selected markers
    cmap       = colormap('jet');
    marker_idx = [13, 15, 41, 46];
    cidx       = ceil(linspace(1, size(cmap,1), numel(marker_idx)));
    scatter(coords(marker_idx,1), coords(marker_idx,2), ...
            circle_size, cmap(cidx,:), 's', 'LineWidth', 2);

    % Connect highlighted markers
    path_idx = [15, 41, 46, 13, 15];
    plot(coords(path_idx,1), coords(path_idx,2), 'LineWidth', 1.5, ...
         'Color', [0 0.4470 0.7410]);

    hold off;
end
