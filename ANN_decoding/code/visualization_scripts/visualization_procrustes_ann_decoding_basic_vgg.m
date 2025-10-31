%{
% FILENAME: visualization_procrustes_ann_decoding_basic_vgg.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Visualizes the Procrustes decoding results specifically for the VGG model.
%   The script loads pre-computed accuracies from all 5 VGG pooling layers,
%   combines data from reciprocal stimulus pairs, calculates statistics, and
%   generates a comprehensive bar plot comparing performance across all
%   layers and conditions.
%}

%% INITIALIZATION
clc;
clear;

%% CONFIGURATION
% Define the network to be analyzed.
network_name = 'vgg';

%% LOAD DATA
% Load all result files for the specified layers and stimulus pairs.
% The variables V1, V2, etc., correspond to VGG pooling layers 1, 2, etc.
V1_acec_result = result_loader(network_name,'1','acec');
V2_acec_result = result_loader(network_name,'2','acec');
V3_acec_result = result_loader(network_name,'3','acec');
V4_acec_result = result_loader(network_name,'4','acec');
V5_acec_result = result_loader(network_name,'5','acec');

V1_ecac_result = result_loader(network_name,'1','ecac');
V2_ecac_result = result_loader(network_name,'2','ecac');
V3_ecac_result = result_loader(network_name,'3','ecac');
V4_ecac_result = result_loader(network_name,'4','ecac');
V5_ecac_result = result_loader(network_name,'5','ecac');

V1_acex_result = result_loader(network_name,'1','acex');
V2_acex_result = result_loader(network_name,'2','acex');
V3_acex_result = result_loader(network_name,'3','acex');
V4_acex_result = result_loader(network_name,'4','acex');
V5_acex_result = result_loader(network_name,'5','acex');

V1_exac_result = result_loader(network_name,'1','exac');
V2_exac_result = result_loader(network_name,'2','exac');
V3_exac_result = result_loader(network_name,'3','exac');
V4_exac_result = result_loader(network_name,'4','exac');
V5_exac_result = result_loader(network_name,'5','exac');

V1_ecex_result = result_loader(network_name,'1','ecex');
V2_ecex_result = result_loader(network_name,'2','ecex');
V3_ecex_result = result_loader(network_name,'3','ecex');
V4_ecex_result = result_loader(network_name,'4','ecex');
V5_ecex_result = result_loader(network_name,'5','ecex');

V1_exec_result = result_loader(network_name,'1','exec');
V2_exec_result = result_loader(network_name,'2','exec');
V3_exec_result = result_loader(network_name,'3','exec');
V4_exec_result = result_loader(network_name,'4','exec');
V5_exec_result = result_loader(network_name,'5','exec');

%% PROCESS DATA
% Combine results from reciprocal stimulus pairs (e.g., ac-ec and ec-ac)
% to get an averaged performance measure for each pair type.
V1_acec_combined = [V1_acec_result ; V1_ecac_result];
V2_acec_combined = [V2_acec_result ; V2_ecac_result];
V3_acec_combined = [V3_acec_result ; V3_ecac_result];
V4_acec_combined = [V4_acec_result ; V4_ecac_result];
V5_acec_combined = [V5_acec_result ; V5_ecac_result];

V1_acex_combined = [V1_acex_result ; V1_exac_result];
V2_acex_combined = [V2_acex_result ; V2_exac_result];
V3_acex_combined = [V3_acex_result ; V3_exac_result];
V4_acex_combined = [V4_acex_result ; V4_exac_result];
V5_acex_combined = [V5_acex_result ; V5_exac_result];

V1_ecex_combined = [V1_ecex_result ; V1_exec_result];
V2_ecex_combined = [V2_ecex_result ; V2_exec_result];
V3_ecex_combined = [V3_ecex_result ; V3_exec_result];
V4_ecex_combined = [V4_ecex_result ; V4_exec_result];
V5_ecex_combined = [V5_ecex_result ; V5_exec_result];

%% CALCULATE STATISTICS AND ASSEMBLE PLOTTING MATRICES
% Calculate mean and Standard Error of the Mean (SEM) for each condition.
total_repeat = length(V1_acec_combined);
calc_stats = @(data) [mean(data, 1); std(data, 0, 1) / sqrt(total_repeat)];

V1_acec_meanstd = calc_stats(V1_acec_combined);
V2_acec_meanstd = calc_stats(V2_acec_combined);
V3_acec_meanstd = calc_stats(V3_acec_combined);
V4_acec_meanstd = calc_stats(V4_acec_combined);
V5_acec_meanstd = calc_stats(V5_acec_combined);

V1_acex_meanstd = calc_stats(V1_acex_combined);
V2_acex_meanstd = calc_stats(V2_acex_combined);
V3_acex_meanstd = calc_stats(V3_acex_combined);
V4_acex_meanstd = calc_stats(V4_acex_combined);
V5_acex_meanstd = calc_stats(V5_acex_combined);

V1_ecex_meanstd = calc_stats(V1_ecex_combined);
V2_ecex_meanstd = calc_stats(V2_ecex_combined);
V3_ecex_meanstd = calc_stats(V3_ecex_combined);
V4_ecex_meanstd = calc_stats(V4_ecex_combined);
V5_ecex_meanstd = calc_stats(V5_ecex_combined);

% Define the column index for the desired partial transformation to plot.
% Result columns: 1:self, 2:non-PT, 3:full-PT, 4:PT-control,
% 5:scale-only, 6:rotation-only, 7:translation-only, 8:non-PT-control.
partial_trans_idx = 6; % Use rotation-only for the main PT bar

% Assemble the accuracy data into a matrix for the bar plot.
acc_group = [V1_acec_meanstd(1,1) V2_acec_meanstd(1,1) V3_acec_meanstd(1,1) V4_acec_meanstd(1,1) V5_acec_meanstd(1,1) V1_acec_meanstd(1,2) V2_acec_meanstd(1,2) V3_acec_meanstd(1,2) V4_acec_meanstd(1,2) V5_acec_meanstd(1,2) V1_acec_meanstd(1,partial_trans_idx) V2_acec_meanstd(1,partial_trans_idx) V3_acec_meanstd(1,partial_trans_idx) V4_acec_meanstd(1,partial_trans_idx) V5_acec_meanstd(1,partial_trans_idx) V1_acec_meanstd(1,4) V2_acec_meanstd(1,4) V3_acec_meanstd(1,4) V4_acec_meanstd(1,4) V5_acec_meanstd(1,4) V1_acec_meanstd(1,8) V2_acec_meanstd(1,8) V3_acec_meanstd(1,8) V4_acec_meanstd(1,8) V5_acec_meanstd(1,8);
             V1_ecex_meanstd(1,1) V2_ecex_meanstd(1,1) V3_ecex_meanstd(1,1) V4_ecex_meanstd(1,1) V5_ecex_meanstd(1,1) V1_ecex_meanstd(1,2) V2_ecex_meanstd(1,2) V3_ecex_meanstd(1,2) V4_ecex_meanstd(1,2) V5_ecex_meanstd(1,2) V1_ecex_meanstd(1,partial_trans_idx) V2_ecex_meanstd(1,partial_trans_idx) V3_ecex_meanstd(1,partial_trans_idx) V4_ecex_meanstd(1,partial_trans_idx) V5_ecex_meanstd(1,partial_trans_idx) V1_ecex_meanstd(1,4) V2_ecex_meanstd(1,4) V3_ecex_meanstd(1,4) V4_ecex_meanstd(1,4) V5_ecex_meanstd(1,4) V1_ecex_meanstd(1,8) V2_ecex_meanstd(1,8) V3_ecex_meanstd(1,8) V4_ecex_meanstd(1,8) V5_ecex_meanstd(1,8);
             V1_acex_meanstd(1,1) V2_acex_meanstd(1,1) V3_acex_meanstd(1,1) V4_acex_meanstd(1,1) V5_acex_meanstd(1,1) V1_acex_meanstd(1,2) V2_acex_meanstd(1,2) V3_acex_meanstd(1,2) V4_acex_meanstd(1,2) V5_acex_meanstd(1,2) V1_acex_meanstd(1,partial_trans_idx) V2_acex_meanstd(1,partial_trans_idx) V3_acex_meanstd(1,partial_trans_idx) V4_acex_meanstd(1,partial_trans_idx) V5_acex_meanstd(1,partial_trans_idx) V1_acex_meanstd(1,4) V2_acex_meanstd(1,4) V3_acex_meanstd(1,4) V4_acex_meanstd(1,4) V5_acex_meanstd(1,4) V1_acex_meanstd(1,8) V2_acex_meanstd(1,8) V3_acex_meanstd(1,8) V4_acex_meanstd(1,8) V5_acex_meanstd(1,8)];

% Assemble the corresponding error (SEM) data.
err_group = [V1_acec_meanstd(2,1) V2_acec_meanstd(2,1) V3_acec_meanstd(2,1) V4_acec_meanstd(2,1) V5_acec_meanstd(2,1) V1_acec_meanstd(2,2) V2_acec_meanstd(2,2) V3_acec_meanstd(2,2) V4_acec_meanstd(2,2) V5_acec_meanstd(2,2) V1_acec_meanstd(2,partial_trans_idx) V2_acec_meanstd(2,partial_trans_idx) V3_acec_meanstd(2,partial_trans_idx) V4_acec_meanstd(2,partial_trans_idx) V5_acec_meanstd(2,partial_trans_idx) V1_acec_meanstd(2,4) V2_acec_meanstd(2,4) V3_acec_meanstd(2,4) V4_acec_meanstd(2,4) V5_acec_meanstd(2,4) V1_acec_meanstd(2,8) V2_acec_meanstd(2,8) V3_acec_meanstd(2,8) V4_acec_meanstd(2,8) V5_acec_meanstd(2,8);
             V1_ecex_meanstd(2,1) V2_ecex_meanstd(2,1) V3_ecex_meanstd(2,1) V4_ecex_meanstd(2,1) V5_ecex_meanstd(2,1) V1_ecex_meanstd(2,2) V2_ecex_meanstd(2,2) V3_ecex_meanstd(2,2) V4_ecex_meanstd(2,2) V5_ecex_meanstd(2,2) V1_ecex_meanstd(2,partial_trans_idx) V2_ecex_meanstd(2,partial_trans_idx) V3_ecex_meanstd(2,partial_trans_idx) V4_ecex_meanstd(2,partial_trans_idx) V5_ecex_meanstd(2,partial_trans_idx) V1_ecex_meanstd(2,4) V2_ecex_meanstd(2,4) V3_ecex_meanstd(2,4) V4_ecex_meanstd(2,4) V5_ecex_meanstd(2,4) V1_ecex_meanstd(2,8) V2_ecex_meanstd(2,8) V3_ecex_meanstd(2,8) V4_ecex_meanstd(2,8) V5_ecex_meanstd(2,8);
             V1_acex_meanstd(2,1) V2_acex_meanstd(2,1) V3_acex_meanstd(2,1) V4_acex_meanstd(2,1) V5_acex_meanstd(2,1) V1_acex_meanstd(2,2) V2_acex_meanstd(2,2) V3_acex_meanstd(2,2) V4_acex_meanstd(2,2) V5_acex_meanstd(2,2) V1_acex_meanstd(2,partial_trans_idx) V2_acex_meanstd(2,partial_trans_idx) V3_acex_meanstd(2,partial_trans_idx) V4_acex_meanstd(2,partial_trans_idx) V5_acex_meanstd(2,partial_trans_idx) V1_acex_meanstd(2,4) V2_acex_meanstd(2,4) V3_acex_meanstd(2,4) V4_acex_meanstd(2,4) V5_acex_meanstd(2,4) V1_acex_meanstd(2,8) V2_acex_meanstd(2,8) V3_acex_meanstd(2,8) V4_acex_meanstd(2,8) V5_acex_meanstd(2,8)];

%% GENERATE BAR PLOT
% Create a grouped bar plot to visualize the decoding results.
x_axis_combined = categorical(["EC-AC", "EC-EX", "AC-EX"]);
x_axis_combined = reordercats(x_axis_combined, ["EC-AC", "EC-EX", "AC-EX"]);

figure('Position', [0, 0, 2700, 540]);
b = bar(x_axis_combined, acc_group);

% --- Customize Bar Styles ---
% Set colors for each of the 5 layers' bars.
for i = 1:5
    b(i*5-4).FaceColor = [0 0.4470 0.7410];      % Pool 1 Color
    b(i*5-3).FaceColor = [0.8500 0.3250 0.0980]; % Pool 2 Color
    b(i*5-2).FaceColor = [0.9290 0.6940 0.1250]; % Pool 3 Color
    b(i*5-1).FaceColor = [0.4940 0.1840 0.5560]; % Pool 4 Color
    b(i*5).FaceColor = [0.4660 0.6740 0.1880];   % Pool 5 Color
end

% Set line styles to differentiate decoding types.
b(6).LineStyle = "--";
b(7).LineStyle = "--";
b(8).LineStyle = "--";
b(9).LineStyle = "--";
b(10).LineStyle = "--";
b(16).LineStyle = ":";
b(17).LineStyle = ":";
b(18).LineStyle = ":";
b(19).LineStyle = ":";
b(20).LineStyle = ":";
b(21).LineStyle = "-.";
b(22).LineStyle = "-.";
b(23).LineStyle = "-.";
b(24).LineStyle = "-.";
b(25).LineStyle = "-.";
lwidth = 1.5;
for i = 6:25, b(i).LineWidth = lwidth; end

% --- Add Labels, Title, and Error Bars ---
ylim([0 1.1]);
yline(0.02, '--'); % 2% chance level line

legend('Pool 1 self', 'Pool 2 self', 'Pool 3 self', 'Pool 4 self', 'Pool 5 self', ...
       'Pool 1 non-PT', 'Pool 2 non-PT', 'Pool 3 non-PT', 'Pool 4 non-PT', 'Pool 5 non-PT', ...
       'Pool 1 PT (rotation-only)', 'Pool 2 PT (rotation-only)', 'Pool 3 PT (rotation-only)', 'Pool 4 PT (rotation-only)', 'Pool 5 PT (rotation-only)', ...
       'Pool 1 PT control', 'Pool 2 PT control', 'Pool 3 PT control', 'Pool 4 PT control', 'Pool 5 PT control', ...
       'Pool 1 non-PT control', 'Pool 2 non-PT control', 'Pool 3 non-PT control', 'Pool 4 non-PT control', 'Pool 5 non-PT control', ...
       'Theoretical Chance Level: 2%', 'Location', 'bestoutside');

ylabel('Decoding Accuracy');
xlabel('Stimulus Pair');
title(strcat(upper(network_name), ' Cue-transfer Decoding (PT-Rotation Only)'));

hold on;
% Calculate x-coordinates for error bars.
[ngroups, nbars] = size(acc_group);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i, :) = b(i).XEndPoints;
end

% Plot the error bars using original specified style.
errorbar(x', acc_group, err_group, 'r', 'linestyle', 'none', 'HandleVisibility', 'off', 'LineWidth', 1);
hold off;
%% HELPER FUNCTIONS
function dist_result = result_loader(network_name, layer, stimpair)
    % Loads a specific result file and extracts the relevant data structure.
    file_path = fullfile('..', '..', 'results', 'decoding_outputs', 'procrustes_ann_decoding_basic_results', network_name, ['pool', layer], [stimpair, '_results.mat']);
    dist_result = load(file_path, [stimpair, '_results']);
    dist_result = dist_result.([stimpair, '_results']);
end