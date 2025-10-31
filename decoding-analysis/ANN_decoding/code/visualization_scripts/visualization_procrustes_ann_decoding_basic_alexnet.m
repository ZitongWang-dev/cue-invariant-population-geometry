%{
% FILENAME: visualization_procrustes_ann_decoding_basic_alexnet.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Visualizes the Procrustes decoding results for the AlexNet model. This
%   script loads the pre-computed decoding accuracies from multiple pooling
%   layers (1, 2, and 5), processes the data by combining reciprocal pairs,
%   and generates a grouped bar plot comparing the performance across layers
%   and stimulus conditions.
%}
%%
clc
clear

%% CONFIGURATION
% Define the network to be analyzed.
network_name = 'alexnet';
%% LOAD DATA
% Load all result files for the specified layers and stimulus pairs using a
% helper function. The variable names V1, V2, V5 correspond to the results
% from pooling layers 1, 2, and 5 respectively.
V1_acec_result = result_loader(network_name, '1', 'acec');
V2_acec_result = result_loader(network_name, '2', 'acec');
V5_acec_result = result_loader(network_name, '5', 'acec');

V1_ecac_result = result_loader(network_name, '1', 'ecac');
V2_ecac_result = result_loader(network_name, '2', 'ecac');
V5_ecac_result = result_loader(network_name, '5', 'ecac');

V1_acex_result = result_loader(network_name, '1', 'acex');
V2_acex_result = result_loader(network_name, '2', 'acex');
V5_acex_result = result_loader(network_name, '5', 'acex');

V1_exac_result = result_loader(network_name, '1', 'exac');
V2_exac_result = result_loader(network_name, '2', 'exac');
V5_exac_result = result_loader(network_name, '5', 'exac');

V1_ecex_result = result_loader(network_name, '1', 'ecex');
V2_ecex_result = result_loader(network_name, '2', 'ecex');
V5_ecex_result = result_loader(network_name, '5', 'ecex');

V1_exec_result = result_loader(network_name, '1', 'exec');
V2_exec_result = result_loader(network_name, '2', 'exec');
V5_exec_result = result_loader(network_name, '5', 'exec');

%% PROCESS DATA
% Combine results from reciprocal stimulus pairs (e.g., ac-ec and ec-ac)
% to get an averaged performance measure for each pair type.
V1_acec_combined = [V1_acec_result; V1_ecac_result];
V2_acec_combined = [V2_acec_result; V2_ecac_result];
V5_acec_combined = [V5_acec_result; V5_ecac_result];

V1_acex_combined = [V1_acex_result; V1_exac_result];
V2_acex_combined = [V2_acex_result; V2_exac_result];
V5_acex_combined = [V5_acex_result; V5_exac_result];

V1_ecex_combined = [V1_ecex_result; V1_exec_result];
V2_ecex_combined = [V2_ecex_result; V2_exec_result];
V5_ecex_combined = [V5_ecex_result; V5_exec_result];
%% CALCULATE STATISTICS AND ASSEMBLE PLOTTING MATRICES
% Calculate mean and Standard Error of the Mean (SEM) for each condition.
total_repeat = length(V1_acec_combined);
calc_stats = @(data) [mean(data, 1); std(data, 0, 1) / sqrt(total_repeat)];

V1_acec_meanstd = calc_stats(V1_acec_combined);
V2_acec_meanstd = calc_stats(V2_acec_combined);
V5_acec_meanstd = calc_stats(V5_acec_combined);

V1_acex_meanstd = calc_stats(V1_acex_combined);
V2_acex_meanstd = calc_stats(V2_acex_combined);
V5_acex_meanstd = calc_stats(V5_acex_combined);

V1_ecex_meanstd = calc_stats(V1_ecex_combined);
V2_ecex_meanstd = calc_stats(V2_ecex_combined);
V5_ecex_meanstd = calc_stats(V5_ecex_combined);

% Define the column index for the desired partial transformation to plot.
% Result columns: 1:self, 2:non-PT, 3:full-PT, 4:PT-control,
% 5:scale-only, 6:rotation-only, 7:translation-only, 8:non-PT-control.
partial_trans_idx = 6; % Use rotation-only for the main PT bar

% Assemble the accuracy data into a matrix for the bar plot.
% Each row is a stimulus pair; columns group by decoding type and layer.
acc_group = [V1_acec_meanstd(1,1) V2_acec_meanstd(1,1) V5_acec_meanstd(1,1) V1_acec_meanstd(1,2) V2_acec_meanstd(1,2) V5_acec_meanstd(1,2) V1_acec_meanstd(1,partial_trans_idx) V2_acec_meanstd(1,partial_trans_idx) V5_acec_meanstd(1,partial_trans_idx) V1_acec_meanstd(1,4) V2_acec_meanstd(1,4) V5_acec_meanstd(1,4) V1_acec_meanstd(1,8) V2_acec_meanstd(1,8) V5_acec_meanstd(1,8);
             V1_ecex_meanstd(1,1) V2_ecex_meanstd(1,1) V5_ecex_meanstd(1,1) V1_ecex_meanstd(1,2) V2_ecex_meanstd(1,2) V5_ecex_meanstd(1,2) V1_ecex_meanstd(1,partial_trans_idx) V2_ecex_meanstd(1,partial_trans_idx) V5_ecex_meanstd(1,partial_trans_idx) V1_ecex_meanstd(1,4) V2_ecex_meanstd(1,4) V5_ecex_meanstd(1,4) V1_ecex_meanstd(1,8) V2_ecex_meanstd(1,8) V5_ecex_meanstd(1,8);
             V1_acex_meanstd(1,1) V2_acex_meanstd(1,1) V5_acex_meanstd(1,1) V1_acex_meanstd(1,2) V2_acex_meanstd(1,2) V5_acex_meanstd(1,2) V1_acex_meanstd(1,partial_trans_idx) V2_acex_meanstd(1,partial_trans_idx) V5_acex_meanstd(1,partial_trans_idx) V1_acex_meanstd(1,4) V2_acex_meanstd(1,4) V5_acex_meanstd(1,4) V1_acex_meanstd(1,8) V2_acex_meanstd(1,8) V5_acex_meanstd(1,8)];

% Assemble the corresponding error (SEM) data.
err_group = [V1_acec_meanstd(2,1) V2_acec_meanstd(2,1) V5_acec_meanstd(2,1) V1_acec_meanstd(2,2) V2_acec_meanstd(2,2) V5_acec_meanstd(2,2) V1_acec_meanstd(2,partial_trans_idx) V2_acec_meanstd(2,partial_trans_idx) V5_acec_meanstd(2,partial_trans_idx) V1_acec_meanstd(2,4) V2_acec_meanstd(2,4) V5_acec_meanstd(2,4) V1_acec_meanstd(2,8) V2_acec_meanstd(2,8) V5_acec_meanstd(2,8);
             V1_ecex_meanstd(2,1) V2_ecex_meanstd(2,1) V5_ecex_meanstd(2,1) V1_ecex_meanstd(2,2) V2_ecex_meanstd(2,2) V5_ecex_meanstd(2,2) V1_ecex_meanstd(2,partial_trans_idx) V2_ecex_meanstd(2,partial_trans_idx) V5_ecex_meanstd(2,partial_trans_idx) V1_ecex_meanstd(2,4) V2_ecex_meanstd(2,4) V5_ecex_meanstd(2,4) V1_ecex_meanstd(2,8) V2_ecex_meanstd(2,8) V5_ecex_meanstd(2,8);
             V1_acex_meanstd(2,1) V2_acex_meanstd(2,1) V5_acex_meanstd(2,1) V1_acex_meanstd(2,2) V2_acex_meanstd(2,2) V5_acex_meanstd(2,2) V1_acex_meanstd(2,partial_trans_idx) V2_acex_meanstd(2,partial_trans_idx) V5_acex_meanstd(2,partial_trans_idx) V1_acex_meanstd(2,4) V2_acex_meanstd(2,4) V5_acex_meanstd(2,4) V1_acex_meanstd(2,8) V2_acex_meanstd(2,8) V5_acex_meanstd(2,8)];

%% GENERATE BAR PLOT
% Create a grouped bar plot to visualize the decoding results.
x_axis_combined = categorical(["EC-AC", "EC-EX", "AC-EX"]);
x_axis_combined = reordercats(x_axis_combined, ["EC-AC", "EC-EX", "AC-EX"]);

figure('Position', [100, 100, 2160, 540]);
b = bar(x_axis_combined, acc_group);

% --- Customize Bar Styles ---
% Set colors for each layer's bars.
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]}; % Blue, Orange, Yellow
for i = 1:size(b, 2)
    b(i).FaceColor = colors{mod(i-1, 3) + 1};
end

% Set line styles to differentiate decoding types.
lwidth = 1.5;
for i = 4:6,  b(i).LineStyle = '--'; b(i).LineWidth = lwidth; end % non-PT
for i = 7:9,  b(i).LineWidth = lwidth; end %PT
for i = 10:12, b(i).LineStyle = ':';  b(i).LineWidth = lwidth; end % PT-control
for i = 13:15, b(i).LineStyle = '-.'; b(i).LineWidth = lwidth; end % non-PT control

% --- Add Labels, Title, and Error Bars ---
ylim([0 1.1]);
yline(0.02, '--', 'Color', [0.5 0.5 0.5]); % 2% chance level line

legend('Pool 1 self', 'Pool 2 self', 'Pool 5 self', ...
    'Pool 1 non-PT', 'Pool 2 non-PT', 'Pool 5 non-PT', ...
    'Pool 1 PT (rotation-only)', 'Pool 2 PT (rotation-only)', 'Pool 5 PT (rotation-only)', ...
    'Pool 1 PT control', 'Pool 2 PT control', 'Pool 5 PT control', ...
    'Pool 1 non-PT control', 'Pool 2 non-PT control', 'Pool 5 non-PT control', ...
    'Theoretical Chance Level: 2%');

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

% Plot the error bars.
errorbar(x', acc_group, err_group, 'r', 'linestyle', 'none', 'HandleVisibility', 'off', 'LineWidth', 1);
hold off;

%%
function dist_result = result_loader(monkey,vp,stimpair)
dist_result = load(strcat('..\..\results\decoding_outputs\','procrustes_ann_decoding_basic_results\',monkey,'\pool',vp,'\',stimpair,'_','results','.mat'),[stimpair,'_','results']); 
dist_result = dist_result.([stimpair,'_','results']);

end