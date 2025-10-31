%{
% FILENAME: visualization_procrustes_Gabor_decoding_basic.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Visualizes and compares the Procrustes decoding results from two
%   Gabor filter models: a "simple cell" model ('even and odd combined
%   downsampled' or 'even and odd combined') and a "complex cell" model. The script loads pre-computed
%   accuracies for both models, calculates summary statistics, and generates
%   a grouped bar plot to compare their performance across various stimulus
%   pairs and decoding conditions.
%}

%% INITIALIZATION
clc;
clear;

%% LOAD DATA
% Load the decoding result files for the two Gabor models.
% The 'V1' prefixed variables correspond to the simple cell model.
% The 'V2' prefixed variables correspond to the complex cell model.
V1_acec_result = result_loader('even and odd combined downsampled','acec');
V2_acec_result = result_loader('complex','acec');

V1_ecac_result = result_loader('even and odd combined downsampled','ecac');
V2_ecac_result = result_loader('complex','ecac');

V1_acex_result = result_loader('even and odd combined downsampled','acex');
V2_acex_result = result_loader('complex','acex');

V1_exac_result = result_loader('even and odd combined downsampled','exac');
V2_exac_result = result_loader('complex','exac');

V1_ecex_result = result_loader('even and odd combined downsampled','ecex');
V2_ecex_result = result_loader('complex','ecex');

V1_exec_result = result_loader('even and odd combined downsampled','exec');
V2_exec_result = result_loader('complex','exec');
%% PROCESS DATA
% Combine results from reciprocal stimulus pairs (e.g., ac-ec and ec-ac)
% to get an averaged performance measure for each pair type.
V1_acec_combined = [V1_acec_result ; V1_ecac_result];
V2_acec_combined = [V2_acec_result ; V2_ecac_result];

V1_acex_combined = [V1_acex_result ; V1_exac_result];
V2_acex_combined = [V2_acex_result ; V2_exac_result];

V1_ecex_combined = [V1_ecex_result ; V1_exec_result];
V2_ecex_combined = [V2_ecex_result ; V2_exec_result];


%% GENERATE SINGLE PANEL PLOT
% This section calculates statistics and creates a single, comprehensive
% bar plot comparing the simple and complex Gabor models.

% --- Calculate Statistics ---
% Compute the mean and Standard Error of the Mean (SEM) for each condition.
total_repeat = length(V1_acec_combined);
calc_stats = @(data) [mean(data, 1); std(data, 0, 1) / sqrt(total_repeat)];

V1_acec_meanstd = calc_stats(V1_acec_combined);
V2_acec_meanstd = calc_stats(V2_acec_combined);

V1_acex_meanstd = calc_stats(V1_acex_combined);
V2_acex_meanstd = calc_stats(V2_acex_combined);

V1_ecex_meanstd = calc_stats(V1_ecex_combined);
V2_ecex_meanstd = calc_stats(V2_ecex_combined);

% --- Assemble Plotting Matrices ---
% Define the column index for the desired partial transformation to plot.
% Result columns: 1:self, 2:non-PT, 3:full-PT, 4:PT-control,
% 5:scale-only, 6:rotation-only, 7:translation-only, 8:non-PT-control.
partial_trans_idx = 6; % Use rotation-only for the main PT bar

% Assemble the accuracy data into a matrix for the bar plot.
% Each row is a stimulus pair; columns group by decoding type and model.
acc_group = [V1_acec_meanstd(1,1) V2_acec_meanstd(1,1) V1_acec_meanstd(1,2) V2_acec_meanstd(1,2) V1_acec_meanstd(1,partial_trans_idx) V2_acec_meanstd(1,partial_trans_idx) V1_acec_meanstd(1,4) V2_acec_meanstd(1,4) V1_acec_meanstd(1,8) V2_acec_meanstd(1,8);
             V1_ecex_meanstd(1,1) V2_ecex_meanstd(1,1) V1_ecex_meanstd(1,2) V2_ecex_meanstd(1,2) V1_ecex_meanstd(1,partial_trans_idx) V2_ecex_meanstd(1,partial_trans_idx) V1_ecex_meanstd(1,4) V2_ecex_meanstd(1,4) V1_ecex_meanstd(1,8) V2_ecex_meanstd(1,8);
             V1_acex_meanstd(1,1) V2_acex_meanstd(1,1) V1_acex_meanstd(1,2) V2_acex_meanstd(1,2) V1_acex_meanstd(1,partial_trans_idx) V2_acex_meanstd(1,partial_trans_idx) V1_acex_meanstd(1,4) V2_acex_meanstd(1,4) V1_acex_meanstd(1,8) V2_acex_meanstd(1,8)];

% Assemble the corresponding error (SEM) data.
err_group = [V1_acec_meanstd(2,1) V2_acec_meanstd(2,1) V1_acec_meanstd(2,2) V2_acec_meanstd(2,2) V1_acec_meanstd(2,partial_trans_idx) V2_acec_meanstd(2,partial_trans_idx) V1_acec_meanstd(2,4) V2_acec_meanstd(2,4) V1_acec_meanstd(2,8) V2_acec_meanstd(2,8);
             V1_ecex_meanstd(2,1) V2_ecex_meanstd(2,1) V1_ecex_meanstd(2,2) V2_ecex_meanstd(2,2) V1_ecex_meanstd(2,partial_trans_idx) V2_ecex_meanstd(2,partial_trans_idx) V1_ecex_meanstd(2,4) V2_ecex_meanstd(2,4) V1_ecex_meanstd(2,8) V2_ecex_meanstd(2,8);
             V1_acex_meanstd(2,1) V2_acex_meanstd(2,1) V1_acex_meanstd(2,2) V2_acex_meanstd(2,2) V1_acex_meanstd(2,partial_trans_idx) V2_acex_meanstd(2,partial_trans_idx) V1_acex_meanstd(2,4) V2_acex_meanstd(2,4) V1_acex_meanstd(2,8) V2_acex_meanstd(2,8)];

% --- Create the Bar Plot ---
x_axis_combined = categorical(["EC-AC", "EC-EX", "AC-EX"]);
x_axis_combined = reordercats(x_axis_combined, ["EC-AC", "EC-EX", "AC-EX"]);

figure('Position', [100, 100, 1080, 540]);
b = bar(x_axis_combined, acc_group);

% Set bar colors to distinguish simple (blue) and complex (orange) models.
b(1).FaceColor = [0 0.4470 0.7410];
b(2).FaceColor = [0.8500 0.3250 0.0980];
b(3).FaceColor = [0 0.4470 0.7410];
b(4).FaceColor = [0.8500 0.3250 0.0980];
b(5).FaceColor = [0 0.4470 0.7410];
b(6).FaceColor = [0.8500 0.3250 0.0980];
b(7).FaceColor = [0 0.4470 0.7410];
b(8).FaceColor = [0.8500 0.3250 0.0980];
b(9).FaceColor = [0 0.4470 0.7410];
b(10).FaceColor = [0.8500 0.3250 0.0980];

% Set line styles to differentiate decoding types.
b(3).LineStyle = "--";
b(4).LineStyle = "--";
b(7).LineStyle = ":";
b(8).LineStyle = ":";
b(9).LineStyle = '-.';
b(10).LineStyle = '-.';
lwidth = 1.5;
for i = 3:10, b(i).LineWidth = lwidth; end

% --- Add Labels, Title, and Error Bars ---
ylim([0 1.05]);
yline(0.02, '--'); % 2% chance level line
legend('Simple Model self-decoding', 'Complex Model self-decoding', ...
       'Simple Model non-PT', 'Complex Model non-PT', ...
       'Simple Model PT (rotation-only)', 'Complex Model PT (rotation-only)', ...
       'Simple Model PT control', 'Complex Model PT control', ...
       'Simple Model non-PT control', 'Complex Model non-PT control', ...
       'Theoretical Chance Level: 2%', 'Location', 'best');
ylabel('Decoding Accuracy');
title('Gabor Response Decoding (Z-scored, Rotation-Only PT)');
xlabel('Stimulus Pair');

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
function organized_result = result_loader(vp,stimpair)
decoding_result = load(strcat('..\..','\results\decoding_outputs\',vp,'\',stimpair,'_','results','.mat'),[stimpair,'_','results']); 
decoding_result = decoding_result.([stimpair,'_','results']);
organized_result = decoding_result;
end