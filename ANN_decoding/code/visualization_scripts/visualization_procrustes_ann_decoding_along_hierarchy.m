%{
% FILENAME: visualization_procrustes_ann_decoding_along_hierarchy.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Visualizes the hierarchical (cross-layer) Procrustes decoding results
%   for the AlexNet model, specifically comparing Pool 2 and Pool 5.
%   The script loads pre-computed results and generates a grouped bar plot
%   that directly compares the "feedforward" (e.g., train on Pool 2, test
%   on Pool 5) and "feedback" (train on Pool 5, test on Pool 2) decoding
%   directions for various cue combinations.
%}
%%
clc
clear
% Directory path where the files are stored
model = 'alexnet';
feedforward_folder = ['..\..\results\decoding_outputs\','procrustes_ann_decoding_along_hierarchy_results\',model,'\','alexnet POOL 2&5\','V1 decoding V2'];
feedback_folder = ['..\..\results\decoding_outputs\','procrustes_ann_decoding_along_hierarchy_results\',model,'\','alexnet POOL 2&5\','V2 decoding V1'];
feedforward_results = readDataFolder(feedforward_folder);
feedback_results = readDataFolder(feedback_folder);

%%
% Assuming feedforward_results and feedback_results are loaded

% Initialize variables to store mean and SEM values for feedforward and feedback
selectedCols = [1, 2, 6, 4, 8];  % Columns to use (in order)
numCols = length(selectedCols);
numPairs = 9;

feedforward_means = zeros(numPairs, numCols);
feedback_means = zeros(numPairs, numCols);
feedforward_sems = zeros(numPairs, numCols);
feedback_sems = zeros(numPairs, numCols);

% Get the field names dynamically from feedforward_results
ff_field = fieldnames(feedforward_results);
fb_field = fieldnames(feedback_results);

% Compute row-wise mean and SEM for each pair
for i = 1:numPairs
    % Feedforward data
    ff_data = feedforward_results.(ff_field{i})(:, selectedCols);
    feedforward_means(i, :) = mean(ff_data, 1);  % Row-wise mean
    feedforward_sems(i, :) = std(ff_data, 0, 1) / sqrt(size(ff_data, 1));  % SEM

    % Feedback data
    fb_data = feedback_results.(fb_field{i})(:, selectedCols);
    feedback_means(i, :) = mean(fb_data, 1);  % Row-wise mean
    feedback_sems(i, :) = std(fb_data, 0, 1) / sqrt(size(fb_data, 1));  % SEM
end

% Combine the data into alternating feedforward and feedback pairs
combined_means = zeros(numPairs, 2 * numCols);
combined_sems = zeros(numPairs, 2 * numCols);

for i = 1:numCols
    combined_means(:, 2*i-1) = feedforward_means(:, i);  % Feedforward
    combined_means(:, 2*i) = feedback_means(:, i);  % Feedback
    
    combined_sems(:, 2*i-1) = feedforward_sems(:, i);  % Feedforward SEM
    combined_sems(:, 2*i) = feedback_sems(:, i);  % Feedback SEM
end

% Create the bar graph
width = 2400;
height = 400;
figure('Position', [100, 100, width, height]);
hold on;

% Plot the bars for the combined feedforward and feedback means
b = bar(combined_means, 'grouped');

% Customize the colors
for i = 1:2:numCols*2
    b(i).FaceColor = [0 0.4470 0.7410];  % Blue for feedforward
    b(i+1).FaceColor = [0.8500 0.3250 0.0980];  % Orange for feedback
end

% Adding error bars for each pair
numBars = size(combined_means, 2);  % Total number of bars
groupwidth = min(0.8, numBars/(numBars + 1.5));

for i = 1:numBars
    x = b(i).XEndPoints;
    errorbar(x, combined_means(:, i), combined_sems(:, i), 'r','linestyle','none','HandleVisibility','off','LineWidth',1);
end

% Custom xtick labels with the format 'Train [condition] Test [condition]'
xtick_labels = cell(numPairs, 1);

for i = 1:numPairs
    field_name = ff_field{i};
    
    % Extract the relevant parts: train and test conditions
    train_condition = field_name(7:8);  % Last two characters are train condition
    test_condition = field_name(3:4);  % Characters 4 and 5 are the test condition
    
    % Create the label in the desired format
    xtick_labels{i} = ['Train ', upper(train_condition), ' Test ', upper(test_condition)];
end
            
set(gca, 'xtick', 1:numPairs, 'xticklabel', xtick_labels);
b(3).LineStyle = "--";
b(4).LineStyle = "--";
b(7).LineStyle = ":";
b(8).LineStyle = ":";
b(9).LineStyle = '-.';
b(10).LineStyle = '-.';

lwidth = 1.5;
for i = 3:10 b(i).LineWidth = lwidth; end
% Customize appearance
ylim([0 1.1])
ylabel('Decoding accuracy');
legend('V1 self-decoding','V2 self-decoding', ...
    'V1 non-PT-decoding V2','V2 non-PT-decoding V1', ...
    'V1 PT-decoding(only rotation) V2','V2 PT-decoding(only rotation) V1', ...
    'V1 PT control','V2 PT control', ...
    'V1 non-PT control','V2 non-PT control','Location','best')
title('Feedforward and Feedback Decoding (ALEXNET POOL 2 & 5)');

% Rotate xtick labels for clarity
% xtickangle(45);

% Add grid if needed

hold off;

%%

function results = readDataFolder(folder)

% Get a list of all .mat files in the folder
fileList = dir(fullfile(folder, '*.mat'));

% Initialize a structure to store the data from each file
results = struct();

% Loop through the files and load them
for i = 1:length(fileList)
    fileName = fileList(i).name;  % Get the file name
    filePath = fullfile(folder, fileName);  % Get the full path
    data = load(filePath);  % Load the .mat file
    fieldName = fieldnames(data);  
    % Extract the value inside the 1x1 cell
    cellData = data.(fieldName{1});
    value = cellData{1};
    % Store in a structure using the file name (without extension) as the key
    [~, name, ~] = fileparts(fileName);
    results.(name) = value;
end
end