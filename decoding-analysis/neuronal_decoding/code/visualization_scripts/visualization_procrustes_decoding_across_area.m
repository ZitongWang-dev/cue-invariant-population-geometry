%{
Filename: visualization_procrustes_decoding_across_area.m
Author: Zitong Wang
Date: 2025-07-10

Description:
    Visualizes feed‑forward (V1→V2) and feedback (V2→V1) Procrustes decoding
    across cortical areas. Loads all nine stimulus‑pair result files, groups
    them by training cue (v1ac, v1ex, v1ec / v2ac, v2ex, v2ec), computes
    mean ± SEM for selected accuracy columns, and plots paired bar charts with
    error bars.

%}
%%
clc
clear
% Directory path where the files are stored
monkey = 'FR';
feedforward_folder = ['..\..\results\decoding_outputs\','procrustes_decoding_across_area_results','\', monkey,'\','V1 decoding V2','\']
feedback_folder = ['..\..\results\decoding_outputs\','procrustes_decoding_across_area_results','\', monkey,'\','V2 decoding V1','\']
feedforward_results = readDataFolder(feedforward_folder);
feedback_results = readDataFolder(feedback_folder);
%% group results with same training cue
% Assuming feedforward_results and feedback_results are loaded

% Initialize variables to store the combined data for feedforward and feedback
grouped_feedforward = struct();
grouped_feedback = struct();
numPairs = 9;
% Get the field names from feedforward_results
ff_field = fieldnames(feedforward_results);
fb_field = fieldnames(feedback_results);

% Group and vertically combine arrays in feedforward_results
for i = 1:numPairs
    field_name = ff_field{i};
    
    % Extract the second name (e.g., 'v1ac' from 'v2acv1ac_results')
    second_name = field_name(5:8);  
    
    % Check if this second name already exists in the grouped structure
    if isfield(grouped_feedforward, second_name)
        % If it exists, concatenate the new data
        grouped_feedforward.(second_name) = [grouped_feedforward.(second_name); feedforward_results.(field_name)];
    else
        % If it doesn't exist, initialize with the current data
        grouped_feedforward.(second_name) = feedforward_results.(field_name);
    end
end

% Repeat the same process for feedback_results
for i = 1:numPairs
    field_name = fb_field{i};
    
    % Extract the second name (e.g., 'v1ac' from 'v2acv1ac_results')
    second_name = field_name(5:8);  
    
    % Check if this second name already exists in the grouped structure
    if isfield(grouped_feedback, second_name)
        % If it exists, concatenate the new data
        grouped_feedback.(second_name) = [grouped_feedback.(second_name); feedback_results.(field_name)];
    else
        % If it doesn't exist, initialize with the current data
        grouped_feedback.(second_name) = feedback_results.(field_name);
    end
end


%%
% Assuming the grouped_feedforward and grouped_feedback are already created
selectedCols = [1, 2, 6,4,8];  % Only using columns 1, 2, and 6
numCols = length(selectedCols);

% Desired group order: AC, EX, EC
v1GroupOrder = {'v1ac', 'v1ex', 'v1ec'};
v2GroupOrder = {'v2ac', 'v2ex', 'v2ec'};
numGroups = length(v1GroupOrder);

feedforward_means = zeros(numGroups, numCols);
feedback_means = zeros(numGroups, numCols);
feedforward_sems = zeros(numGroups, numCols);
feedback_sems = zeros(numGroups, numCols);

% Compute row-wise mean and SEM for each group in the desired order
for i = 1:numGroups
    % Feedforward data
    ff_data = grouped_feedforward.(v1GroupOrder{i})(:, selectedCols);
    feedforward_means(i, :) = mean(ff_data, 1);  % Row-wise mean
    feedforward_sems(i, :) = std(ff_data, 0, 1) / sqrt(size(ff_data, 1));  % SEM

    % Feedback data
    fb_data = grouped_feedback.(v2GroupOrder{i})(:, selectedCols);
    feedback_means(i, :) = mean(fb_data, 1);  % Row-wise mean
    feedback_sems(i, :) = std(fb_data, 0, 1) / sqrt(size(fb_data, 1));  % SEM
end

% Combine the data into alternating feedforward and feedback pairs
combined_means = zeros(numGroups, 2 * numCols);
combined_sems = zeros(numGroups, 2 * numCols);

for i = 1:numCols
    combined_means(:, 2*i-1) = feedforward_means(:, i);  % Feedforward
    combined_means(:, 2*i) = feedback_means(:, i);  % Feedback
    
    combined_sems(:, 2*i-1) = feedforward_sems(:, i);  % Feedforward SEM
    combined_sems(:, 2*i) = feedback_sems(:, i);  % Feedback SEM
end

% Now, create the bar graph
width = 1080;
height = 540;
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
    errorbar(x, combined_means(:, i), combined_sems(:, i), 'r', 'linestyle', 'none','HandleVisibility','off',LineWidth=1);
end
b(3).LineStyle = '--'; b(4).LineStyle = '--';
b(7).LineStyle = ':';  b(8).LineStyle = ':';
b(9).LineStyle = '-.';  b(10).LineStyle = '-.';
for i = 3:10, b(i).LineWidth = 1.5; end
% Custom xtick labels with the format 'Train [condition] Test [condition]'
xtick_labels = {'AC → AC|EX|EC', 'EX → AC|EX|EC', 'EC → AC|EX|EC'};
            
set(gca, 'xtick', 1:numGroups, 'xticklabel', xtick_labels);
% Add a horizontal dashed line at y = 0.02
yline(0.02, '--k', 'LineWidth', 1);
% Customize appearance
ylabel('Decoding accuracy');
title([monkey, ' Feedforward and Feedback Decoding Accuracy']);
legend({'V1 self','V2 self','V1-to-V2 non‑PT','V2-to-V1 non‑PT','V1-to-V2 rot‑only','V2-to-V1 rot‑only','V1-to-V2 PT‑ctrl','V2-to-V1 PT‑ctrl','V1-to-V2 non‑PT‑ctrl','V2-to-V1 non‑PT‑ctrl','2 % chance'},'Location','best');
ylim([0 0.7])


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