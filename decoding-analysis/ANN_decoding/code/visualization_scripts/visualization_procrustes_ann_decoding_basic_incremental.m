%{
% FILENAME: visualization_procrustes_ann_decoding_basic_incremental.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Visualizes the incremental Procrustes decoding results for ANNs.
%   This script plots how decoding accuracy scales with the number of units
%   in the population. It loads the incremental decoding results, aggregates
%   them across all stimulus pairs, and generates line plots with error bars
%   for the selected network (AlexNet or VGG).
%}

%% INITIALIZATION
clc;
clear;

%% CONFIGURATION
% Select the network to visualize: 'alexnet' or 'vgg'.
network_name = 'vgg';

%% PLOTTING
% This block uses a switch statement to generate the appropriate plot based
% on the selected network_name.
switch lower(network_name)
    case 'alexnet'
        %% ALEXNET VISUALIZATION
        % --- Load and Aggregate Results for Each Layer ---
        % The helper functions generate the x-axis (neuron counts) and
        % aggregate the average decoding performance across all stimulus pairs.
        layers_to_plot = {'pool1', 'pool2', 'pool5'};
        results = cell(length(layers_to_plot), 1);
        neuron_lists = cell(length(layers_to_plot), 1);
        
        for i = 1:length(layers_to_plot)
            neuron_lists{i} = gen_neuron_list(network_name, layers_to_plot{i});
            results{i} = total_results(network_name, layers_to_plot{i});
        end
        
        % --- Create the Plot ---
        figure('Position', [100, 100, 1440, 720]);
        hold on;
        
        colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
        line_widths = [1.5, 1.5, 1];

        % Plot data for each layer with error bars.
        % Solid lines: self-decoding; Dashed lines: PT-decoding.
        for i = 1:length(layers_to_plot)
            repeat_num = size(results{i}{1}, 1);
            lw = line_widths(i);
            % Self-decoding plot (solid line)
            errorbar(neuron_lists{i}, mean(results{i}{1}, 1), std(results{i}{1}, 0, 1) / sqrt(repeat_num), 'Color', colors{i}, 'LineWidth', lw);
            % PT-decoding plot (dashed line)
            errorbar(neuron_lists{i}, mean(results{i}{2}, 1), std(results{i}{2}, 0, 1) / sqrt(repeat_num), 'Color', colors{i}, 'LineStyle', '--', 'LineWidth', lw);
        end
        
        hold off;
        
        % --- Add Labels and Title ---
        legend('Pool 1 self', 'Pool 1 PT', ...
               'Pool 2 self', 'Pool 2 PT', ...
               'Pool 5 self', 'Pool 5 PT', ...
               'Location', 'southeast');
        title([upper(network_name), ' Averaged Decoding Accuracy (Rotation-only PT)']);
        xlabel('Number of Units');
        ylabel('Decoding Accuracy');
        ylim([0 1.1]);
        xlim([0 max(neuron_lists{2}) + 5]);

    case 'vgg'
        %% VGG VISUALIZATION
        % --- Load and Aggregate Results for Each Layer ---
        layers_to_plot = {'pool1', 'pool2', 'pool3', 'pool4', 'pool5'};
        results = cell(length(layers_to_plot), 1);
        neuron_lists = cell(length(layers_to_plot), 1);
        
        for i = 1:length(layers_to_plot)
            neuron_lists{i} = gen_neuron_list(network_name, layers_to_plot{i});
            results{i} = total_results(network_name, layers_to_plot{i});
        end

        % --- Create the Plot ---
        figure('Position', [100, 100, 1440, 720]);
        hold on;
        
        colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], [0.4660 0.6740 0.1880]};
        lw = 1.5;

        % Plot data for each layer with error bars.
        % Solid lines: self-decoding; Dashed lines: PT-decoding.
        for i = 1:length(layers_to_plot)
            repeat_num = size(results{i}{1}, 1);
            % Self-decoding plot (solid line)
            errorbar(neuron_lists{i}, mean(results{i}{1}, 1), std(results{i}{1}, 0, 1) / sqrt(repeat_num), 'Color', colors{i}, 'LineWidth', lw);
            % PT-decoding plot (dashed line)
            errorbar(neuron_lists{i}, mean(results{i}{2}, 1), std(results{i}{2}, 0, 1) / sqrt(repeat_num), 'Color', colors{i}, 'LineStyle', '--', 'LineWidth', lw);
        end
        
        hold off;

        % --- Add Labels and Title ---
        legend('Pool 1 self', 'Pool 1 PT', ...
               'Pool 2 self', 'Pool 2 PT', ...
               'Pool 3 self', 'Pool 3 PT', ...
               'Pool 4 self', 'Pool 4 PT', ...
               'Pool 5 self', 'Pool 5 PT', ...
               'Location', 'southeast');
        title([upper(network_name), ' Averaged Decoding Accuracy (Rotation-only PT)']);
        xlabel('Number of Units');
        ylabel('Decoding Accuracy');
        ylim([0 1.1]);
        xlim([0 max(neuron_lists{4}) + 5]);
        
    otherwise
        error("Invalid network name. Please choose 'alexnet' or 'vgg'.");
end


%% HELPER FUNCTIONS

function organized_result = result_loader(filePath, fileName)
    % Loads a single incremental result file and extracts the self-decoding
    % and rotation-only PT decoding accuracies.
    decoding_result = load(fullfile(filePath, fileName));
    decoding_result = struct2cell(decoding_result);
    decoding_result = decoding_result{1}; % Result for one stim-pair
    
    neuron_length = length(decoding_result);
    % Result columns: 1:self, 2:non-PT, ..., 6:rotation-only PT, ...
    organized_result = cell(1, 2);
    for i = 1:neuron_length
       organized_result{1} = [organized_result{1}, decoding_result{i}(:, 1)]; % Self-decoding
       organized_result{2} = [organized_result{2}, decoding_result{i}(:, 6)]; % Rotation-only PT
    end
end

function neuron_num_list = gen_neuron_list(monkey, vp)
    % Re-generates the list of neuron counts used in the incremental analysis
    % to serve as the x-axis for plotting.
    neuron_num_coode = {[96]; [250]; [249]; [64]; [128]; [256]; [511]; [507]};
    name2neuronidx = struct('alexnetpool1',1,'alexnetpool2',2,'alexnetpool5',3,'vggpool1',4,'vggpool2',5,'vggpool3',6,'vggpool4',7,'vggpool5',8);
    
    max_neuron_num = neuron_num_coode{name2neuronidx.(strcat(monkey,vp))};
    
    neuron_num_list = [];
    n = 1;
    while (5 * 2^(n-1)) < max_neuron_num
        neuron_num_list = [neuron_num_list, 5 * 2^(n-1)];
        n = n + 1;
    end
    neuron_num_list = unique([neuron_num_list, max_neuron_num]);
end

function sum_decoding_results = total_results(monkey, vp)
    % Aggregates results from all stimulus pair files for a given layer.
    % This provides an average performance for the layer.
    filePath = fullfile('..', '..', 'results', 'decoding_outputs', 'procrustes_ann_decoding_basic_incremental_results', monkey, vp);
    fileList = dir(fullfile(filePath, '*.mat'));
    fileList = {fileList.name};
    
    % Initialize cell array to store summed results for self- and PT-decoding.
    sum_decoding_results = cell(1, 2);
    for i = 1:length(fileList)
        fileName = fileList{i};
        decoding_result = result_loader(filePath, fileName);
        % Concatenate results from each file.
        for j = 1:length(sum_decoding_results)
            sum_decoding_results{j} = [sum_decoding_results{j}; decoding_result{j}];
        end
    end
end