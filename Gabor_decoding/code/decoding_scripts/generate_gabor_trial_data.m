%{
% FILENAME: generate_gabor_trial_data.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Generates synthetic trial-by-trial data from Gabor filter responses to
%   mimic neuronal recordings, preparing the data for subsequent decoding
%   analyses. The script provides distinct processing pipelines for
%   different Gabor models based on the 'gabor_type' variable.
%
%   - 'complex': Models complex cells by applying Poisson noise directly
%     to the filter energy responses.
%   - 'odd'/'even': Models the ON/OFF components of simple cells. It
%     applies half-wave rectification to the filter responses before
%     generating Poisson-distributed trials.
%   - 'even and odd combined': A final assembly step that loads the
%     previously generated 'odd' and 'even' trial data and concatenates
%     them to create a complete simple cell model dataset.
%}

%% INITIALIZATION
clc;
clear;

%% CONFIGURATION
% Select the Gabor model type to process.
% Options: 'complex', 'odd', 'even', 'even and odd combined'.
gabor_type = 'even and odd combined';

%% DATA GENERATION
% This switch block directs processing based on the selected gabor_type.
switch lower(gabor_type)
    case {'complex', 'odd', 'even'}
        %% Load and Preprocess Raw Gabor Responses
        fprintf('Loading and preprocessing raw data for "%s" model...\n', gabor_type);
        
        % --- Load Data ---
        data_folder = '..\..\Gabor_filters_response_data\gabor_responses_3_scale_16_orientation';
        orginal_response = cell(1, 3);
        
        ac_data = load(fullfile(data_folder, gabor_type, 'ac.mat'));
        orginal_response{1} = ac_data.center_resps;
        ec_data = load(fullfile(data_folder, gabor_type, 'ec.mat'));
        orginal_response{2} = ec_data.center_resps;
        ex_data = load(fullfile(data_folder, gabor_type, 'ex.mat'));
        orginal_response{3} = ex_data.center_resps;
        
        bad_channels = load(fullfile(data_folder, gabor_type, 'bad_channels.mat'));
        bad_channels = bad_channels.all_bad_channels;
        
        % --- Remove Bad Channels ---
        for i = 1:3
            orginal_response{i}(:, bad_channels) = [];
        end
        
        %% Process Responses and Generate Trials
        processed_response = cell(1, 3);
        if any(strcmp(gabor_type, {'odd', 'even'}))
            % --- Half-Wave Rectification for Odd/Even Filters ---
            % This step models the separate ON (positive) and OFF (negative)
            % subfields of simple cells by splitting each filter's response
            % into two non-negative channels.
            fprintf('Applying half-wave rectification...\n');
            for r = 1:3
                oneRenderingData = orginal_response{r};
                [~, neuron_num] = size(oneRenderingData);
                hwrOneRenderingData = [];
                for n = 1:neuron_num
                    pos_oneNeuronData = oneRenderingData(:, n);
                    neg_oneNeuronData = oneRenderingData(:, n);
                    
                    pos_oneNeuronData(pos_oneNeuronData < 0) = 0; % Positive (ON) part
                    
                    neg_oneNeuronData(neg_oneNeuronData > 0) = 0; % Negative (OFF) part
                    neg_oneNeuronData = -neg_oneNeuronData;
                    
                    hwrOneRenderingData = [hwrOneRenderingData, pos_oneNeuronData, neg_oneNeuronData];
                end
                processed_response{r} = hwrOneRenderingData;
            end
        else % 'complex' case
            processed_response = orginal_response;
        end
        
        % --- Generate Trials using Univariate Poisson Process ---
        fprintf('Generating synthetic trials using Poisson process...\n');
        rng(1); % for reproducibility
        repeat = 10;
        response_array = cell(1, 3);
        for i = 1:3
            temp = processed_response{i};
            temp = num2cell(temp);
            temp_array = cellfun(@(x) poissrnd(x, repeat, 1), temp, 'UniformOutput', false);
            response_array{i} = temp_array;
        end
        
        %% Save Generated Trial Data
        save_path = fullfile('..','..','Gabor_filters_response_data', 'sampled_trials', gabor_type);
        if ~exist(save_path, 'dir'), mkdir(save_path); end
        save(fullfile(save_path, [gabor_type, '.mat']), "response_array");
        fprintf('Saved trial data for "%s" model to %s.\n', gabor_type, save_path);

    case 'even and odd combined'
        %% Assemble Combined Model from Odd and Even Components
        % This step loads the previously generated 'odd' and 'even' trial
        % data and concatenates them to form the complete combined model.
        fprintf('Assembling "even and odd combined" model...\n');
        
        % --- Load Odd and Even Trial Data ---
        % Using an absolute path as in the original script.
        odd_trials_path = '..\..\Gabor_filters_response_data\sampled_trials\odd\odd.mat';
        even_trials_path = '..\..\Gabor_filters_response_data\sampled_trials\even\even.mat';
        
        odd_trials = load(odd_trials_path);
        odd_trials = odd_trials.response_array;
        even_trials = load(even_trials_path);
        even_trials = even_trials.response_array;
        
        % --- Concatenate Odd and Even Filters ---
        % The filter responses are combined horizontally to create a single,
        % larger population.
        response_array = cellfun(@horzcat, odd_trials, even_trials, 'UniformOutput', false);
        
        %% Save Assembled Combined Data
        save_path = fullfile('..','..','Gabor_filters_response_data', 'sampled_trials', gabor_type);
        if ~exist(save_path, 'dir'), mkdir(save_path); end
        save(fullfile(save_path, 'even and odd combined.mat'), "response_array");
        fprintf('Saved assembled "even and odd combined" model data to %s.\n', save_path);
        
    otherwise
        error("Invalid gabor_type. Please choose 'complex', 'odd', 'even', or 'even and odd combined'.");
end