%{
% FILENAME: generate_ann_trial_data.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-10
%
% DESCRIPTION:
%   Processes mean activation data from specified ANN layers (VGG/AlexNet).
%   It loads raw activations for AC, EC, and EX stimuli, removes unresponsive
%   "bad channels", and generates synthetic trial-by-trial data using a
%   Poisson process to mimic neuronal variability. The output is a MAT-file
%   ready for decoding analysis.
% INPUTS:
%   - .mat files containing mean ANN activations for each stimulus type,
%     located in a path like '.\<network>\pool<pool_layer>\'.
%   - 'bad_channels.mat': A file listing indices of units to be excluded.
%
% OUTPUTS:
%   - A raw data file: '\<network>\<network>_layer<pool_layer>_raw.mat'
%     containing the cleaned mean activations.
%   - A final data file: '..\..\ANN_response_data\<network>\<network>_layer<pool_layer>.mat'
%     containing the generated trial-by-trial data.
%}
%% INITIALIZATION
clc;
clear;

network = 'vgg_pooling_resps'; % Options: 'alexnet_pooling_resps', 'vgg_pooling_resps'
pool_layer = '2';              % Specify the pooling layer number as a string

%% LOAD AND PREPROCESS RAW ANN ACTIVATIONS
% =========================================================================
% This section loads the mean activation data for each stimulus type,
% removes specified "bad channels" (e.g., non-responsive units), and
% saves the result as an intermediate "raw" file.
% =========================================================================

% --- Load the original mean activation data ---
fprintf('Loading and preprocessing data for %s, layer %s...\n', network, pool_layer);
orginal_response = cell(1, 3);

% Load Appearance Contour (AC) data
ac_data_path = fullfile('..','..','ANN_response_data', network, ['pool' pool_layer], 'ac.mat');
ac_data = load(ac_data_path);
orginal_response{1} = ac_data.center_resps;

% Load Edge Contour (EC) data
ec_data_path = fullfile('..','..','ANN_response_data', network, ['pool' pool_layer], 'ec.mat');
ec_data = load(ec_data_path);
orginal_response{2} = ec_data.center_resps;

% Load Example (EX) data
ex_data_path = fullfile('..','..','ANN_response_data', network, ['pool' pool_layer], 'ex.mat');
ex_data = load(ex_data_path);
orginal_response{3} = ex_data.center_resps;

% Load the indices of bad channels to be removed
bad_channels_path = fullfile('..','..','ANN_response_data', network, ['pool' pool_layer], 'bad_channels.mat');
bad_channels = load(bad_channels_path);
bad_channels = bad_channels.all_bad_channels;

% --- Remove the bad channels from each stimulus condition ---
for i = 1:3
    orginal_response{i}(:, bad_channels) = [];
end
response_array = orginal_response;

% --- Save the cleaned, preprocessed data ---
raw_save_path = fullfile('..','..','ANN_response_data', network, [network, '_layer', pool_layer, '_raw.mat']);
save(raw_save_path, "response_array");
fprintf('Saved preprocessed raw data to %s\n', raw_save_path);


%% GENERATE SYNTHETIC TRIAL DATA USING A POISSON PROCESS
% =========================================================================
% This section simulates trial-by-trial variability from the mean ANN
% activations. This is a standard method to make ANN data comparable to
% real neuronal data for decoding analyses.
% =========================================================================

clc;
rng(1); % Set random seed for reproducibility

% --- Load the combined raw data ---
% This step is slightly redundant but ensures the next block is independent.
raw_data_path = ['..\..\ANN_response_data\', network, '\', network, '_layer', pool_layer, '_raw.mat'];
orginal_response_data = load(raw_data_path);
orginal_response = orginal_response_data.response_array;

% --- Generate responses using a univariate Poisson distribution ---
fprintf('Generating synthetic trials using Poisson process...\n');
repeat = 10; % Number of synthetic trials to generate for each stimulus
response_array = cell(1, 3);

for i = 1:3 % Loop through each stimulus condition (AC, EC, EX)
    % Get the mean activations for the current condition
    temp = orginal_response{i};
    temp = num2cell(temp);
    
    % For each unit's mean activation (x), generate 'repeat' number of
    % trials by drawing from a Poisson distribution: poissrnd(x, repeat, 1).
    temp_array = cellfun(@(x) poissrnd(x, repeat, 1), temp, 'UniformOutput', false);
    response_array{i} = temp_array;
end

% --- Save the final dataset with synthetic trials ---
final_save_path = ['..\..\ANN_response_data\', network, '\', network, '_layer', pool_layer, '.mat'];
save(final_save_path, "response_array");
fprintf('Saved final data with synthetic trials to %s\n', final_save_path);
