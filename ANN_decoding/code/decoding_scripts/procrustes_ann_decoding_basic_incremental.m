%{
% FILENAME: procrustes_ann_decoding_basic_incremental.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Performs an incremental Procrustes-based transfer decoding on synthetic
%   ANN data (AlexNet/VGG). This version tests how decoding accuracy scales
%   with the number of units ('neurons') in the population. It runs the
%   analysis on progressively larger subsets of units for all stimulus
%   pairs and saves the results.
%}

%% INITIALIZATION
clc;
clear;

%% CONFIGURATION
% Specify the network and layer for analysis.
net = 'vgg'; % 'alexnet' or 'vgg'
layer = '4'; % Pooling layer number as a string

%% LOAD AND PREPARE DATA
% Load the synthetic ANN data, prepare save directories, and format data
% for the incremental decoding functions.
fprintf('Loading data for %s, layer %s...\n', net, layer);

% --- Load synthetic trial data ---
data_path = fullfile('..', '..', 'ANN_response_data', [net, '_pooling_resps'], [net, '_pooling_resps_layer', layer, '.mat']);
spike_data = load(data_path, 'response_array');
spike_data = spike_data.response_array;

% --- Prepare output directory ---
save_path = fullfile('..', '..', 'results', 'decoding_outputs', 'procrustes_ann_decoding_basic_incremental_results', net, ['pool', layer]);
if ~exist(save_path, 'dir'), mkdir(save_path); end

% --- Format data for decoding ---
% Convert the cell array of trials into a standard [trial x neuron] matrix.
data_trial = cellfun(@(x) cell2mat(x), spike_data, 'UniformOutput', false);

% Define stimulus labels (50 stimuli, 10 trials each).
labels = reshape(repmat((1:50)', 1, 10)', [], 1);

% --- Create an incremental list of neuron counts for analysis ---
% This list determines the different population sizes that will be tested.
max_neuron_num = size(data_trial{1}, 2);
neuron_num_list = [];
n = 1;
% The list grows exponentially (5, 10, 20, 40...) to efficiently sample
% the effect of population size on decoding accuracy.
while (5 * 2^(n-1)) < max_neuron_num
    neuron_num_list = [neuron_num_list, 5 * 2^(n-1)];
    n = n + 1;
end
% Ensure the full population size is always included as the last step.
neuron_num_list = unique([neuron_num_list, max_neuron_num]);

%% DECODING PARAMETERS
% Set parameters for the decoding analysis, including the number of
% repetitions for robust sampling.
neuron_sample_repeat = 20; % Number of times to resample neuron subsets
trial_sample_repeat = 15;  % Number of times to resample trial partitions
rng(1); % Set the random seed for reproducibility

%% RUN INCREMENTAL DECODING ANALYSIS
% This section runs the incremental Procrustes decoding for all six ordered
% stimulus pairs. The main function will loop through the neuron_num_list.
tic;
fprintf('Starting incremental decoding analysis...\n');

% --- Run decoding for each stimulus pair ---
acec_results = ann_incremental_decoding('ac', 'ec', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'acec_results.mat'), 'acec_results');
fprintf('Saved acec_results.mat\n');

ecex_results = ann_incremental_decoding('ec', 'ex', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'ecex_results.mat'), 'ecex_results');
fprintf('Saved ecex_results.mat\n');

acex_results = ann_incremental_decoding('ac', 'ex', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'acex_results.mat'), 'acex_results');
fprintf('Saved acex_results.mat\n');

ecac_results = ann_incremental_decoding('ec', 'ac', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'ecac_results.mat'), 'ecac_results');
fprintf('Saved ecac_results.mat\n');

exec_results = ann_incremental_decoding('ex', 'ec', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'exec_results.mat'), 'exec_results');
fprintf('Saved exec_results.mat\n');

exac_results = ann_incremental_decoding('ex', 'ac', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'exac_results.mat'), 'exac_results');
fprintf('Saved exac_results.mat\n');

toc;
fprintf('Analysis complete for %s - layer %s.\n', net, layer);

%%
function results = ann_incremental_decoding(stim1,stim2,data_trial,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat)
% stim1 is transformed; stim2 is the training data

pair_wise_data_trial = pair_pcaloader(data_trial,stim1,stim2);

%load data
stim1_data = pair_wise_data_trial(1:500,:);
stim2_data = pair_wise_data_trial(501:1000,:);
% zscore
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);
results = cell(1,length(neuron_num_list)); % each cell stores results for one neuron-number decoding

for neuron_squence = 1:length(neuron_num_list)
    disp([stim1,stim2,num2str(neuron_squence)])
    result_one_cell = [];
    parfor neuron_repeat = 1:neuron_sample_repeat
        % sample different neuron, for neuron_sample_repeat times
        sampled_neuron_idx = sample_neuron(neuron_num_list,neuron_squence);
        stim1_data_sample = stim1_data(:,sampled_neuron_idx);
        stim2_data_sample = stim2_data(:,sampled_neuron_idx);
        % sample trials multiple times
        trial_result = zeros(trial_sample_repeat,8);
        for trial_repeat =1:trial_sample_repeat
            [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1] = data_trial_sampler(stim1_data_sample,stim2_data_sample,labels);
            [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control] = pro_decoding(training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1);
            
            trial_result(trial_repeat,:) = [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control];
    
        end
        result_one_cell = [result_one_cell;trial_result];
    end
    results{neuron_squence} = result_one_cell;
end

end

function stim_data_trial_averged = take_average(stim_data,number_of_average)
[trial,neuron] = size(stim_data);
stim_data_trial_averged = zeros(50,neuron);
for i = 1:50
    temp = stim_data(i*number_of_average-(number_of_average-1):i*number_of_average,:);
    temp_mean = mean(temp,1);
    stim_data_trial_averged(i,:) = temp_mean;
end
end

function [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1] = data_trial_sampler(stim1_data_sample,stim2_data_sample,labels)
% sample trial
tial_out_of_ten = 1;
valid_trial_idx = sample_trial(tial_out_of_ten); % 10 fold
training_trial_idx = setdiff([1:500],valid_trial_idx);

% SVM training data, all stim2 data
training_data = stim2_data_sample;
training_label = labels;

% test data 
stim1_test_data = stim1_data_sample(valid_trial_idx,:);
stim1_test_label = labels(valid_trial_idx,:);

% transformed data
stim2_target_data = stim2_data_sample;
stim1_2transform_data = stim1_data_sample(training_trial_idx,:);
% take the mean
stim2_target_data_averged = take_average(stim2_target_data,10);
stim1_2transform_data_averged = take_average(stim1_2transform_data,10 - tial_out_of_ten);

[d,Z,transform] = procrustes(stim2_target_data_averged,stim1_2transform_data_averged);
t_matrix = transform.T;
b = transform.b;
c = transform.c;
transformed_stim1 = b*stim1_test_data*t_matrix+c; 

% partially transformed data - without 1 transofrm
% without_scaling_transformed_stim1 = stim1_test_data*t_matrix+c;
% without_rotation_transformed_stim1 = b*stim1_test_data+c;
% without_translation_transformed_stim1 = b*stim1_test_data*t_matrix;
% partially_transformed_stim1 = {without_scaling_transformed_stim1,without_rotation_transformed_stim1,without_translation_transformed_stim1};

% partially transformed data - only 1 transofrm
scaling_only_transformed_stim1 = b*stim1_test_data;
rotation_only_transformed_stim1 = stim1_test_data*t_matrix;
translation_only_transformed_stim1 = stim1_test_data + c;
partially_transformed_stim1 = {scaling_only_transformed_stim1,rotation_only_transformed_stim1,translation_only_transformed_stim1};

% randomlized data
shuffle_idx = randperm(size(stim1_2transform_data_averged,1));
toy_averaged = stim1_2transform_data_averged(shuffle_idx,:);
[d_fake,Z_fake,t_fake] = procrustes(stim2_target_data_averged,toy_averaged);
t_fakematrix = t_fake.T;
toy_b = t_fake.b;
toy_c = t_fake.c;
toy_transd = toy_b * stim1_test_data*t_fakematrix + toy_c;


% fake, transform, validation has the same label
end


function sample_trial_idx = sample_trial(trialnum_outof_ten)
% trial sampling
sample_trial_matrix = zeros(50,10); % same for all cases
for i = 1:50
    sample_trial_matrix(i,:) = randperm(10);
end
sample_trial_4condition = sample_trial_matrix(:,1:trialnum_outof_ten); % sampled trial num for each condition
base = [0:10:490]';
sample_trial_idx = sample_trial_4condition+base;
sample_trial_idx = reshape(sample_trial_idx,[],1);
sample_trial_idx = sort(sample_trial_idx);
end

function sampled_neuron_idx = sample_neuron(neuron_num_list,sequence)
% neuron sampling
% sequence: the idx of how many neurons to be sampled
origin_neuron_num = max(neuron_num_list);
sampled_neuron_num = neuron_num_list(sequence);
sampled_neuron_idx = randperm(origin_neuron_num,sampled_neuron_num);
disp([origin_neuron_num,sampled_neuron_num])
end

function pair_wise_pca_data = pair_pcaloader(dca_data,stim1,stim2)
% use the input stimulus name to load pair_wise dca data

name2idx = struct('ac',1,'ec',2,'ex',3);

stim1_data = dca_data{name2idx.(stim1)};
stim2_data = dca_data{name2idx.(stim2)};
pair_wise_pca_data = [stim1_data;stim2_data];
end

function [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control] ...
    = pro_decoding(training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1)
% get 4 accs
stim2_model = fitcecoc(training_data,training_label);

CVMdl1 = crossval(stim2_model,'KFold',10);
genError1 = kfoldLoss(CVMdl1);
genAcc = 1-genError1;

% non-transformation transfer decoding
predicted_Label_s1 = predict(stim2_model,stim1_test_data);
accuracy_s1 = sum(stim1_test_label == predicted_Label_s1)/length(predicted_Label_s1);
% with-transformation transfer decoding
predicted_Label_transformed_s1 = predict(stim2_model,transformed_stim1);
accuracy_transformed_s1 = sum(stim1_test_label == predicted_Label_transformed_s1)/length(predicted_Label_transformed_s1);
% random-transformation transfer decoding
predicted_Label_rand = predict(stim2_model,toy_transd);
accuracy_rand = sum(stim1_test_label == predicted_Label_rand)/length(predicted_Label_rand);
% partial-transformation transfer decoding
% scaling only
predicted_Label_scale = predict(stim2_model,partially_transformed_stim1{1,1});
accuracy_onlyscale= sum(stim1_test_label == predicted_Label_scale)/length(predicted_Label_scale);
% rotation only
predicted_Label_rotation = predict(stim2_model,partially_transformed_stim1{1,2});
accuracy_onlyrotation= sum(stim1_test_label == predicted_Label_rotation)/length(predicted_Label_rotation);
% translation only
predicted_Label_translation = predict(stim2_model,partially_transformed_stim1{1,3});
accuracy_onlytranslation= sum(stim1_test_label == predicted_Label_translation)/length(predicted_Label_translation);
% non-transformation transfer decoding control
accuracy_non_transfer_control = sum(randperm(50)' == predicted_Label_s1)/length(predicted_Label_s1);
end