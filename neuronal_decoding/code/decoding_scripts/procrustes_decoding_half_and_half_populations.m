%{
Filename: procrustes_decoding_half_and_half_populations.m
Author: Zitong Wang
Date: 2025-06-24

Description:
    Performs half-and-half Procrustes transfer decoding on neuronal datasets.
    For each stimulus condition (ac, ec, ex), neurons are split into two equal-sized
    groups.The transform parameters (scale, rotation, translation) are computed on averaged
    trial blocks, applied to test trials of group 2, and decoding accuracy is evaluated
    using an SVM trained on group 1. Repeated random splits and trial shuffles quantify
    variability under full, partial, and randomized transformations.

Parameters to edit:
    monkey = 'FR';     % 'FR' or 'KO'
    vp     = 'V1';     % 'V1' or 'V2'
    timewindow = [330 630];  % spike-count window (ms)
    neuron_sample_repeat = 50;  % number of neuron splits
    trial_sample_repeat  = 15;  % trial shuffles per split
    rng(1);  % for reproducibility

Inputs:
    neuronal_data/<monkey>/<vp>/<monkey>_<vp>_allstim.mat
        - variable 'three_stim_array': 1×3 cell of [trials × neurons] matrices

Outputs:
    results/decoding_outputs/procrustes_decoding_half_and_half_populations_results/
        <monkey>/<vp>/{ac,ec,ex}_results.mat
%}

%% Initialize environment
clc; clear;

%% Configuration
monkey = 'FR';
vp = 'V1';
timewindow = [330 630];
neuron_sample_repeat = 50;
trial_sample_repeat = 15;
rng(1);

%% Load spike data for three stimuli
spike_data = load(fullfile('..','..','neuronal_data',monkey,vp,...
    sprintf('%s_%s_allstim.mat',monkey,vp)), 'three_stim_array').three_stim_array;

%% Prepare results directory
save_path = fullfile('..','..','results','decoding_outputs', ...
    'procrustes_decoding_half_and_half_populations_results', monkey, vp);
if ~exist(save_path,'dir'), mkdir(save_path); end

%% Convert spike counts to trial-by-neuron matrices and labels
for i = 1:numel(spike_data)
    [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, timewindow); %#ok<SAGROW>
end
labels = label_trial{1};

%% Determine neuron counts for this monkey/area
neuron_num_code = {[48]; [48 112]; [109]; [109 146]};
map = struct('FRV1',1,'FRV2',2,'KOV1',3,'KOV2',4);
neuron_num_list = neuron_num_code{map.(strcat(monkey,vp))};

%% Run half-and-half decoding for each stimulus
for sti = {'ac','ec','ex'}
    stim_name = sti{1};
    fprintf('Decoding %s for %s %s\n', stim_name, monkey, vp);
    results = half_decoding(stim_name, data_trial, labels, neuron_num_list, ...
                             neuron_sample_repeat, trial_sample_repeat);
    save(fullfile(save_path, [stim_name '_results.mat']), 'results');
end


%%
function result_one_cell = half_decoding(stim_name,data_trial,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat)
% stim1 is transformed stim2 is the training data

name2idx = struct('ac',1,'ec',2,'ex',3);
data_onestim = data_trial{name2idx.(stim_name)};
neuron_num = max(neuron_num_list);
% sampled neuron numbers are set to be half of the V1's neuron number for
% comparison
group_num = floor(min(neuron_num_list)/2);

result_one_cell = [];
parfor neuron_repeat = 1:neuron_sample_repeat
    % sample different neurons, for neuron_sample_repeat times
    [sampled_neuron_idx1,sampled_neuron_idx2] = sample_neuron(neuron_num,group_num);
    data_sample1 = data_onestim(:,sampled_neuron_idx1);
    data_sample2 = data_onestim(:,sampled_neuron_idx2);
    data_sample1 = zscore(data_sample1);
    data_sample2 = zscore(data_sample2);
    % sample trials multiple times
    trial_result = zeros(trial_sample_repeat,8);
    for trial_repeat =1:trial_sample_repeat
        [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1] = data_trial_sampler(data_sample1,data_sample2,labels);
        [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_noscale,accuracy_norotation,accuracy_notranslation,accuracy_non_transfer_control] = ...
            pro_decoding(training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1);
        
        trial_result(trial_repeat,:) = [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_noscale,accuracy_norotation,accuracy_notranslation,accuracy_non_transfer_control];

    end
    result_one_cell = [result_one_cell;trial_result];
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

function [sampled_neuron_idx1,sampled_neuron_idx2] = sample_neuron(neuron_num,group_num)
% sample two equal groups of neurons

original_idx = 1:neuron_num;
sampled_neuron_idx = randperm(neuron_num,group_num*2);
sampled_neuron_idx1 = original_idx(sampled_neuron_idx(1:group_num));
sampled_neuron_idx2 = original_idx(sampled_neuron_idx(group_num+1:end));

end

function [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_noscale,accuracy_norotation,accuracy_notranslation,accuracy_non_transfer_control] ...
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
accuracy_noscale= sum(stim1_test_label == predicted_Label_scale)/length(predicted_Label_scale);
% rotation only
predicted_Label_rotation = predict(stim2_model,partially_transformed_stim1{1,2});
accuracy_norotation= sum(stim1_test_label == predicted_Label_rotation)/length(predicted_Label_rotation);
% translation only
predicted_Label_translation = predict(stim2_model,partially_transformed_stim1{1,3});
accuracy_notranslation= sum(stim1_test_label == predicted_Label_translation)/length(predicted_Label_translation);
% non-transformation transfer decoding control
accuracy_non_transfer_control = sum(randperm(50)' == predicted_Label_s1)/length(predicted_Label_s1);
end