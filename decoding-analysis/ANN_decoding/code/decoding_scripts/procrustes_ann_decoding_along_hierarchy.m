%{
% FILENAME: procrustes_ann_decoding_along_hierarchy.m
% AUTHOR:   Zitong Wang
% DATE:     2025-07-13
%
% DESCRIPTION:
%   Performs a cross-layer Procrustes transfer decoding on ANN data.
%   This analysis tests whether the geometric structure of activations in one
%   layer can be aligned with the structure in another, higher layer. This
%   mimics the cross-cortical area (e.g., V1-to-V2) decoding performed on
%   real neuronal data. Results for all cue combinations are saved.
%}

%% INITIALIZATION
clc;
clear;

%% CONFIGURATION
% Specify the network and the two layers for hierarchical analysis.
net = 'alexnet';            % 'alexnet' or 'vgg'
lower_layer_idx = '2';      % Index for the lower-level ANN layer
higher_layer_idx = '5';     % Index for the higher-level ANN layer

%% LOAD AND PREPARE DATA
% Load synthetic ANN data from two different layers and prepare it for
% the cross-layer decoding analysis.
fprintf('Loading data for %s, layers %s and %s...\n', net, lower_layer_idx, higher_layer_idx);

% --- Load synthetic trial data for both layers ---
spike_data_lower = load(fullfile('..','..','ANN_response_data',[net,'_pooling_resps'],[net,'_pooling_resps_layer',lower_layer_idx,'.mat']),'response_array');
spike_data_lower = spike_data_lower.response_array;
spike_data_higher = load(fullfile('..','..','ANN_response_data',[net,'_pooling_resps'],[net,'_pooling_resps_layer',higher_layer_idx,'.mat']),'response_array');
spike_data_higher = spike_data_higher.response_array;

% --- Prepare output directory ---
save_path = fullfile('..','..','results','decoding_outputs','procrustes_ann_decoding_along_hierarchy_results',net,[net,' POOL ',lower_layer_idx,'&',higher_layer_idx],'Higher decoding Lower');
if ~exist(save_path, 'dir'), mkdir(save_path); end

% --- Format data for decoding ---
data_trial_lower = cellfun(@(x) cell2mat(x),spike_data_lower,'UniformOutput',false);
data_trial_higher = cellfun(@(x) cell2mat(x),spike_data_higher,'UniformOutput',false);

% --- Define neuron count for analysis ---
% To ensure a fair comparison, the population size is set to the minimum
% of the two layers' unit counts.
neuron_num_list = min(size(data_trial_lower{1},2),size(data_trial_higher{1},2));

% Define stimulus labels (50 stimuli, 10 trials each).
labels = reshape(repmat((1:50)',1,10)',[],1);

%% DECODING PARAMETERS
% Set parameters for the decoding analysis.
neuron_sample_repeat = 50;
trial_sample_repeat = 15;
rng(1); % Set random seed for reproducibility

% NOTE: The following section uses strings like 'v1ac' and 'v2ec' as
% arguments. In this context, 'v1' refers to the 'lower_layer' data and
% 'v2' refers to the 'higher_layer' data.

%% RUN HIERARCHICAL DECODING (Higher-to-Lower)
% This section trains a decoder on the lower layer and tests its ability
% to decode activity from the higher layer after Procrustes alignment.
fprintf('Running Higher-to-Lower decoding...\n');
tic;

% --- Train on Lower Layer AC ---
v2acv1ac_results = hier_decoding('v2ac','v1ac',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2acv1ac_results.mat'),'v2acv1ac_results');
v2ecv1ac_results = hier_decoding('v2ec','v1ac',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2ecv1ac_results.mat'),'v2ecv1ac_results');
v2exv1ac_results = hier_decoding('v2ex','v1ac',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2exv1ac_results.mat'),'v2exv1ac_results');

% --- Train on Lower Layer EC ---
v2acv1ec_results = hier_decoding('v2ac','v1ec',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2acv1ec_results.mat'),'v2acv1ec_results');
v2ecv1ec_results = hier_decoding('v2ec','v1ec',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2ecv1ec_results.mat'),'v2ecv1ec_results');
v2exv1ec_results = hier_decoding('v2ex','v1ec',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2exv1ec_results.mat'),'v2exv1ec_results');

% --- Train on Lower Layer EX ---
v2acv1ex_results = hier_decoding('v2ac','v1ex',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2acv1ex_results.mat'),'v2acv1ex_results');
v2ecv1ex_results = hier_decoding('v2ec','v1ex',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2ecv1ex_results.mat'),'v2ecv1ex_results');
v2exv1ex_results = hier_decoding('v2ex','v1ex',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v2exv1ex_results.mat'),'v2exv1ex_results');
toc;

%% RUN HIERARCHICAL DECODING (Lower-to-Higher)
% This section performs the reverse analysis: training on the higher layer
% and testing on the lower layer.
fprintf('Running Lower-to-Higher decoding...\n');
tic;

% --- Train on Higher Layer AC ---
v1acv2ac_results = hier_decoding('v1ac','v2ac',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1acv2ac_results.mat'),'v1acv2ac_results');
v1ecv2ac_results = hier_decoding('v1ec','v2ac',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1ecv2ac_results.mat'),'v1ecv2ac_results');
v1exv2ac_results = hier_decoding('v1ex','v2ac',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1exv2ac_results.mat'),'v1exv2ac_results'); 

% --- Train on Higher Layer EC ---
v1acv2ec_results = hier_decoding('v1ac','v2ec',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1acv2ec_results.mat'),'v1acv2ec_results');
v1ecv2ec_results = hier_decoding('v1ec','v2ec',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1ecv2ec_results.mat'),'v1ecv2ec_results');
v1exv2ec_results = hier_decoding('v1ex','v2ec',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1exv2ec_results.mat'),'v1exv2ec_results');

% --- Train on Higher Layer EX ---
v1acv2ex_results = hier_decoding('v1ac','v2ex',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1acv2ex_results.mat'),'v1acv2ex_results');
v1ecv2ex_results = hier_decoding('v1ec','v2ex',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1ecv2ex_results.mat'),'v1ecv2ex_results');
v1exv2ex_results = hier_decoding('v1ex','v2ex',data_trial_lower,data_trial_higher,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
save(fullfile(save_path,'v1exv2ex_results.mat'),'v1exv2ex_results');
toc;

%%
function results = hier_decoding(stim1,stim2,data_trial_v1,data_trial_v2,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat)
% stim1 is transformed; stim2 is the target
%load data
[stim1_data,stim2_data] = stim_data_loader(data_trial_v1,data_trial_v2,stim1,stim2);

% zscore
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);

results = cell(1,length(neuron_num_list)); % each cell stores results for one neuron-number decoding

disp([stim1,stim2])
result_one_cell = [];
parfor neuron_repeat = 1:neuron_sample_repeat
    neuron_repeat
    % sample different neuron, for neuron_sample_repeat times
    sampled_neuron_idx1 = randperm(size(stim1_data,2),neuron_num_list);
    sampled_neuron_idx2 = randperm(size(stim2_data,2),neuron_num_list);
    stim1_data_sample = stim1_data(:,sampled_neuron_idx1);
    stim2_data_sample = stim2_data(:,sampled_neuron_idx2);
    % sample trials multiple times
    trial_result = zeros(trial_sample_repeat,8);
    for trial_repeat =1:trial_sample_repeat
        [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1] = data_trial_sampler(stim1_data_sample,stim2_data_sample,labels);
        [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_noscale,accuracy_norotation,accuracy_notranslation,accuracy_non_transfer_control] = pro_decoding(training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,toy_transd,partially_transformed_stim1);
            
        trial_result(trial_repeat,:) = [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_noscale,accuracy_norotation,accuracy_notranslation,accuracy_non_transfer_control];

    end
    result_one_cell = [result_one_cell;trial_result];
end
results{1} = result_one_cell;

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

function [stim1_data,stim2_data] = stim_data_loader(data_v1,data_v2,stim1,stim2)
% use the input stimulus name to load pair_wise dca data

name2idx = struct('ac',1,'ec',2,'ex',3);
% load data according to the name
if stim1(1:2) == 'v1'
    stim1_data = data_v1{name2idx.(stim1(3:4))};
elseif stim1(1:2) == 'v2'
    stim1_data = data_v2{name2idx.(stim1(3:4))};
end

if stim2(1:2) == 'v1'
    stim2_data = data_v1{name2idx.(stim2(3:4))};
elseif stim2(1:2) == 'v2'
    stim2_data = data_v2{name2idx.(stim2(3:4))};
end

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

function stim_data_trial_averged = take_average(stim_data,number_of_average)
[trial,neuron] = size(stim_data);
stim_data_trial_averged = zeros(50,neuron);
for i = 1:50
    temp = stim_data(i*number_of_average-(number_of_average-1):i*number_of_average,:);
    temp_mean = mean(temp,1);
    stim_data_trial_averged(i,:) = temp_mean;
end
end