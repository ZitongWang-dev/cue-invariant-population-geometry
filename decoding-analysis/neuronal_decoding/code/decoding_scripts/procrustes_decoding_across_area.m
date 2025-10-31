%{
Filename: procrustes_decoding_across_area.m
Author: Zitong Wang
Date: 2025-07-09

Description:
    Hierarchical Procrustes transfer decoding between areas V1 and V2. A
    user‑selectable flag (`decode_direction`) chooses V2→V1 or V1→V2.
    The script loads V1/V2 data, prepares trials, then runs decoding for each
    stimulus pair.

Usage:
    >> decode_direction = 'V2toV1'; run procrustes_decoding_across_area  % default
    >> decode_direction = 'V1toV2'; run procrustes_decoding_across_area  % flip direction

User-configurable parameters:
    decode_direction       – 'V2toV1' | 'V1toV2'
    monkey                 – 'FR' | 'KO'
    timewindow             – [startMS endMS] spike-count window
    neuron_sample_repeat   – # random neuron samples per neuron pool
    trial_sample_repeat    – # trial shuffles per sample
    rng(1)                 – seed for reproducibility
%}

%% House‑keeping & parameters
clc; clear;

decode_direction     = 'V2toV1';   % 'V2toV1' or 'V1toV2'
monkey               = 'KO';       % 'FR' or 'KO'
timewindow           = [330 630];  % spike-count window (ms)

neuron_sample_repeat = 50;         % # random neuron samples
trial_sample_repeat  = 15;         % # trial shuffles per sample

rng(1);  % reproducibility

%% Paths & data loading
root_data = fullfile('..','..','neuronal_data',monkey);
root_out  = fullfile('..','..','results','decoding_outputs', ...
    'procrustes_decoding_across_area_results',monkey);
load_v1 = load(fullfile(root_data,'V1',[monkey '_V1_allstim.mat']),'three_stim_array');
load_v2 = load(fullfile(root_data,'V2',[monkey '_V2_allstim.mat']),'three_stim_array');
spike_data_v1 = load_v1.three_stim_array;
spike_data_v2 = load_v2.three_stim_array;
clear load_v1 load_v2;

%% Prepare trial data
nStim = numel(spike_data_v1);
[data_trial_v1,label_trial_v1] = deal(cell(1,nStim));
[data_trial_v2,label_trial_v2] = deal(cell(1,nStim));
for i = 1:nStim
    [data_trial_v1{i},label_trial_v1{i}] = multiclass_svmloader_PT(spike_data_v1{i},timewindow);
    [data_trial_v2{i},label_trial_v2{i}] = multiclass_svmloader_PT(spike_data_v2{i},timewindow);
end
labels = label_trial_v1{1};  % same labels

%% Neuron‑count list (incremental analysis)
neuron_num_coode = { [48]; [109] };
name2neuronidx   = struct('FR',1,'KO',2);
neuron_num_list  = neuron_num_coode{name2neuronidx.(monkey)};

%% Decode based on direction
switch decode_direction
    case 'V2toV1'
        save_path = fullfile(root_out,'V1 decoding V2');
        if ~exist(save_path,'dir'), mkdir(save_path); end
        % V2→V1: source=v2, target=v1
        fprintf('Running V2→V1 decoding...\n');
        % target templates: V1-ac, V1-ec, V1-ex
        targets = {'v1ac','v1ec','v1ex'};
        sources = {'v2ac','v2ec','v2ex'};
        for t = 1:3
            tgt = targets{t};
            for s = 1:3
                src = sources{s};
                fname = sprintf('%s_%s_results.mat',src,tgt);
                fprintf('Decoding %s→%s\n',upper(src),upper(tgt));
                results = hier_decoding(src,tgt,data_trial_v1,data_trial_v2,labels, ...
                            neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
                save(fullfile(save_path,fname),'results');
            end
        end

    case 'V1toV2'
        save_path = fullfile(root_out,'V2 decoding V1');
        if ~exist(save_path,'dir'), mkdir(save_path); end
        % V1→V2: source=v1, target=v2
        fprintf('Running V1→V2 decoding...\n');
        targets = {'v2ac','v2ec','v2ex'};
        sources = {'v1ac','v1ec','v1ex'};
        for t = 1:3
            tgt = targets{t};
            for s = 1:3
                src = sources{s};
                fname = sprintf('%s_%s_results.mat',src,tgt);
                fprintf('Decoding %s→%s\n',upper(src),upper(tgt));
                results = hier_decoding(src,tgt,data_trial_v1,data_trial_v2,labels, ...
                            neuron_num_list,neuron_sample_repeat,trial_sample_repeat);
                save(fullfile(save_path,fname),'results');
            end
        end

    otherwise
        error('decode_direction must be ''V2toV1'' or ''V1toV2''');
end
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