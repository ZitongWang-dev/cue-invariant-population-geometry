%{
Filename: procrustes_decoding_basic_trialShuffledWithMean.m
Author: Zitong Wang
Date: 2026-07-09

Description:
    Structure-destroying (correspondence-shuffled) control for Procrustes-based
    transfer decoding. Unlike procrustes_decoding_basic.m, whose randomized
    control learns a Procrustes transform from a *shuffled* mean matrix and then
    applies it to the *untouched* test trials (which only asks "does a random
    rotation work" and trivially fails), this script destroys the stimulus
    correspondence in a self-consistent way:

        1. The averaged (mean) transform matrix rows are shuffled by a permutation pi.
        2. A ROTATION-only Procrustes transform is learned from this shuffled matrix
           onto the target (stim2) mean matrix. So row i of the shuffled matrix,
           which is the original condition pi(i), is aligned onto target condition i.
        3. The held-out test trials are shuffled by the SAME permutation pi and the
           learned rotation is applied to them. Thus test trial i (originally
           condition pi(i)) is mapped toward target condition i.

    Decoding is then evaluated against two label sets:
        - target labels (1:50): does PT succeed in mapping the shuffled config onto
          the target config? (i.e. "how well can PT achieve the mapping even though
          the corresponding structure is destroyed"). Expected to remain HIGH.
        - shuffled/identity labels (pi): is each trial's true pre-rotation condition
          identity preserved? Matches only at permutation fixed points -> CHANCE.

    We use ROTATION-only (T matrix only, no scaling b, no translation c) because the
    rotation-only partial transform gave the highest transfer decoding on the real
    (unshuffled) data. The true-correspondence rotation-only decoding is carried in
    the same output for a matched comparison.

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat
        Contains variable 'three_stim_array': cell array of [trial x neuron] matrices for each stimulus.
    - timewindow: two-element vector [start_ms end_ms] defining spike count window.

Outputs:
    - MAT-files saved to:
      results/decoding_outputs/procrustes_decoding_basic_trialShuffledWithMean/<monkey>/<area>/
      Named: acec_results.mat, ecex_results.mat, acex_results.mat, ecac_results.mat, exec_results.mat, exac_results.mat

    Each result cell is an [N x 6] matrix, one row per (neuron_repeat, trial_repeat).
    Columns:
      1 genAcc                 - stim2 within-condition 10-fold CV accuracy (ceiling)
      2 acc_true_rot           - TRUE correspondence, rotation-only transfer decoding (vs 1:50)
      3 acc_shuffle_target     - shuffled correspondence, rotation-only, decoded vs TARGET labels (1:50)
      4 acc_shuffle_identity   - shuffled correspondence, rotation-only, decoded vs SHUFFLED identity labels (pi)
      5 acc_notransform        - raw stim1 test into stim2 model, no transform (vs 1:50)
      6 acc_chance             - random-permutation chance baseline
%}

%%
clc; clear;

%% Configuration
monkey = 'FR';      % Monkey ID ('FR' or 'KO')
vp = 'V1';          % Visual area ('V1' or 'V2')

% Load neuronal data for all stimuli
data_file = fullfile('..','..','neuronal_data', monkey, vp, sprintf('%s_%s_allstim.mat', monkey, vp));
tmp = load(data_file, 'three_stim_array');
spike_data = tmp.three_stim_array;

% Prepare output directory
save_path = fullfile('..','..','results','decoding_outputs','procrustes_decoding_basic_trialShuffledWithMean', monkey, vp);
if ~exist(save_path, 'dir'), mkdir(save_path); end

%% Define decoding time window (ms)
timewindow = [330 630];  % e.g., early: [340 410], late: [410 480], long-late: [410 550]

%% Load trial-by-trial data and labels
num_conditions = numel(spike_data);
data_trial = cell(1, num_conditions);
label_trial = cell(1, num_conditions);
for i = 1:num_conditions
    [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, timewindow);
end

%% Determine neuron sampling list based on monkey & area
neuron_num_code = {[48]; [48 112]; [109]; [109 146]};
map = struct('FRV1',1,'FRV2',2,'KOV1',3,'KOV2',4);
idx = map.(strcat(monkey, vp));
neuron_num_list = neuron_num_code{idx};
labels = label_trial{1};

%% Decoding parameters
neuron_sample_repeat = 50;
trial_sample_repeat = 15;
rng(1);

%% Run structure-shuffled Procrustes control for all stimulus pairs
tic;
acec_results = procrustes_decoding('ac','ec', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'acec_results.mat'), 'acec_results');

ecex_results = procrustes_decoding('ec','ex', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'ecex_results.mat'), 'ecex_results');

acex_results = procrustes_decoding('ac','ex', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'acex_results.mat'), 'acex_results');

ecac_results = procrustes_decoding('ec','ac', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'ecac_results.mat'), 'ecac_results');

exec_results = procrustes_decoding('ex','ec', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'exec_results.mat'), 'exec_results');

exac_results = procrustes_decoding('ex','ac', data_trial, labels, neuron_num_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path, 'exac_results.mat'), 'exac_results');
toc;
%%
function results = procrustes_decoding(stim1,stim2,data_trial,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat)
% stim1 is transformed(comparison); stim2 is the training data(target)

pair_wise_data_trial = pair_pcaloader(data_trial,stim1,stim2);

%load data
stim1_data = pair_wise_data_trial(1:500,:);
stim2_data = pair_wise_data_trial(501:1000,:);

% zscore
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);

results = cell(1,length(neuron_num_list)); % each cell stores results for one neuron-number decoding

for neuron_squence = 1
    disp([stim1,stim2,num2str(neuron_squence)])
    result_one_cell = [];
    parfor neuron_repeat = 1:neuron_sample_repeat
        % sample different neuron, for neuron_sample_repeat times
        sampled_neuron_idx = sample_neuron(neuron_num_list,neuron_squence);
        stim1_data_sample = stim1_data(:,sampled_neuron_idx);
        stim2_data_sample = stim2_data(:,sampled_neuron_idx);
        % sample trials multiple times
        trial_result = zeros(trial_sample_repeat,6);
        for trial_repeat =1:trial_sample_repeat
            [training_data,training_label,stim1_test_data,target_label, ...
             true_rot_transformed,shuffle_transformed,shuffle_identity_label] = ...
                data_trial_sampler(stim1_data_sample,stim2_data_sample,labels);

            [genAcc,acc_true_rot,acc_shuffle_target,acc_shuffle_identity,acc_notransform,acc_chance] = ...
                pro_decoding(training_data,training_label,stim1_test_data,target_label, ...
                             true_rot_transformed,shuffle_transformed,shuffle_identity_label);

            trial_result(trial_repeat,:) = [genAcc,acc_true_rot,acc_shuffle_target,acc_shuffle_identity,acc_notransform,acc_chance];
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

function [training_data,training_label,stim1_test_data,target_label, ...
          true_rot_transformed,shuffle_transformed,shuffle_identity_label] = ...
          data_trial_sampler(stim1_data_sample,stim2_data_sample,labels)
% sample trial
tial_out_of_ten = 1;
valid_trial_idx = sample_trial(tial_out_of_ten); % 10 fold
training_trial_idx = setdiff([1:500],valid_trial_idx);

% SVM training data, all stim2 data
training_data = stim2_data_sample;
training_label = labels;

% test data (one held-out trial per condition, ordered condition 1..50)
stim1_test_data = stim1_data_sample(valid_trial_idx,:);
target_label     = labels(valid_trial_idx,:);   % == (1:50)'

% averaged data used to learn the transform
stim2_target_data = stim2_data_sample;
stim1_2transform_data = stim1_data_sample(training_trial_idx,:);
stim2_target_data_averged     = take_average(stim2_target_data,10);
stim1_2transform_data_averged = take_average(stim1_2transform_data,10 - tial_out_of_ten);

% ---- TRUE correspondence transform (rotation only, matched baseline) ----
[d,Z,transform] = procrustes(stim2_target_data_averged,stim1_2transform_data_averged);
t_matrix = transform.T;
true_rot_transformed = stim1_test_data*t_matrix;   % rotation only (no b, no c)

% ---- Structure-destroying (correspondence-shuffled) control ----
% One permutation pi is used to shuffle BOTH the mean matrix and the held-out
% test trials, so that within the shuffled ordering the trial<->mean
% correspondence is preserved, while the true stim1<->stim2 condition
% correspondence is destroyed.
shuffle_idx = randperm(size(stim1_2transform_data_averged,1));

% shuffle the mean matrix and learn the (broken-correspondence) rotation
toy_averaged = stim1_2transform_data_averged(shuffle_idx,:);
[d_fake,Z_fake,t_fake] = procrustes(stim2_target_data_averged,toy_averaged);
t_fakematrix = t_fake.T;

% shuffle the held-out test trials the same way, then apply rotation only
stim1_test_data_shuffled = stim1_test_data(shuffle_idx,:);
shuffle_transformed = stim1_test_data_shuffled*t_fakematrix;   % rotation only (no toy_b, no toy_c)

% labels for the shuffled test trials:
%   target_label (1:50) is the target position each shuffled trial was aligned to.
%   shuffle_identity_label (pi) is the trial's true condition identity BEFORE rotation.
shuffle_identity_label = target_label(shuffle_idx);

% To use the FULL transform instead of rotation-only, replace the two
% "rotation only" lines above with:
%   true_rot_transformed = transform.b*stim1_test_data*t_matrix + transform.c;
%   shuffle_transformed  = t_fake.b*stim1_test_data_shuffled*t_fakematrix + t_fake.c;
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

function [genAcc,acc_true_rot,acc_shuffle_target,acc_shuffle_identity,acc_notransform,acc_chance] = ...
    pro_decoding(training_data,training_label,stim1_test_data,target_label, ...
                 true_rot_transformed,shuffle_transformed,shuffle_identity_label)

stim2_model = fitcecoc(training_data,training_label);

% stim2 within-condition generalization (ceiling)
CVMdl1 = crossval(stim2_model,'KFold',10);
genAcc = 1 - kfoldLoss(CVMdl1);

% no-transform transfer decoding (raw stim1 test into stim2 model)
predicted_notransform = predict(stim2_model,stim1_test_data);
acc_notransform = sum(target_label == predicted_notransform)/length(predicted_notransform);

% TRUE correspondence, rotation-only transfer decoding (matched baseline)
predicted_true = predict(stim2_model,true_rot_transformed);
acc_true_rot = sum(target_label == predicted_true)/length(predicted_true);

% Structure-shuffled control, rotation-only, evaluated against BOTH label sets
predicted_shuffle = predict(stim2_model,shuffle_transformed);
%   (A) target labels (1:50): can PT map the shuffled config onto the target config?
acc_shuffle_target   = sum(target_label == predicted_shuffle)/length(predicted_shuffle);
%   (B) shuffled identity labels (pi): is the trial's true pre-rotation identity preserved?
acc_shuffle_identity = sum(shuffle_identity_label == predicted_shuffle)/length(predicted_shuffle);

% chance baseline
acc_chance = sum(target_label == randperm(50)')/length(target_label);
end
