%{
Filename: procrustes_decoding_per_session.m
Author: Zitong Wang
Date: 2026-04-28

Description:
    Per-session variant of procrustes_decoding_trial_perturbed.m. Instead
    of decoding from the full pseudo-population, this script loops over
    individual recording sessions, takes only the simultaneously-recorded
    neurons from one session at a time, and runs the full Procrustes
    transfer-decoding pipeline on each. Trial perturbation (none /
    shuffle / affine) is applied per-session before decoding so the same
    pipeline can be used to test the contribution of *real* within-
    session noise correlations.

    No neuron downsampling is performed - each session uses all of its
    good channels (typically 5-15 neurons).

    The pseudo-population data file is reused as the source: columns of
    data_trial{i} are sliced by perturb_opts.neuron_session_id (from
    build_affine_shift_config). Within a session, the trial dimension
    reflects genuine simultaneous recording, so per-session noise
    correlations are real.

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat
    - per-session .mat files in neuronal_data/<monkey>/<area>/ac/
      (used by build_affine_shift_config for session membership)
    - timewindow:    [start_ms end_ms] spike-count window
    - perturb_mode:  'none' | 'shuffle' | 'affine'

Outputs:
    MAT-files saved to:
      results/decoding_outputs/procrustes_decoding_per_session_results/<perturb_mode>/<monkey>/<area>/
      Named: acec_results.mat, ecex_results.mat, acex_results.mat,
             ecac_results.mat, exec_results.mat, exac_results.mat

    Each file contains two variables:
      <pair>_results : struct array, length N_sessions, with fields
            .session_idx        (1..N_sessions)
            .session_file       (e.g. 'FR081414_2.mat')
            .n_neurons          (good channels in this session)
            .neuron_pattern_id  (affine pattern indices used; relevant
                                 for affine mode, recorded for all)
            .accuracy           [trial_sample_repeat x 8]
            .accuracy_columns   {1 x 8} cell array naming the columns
      config : struct holding the run parameters (monkey, vp, timewindow,
               perturb_mode, perturb_seed, trial_sample_repeat,
               tial_out_of_ten, N_sessions, session_files, accuracy_columns)
%}

%%
clc; clear;

%% Configuration
monkey = 'KO';      % 'FR' or 'KO'
vp     = 'V1';      % 'V1' or 'V2'

%% Trial perturbation configuration
perturb_mode = 'none';   % 'none' | 'shuffle' | 'affine'
perturb_seed = 1;          % RNG seed applied before perturbation (for 'shuffle')

%% Decoding parameters
trial_sample_repeat = 50;
timewindow          = [330 630];   % e.g., early: [340 410], late: [410 480], long-late: [410 550]

%% Accuracy column names (matches pro_decoding output order)
accuracy_columns = {'genAcc', ...
                    'accuracy_s1', ...
                    'accuracy_transformed_s1', ...
                    'accuracy_rand', ...
                    'accuracy_onlyscale', ...
                    'accuracy_onlyrotation', ...
                    'accuracy_onlytranslation', ...
                    'accuracy_non_transfer_control'};

%% Load neuronal data for all stimuli
data_file = fullfile('..','..','neuronal_data', monkey, vp, ...
    sprintf('%s_%s_allstim.mat', monkey, vp));
tmp = load(data_file, 'three_stim_array');
spike_data = tmp.three_stim_array;

%% Prepare output directory (tagged by perturb_mode)
save_path = fullfile('..','..','results','decoding_outputs', ...
    'procrustes_decoding_per_session_results', perturb_mode, monkey, vp);
if ~exist(save_path, 'dir'), mkdir(save_path); end

fprintf('Monkey %s | %s | perturb_mode = %s | per-session\n', monkey, vp, perturb_mode);

%% Load trial-by-trial data and labels
num_conditions = numel(spike_data);
data_trial  = cell(1, num_conditions);
label_trial = cell(1, num_conditions);
for i = 1:num_conditions
    [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, timewindow);
end

%% Build session / affine config (always, since session_id is needed for slicing)
N_neurons    = size(data_trial{1}, 2);
perturb_opts = build_affine_shift_config(monkey, vp, N_neurons);
session_files          = perturb_opts.session_files;
neuron_session_id_full = perturb_opts.neuron_session_id;
neuron_pattern_id_full = perturb_opts.neuron_pattern_id;
N_sessions             = numel(session_files);

fprintf('Found %d sessions, %d total neurons in the pseudo-population.\n', ...
    N_sessions, N_neurons);

%% Initialize per-pair result struct arrays
pair_names = {'acec','ecex','acex','ecac','exec','exac'};
pair_stim1 = {'ac' ,'ec' ,'ac' ,'ec' ,'ex' ,'ex' };
pair_stim2 = {'ec' ,'ex' ,'ex' ,'ac' ,'ec' ,'ac' };

empty_entry = struct( ...
    'session_idx',       [], ...
    'session_file',      '', ...
    'n_neurons',         [], ...
    'neuron_pattern_id', [], ...
    'accuracy',          [], ...
    'accuracy_columns',  {accuracy_columns});

pair_results = struct();
for p = 1:numel(pair_names)
    pair_results.(pair_names{p}) = repmat(empty_entry, N_sessions, 1);
end

%% Per-session decoding loop
rng(perturb_seed);   % reproducible shuffles; no-op for deterministic affine / none
labels = label_trial{1};

tic;
for s = 1:N_sessions
    session_mask        = (neuron_session_id_full == s);
    n_session_neurons   = sum(session_mask);
    session_pattern_ids = neuron_pattern_id_full(session_mask);

    fprintf('\n--- Session %d/%d (%s): %d neurons ---\n', ...
        s, N_sessions, session_files{s}, n_session_neurons);

    % Slice each rendering's data to this session's columns
    session_data_trial = cell(1, num_conditions);
    for i = 1:num_conditions
        session_data_trial{i} = data_trial{i}(:, session_mask);
    end

    % Apply per-session perturbation
    if ~strcmp(perturb_mode, 'none')
        if strcmp(perturb_mode, 'affine')
            session_perturb_opts = struct( ...
                'affine_patterns',   perturb_opts.affine_patterns, ...
                'neuron_pattern_id', session_pattern_ids);
        else
            session_perturb_opts = struct();
        end
        for i = 1:num_conditions
            session_data_trial{i} = apply_trial_perturbation( ...
                session_data_trial{i}, label_trial{i}, ...
                perturb_mode, session_perturb_opts);
        end
        fprintf('  Applied %s perturbation.\n', perturb_mode);
    end

    % Run all 6 stimulus pairs for this session
    for p = 1:numel(pair_names)
        accuracy_matrix = procrustes_decoding( ...
            pair_stim1{p}, pair_stim2{p}, ...
            session_data_trial, labels, trial_sample_repeat);

        pair_results.(pair_names{p})(s).session_idx       = s;
        pair_results.(pair_names{p})(s).session_file      = session_files{s};
        pair_results.(pair_names{p})(s).n_neurons         = n_session_neurons;
        pair_results.(pair_names{p})(s).neuron_pattern_id = session_pattern_ids;
        pair_results.(pair_names{p})(s).accuracy          = accuracy_matrix;
        pair_results.(pair_names{p})(s).accuracy_columns  = accuracy_columns;
    end
end
toc;

%% Build config struct
config = struct( ...
    'monkey',              monkey, ...
    'vp',                  vp, ...
    'timewindow',          timewindow, ...
    'perturb_mode',        perturb_mode, ...
    'perturb_seed',        perturb_seed, ...
    'trial_sample_repeat', trial_sample_repeat, ...
    'tial_out_of_ten',     1, ...
    'N_sessions',          N_sessions, ...
    'session_files',       {session_files}, ...
    'accuracy_columns',    {accuracy_columns});

%% Save one file per stimulus pair
for p = 1:numel(pair_names)
    pair_name = pair_names{p};
    save_struct = struct();
    save_struct.(sprintf('%s_results', pair_name)) = pair_results.(pair_name);
    save_struct.config = config;
    save(fullfile(save_path, sprintf('%s_results.mat', pair_name)), ...
        '-struct', 'save_struct');
end
fprintf('\nSaved 6 per-pair files to: %s\n', save_path);

%% ==================== LOCAL FUNCTIONS ====================
function results = procrustes_decoding(stim1, stim2, data_trial, labels, trial_sample_repeat)
% Per-session decoding: uses ALL neurons in data_trial, no neuron downsampling.
% stim1 is transformed (comparison); stim2 is the training data (target).
% Returns: [trial_sample_repeat x 8] accuracy matrix.

pair_wise_data_trial = pair_pcaloader(data_trial, stim1, stim2);

% load data
stim1_data = pair_wise_data_trial(1:500,   :);
stim2_data = pair_wise_data_trial(501:1000,:);

% zscore
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);

disp([stim1, stim2]);

results = zeros(trial_sample_repeat, 8);
parfor trial_repeat = 1:trial_sample_repeat
    [training_data, training_label, ...
     stim1_test_data, stim1_test_label, ...
     transformed_stim1, toy_transd, partially_transformed_stim1] = ...
        data_trial_sampler(stim1_data, stim2_data, labels);

    [genAcc, accuracy_s1, accuracy_transformed_s1, accuracy_rand, ...
     accuracy_onlyscale, accuracy_onlyrotation, accuracy_onlytranslation, ...
     accuracy_non_transfer_control] = ...
        pro_decoding(training_data, training_label, ...
                     stim1_test_data, stim1_test_label, ...
                     transformed_stim1, toy_transd, partially_transformed_stim1);

    results(trial_repeat, :) = [genAcc, accuracy_s1, accuracy_transformed_s1, ...
                                accuracy_rand, accuracy_onlyscale, ...
                                accuracy_onlyrotation, accuracy_onlytranslation, ...
                                accuracy_non_transfer_control];
end
end

function stim_data_trial_averged = take_average(stim_data, number_of_average)
[trial, neuron] = size(stim_data); %#ok<ASGLU>
stim_data_trial_averged = zeros(50, neuron);
for i = 1:50
    temp = stim_data(i*number_of_average-(number_of_average-1):i*number_of_average, :);
    temp_mean = mean(temp, 1);
    stim_data_trial_averged(i, :) = temp_mean;
end
end

function [training_data, training_label, stim1_test_data, stim1_test_label, ...
          transformed_stim1, toy_transd, partially_transformed_stim1] = ...
    data_trial_sampler(stim1_data_sample, stim2_data_sample, labels)
% sample trial
tial_out_of_ten = 1;
valid_trial_idx     = sample_trial(tial_out_of_ten); % 10 fold
training_trial_idx  = setdiff(1:500, valid_trial_idx);

% SVM training data, all stim2 data
training_data  = stim2_data_sample;
training_label = labels;

% test data
stim1_test_data  = stim1_data_sample(valid_trial_idx, :);
stim1_test_label = labels(valid_trial_idx, :);

% transformed data
stim2_target_data       = stim2_data_sample;
stim1_2transform_data   = stim1_data_sample(training_trial_idx, :);
% take the mean
stim2_target_data_averged     = take_average(stim2_target_data, 10);
stim1_2transform_data_averged = take_average(stim1_2transform_data, 10 - tial_out_of_ten);

[d, Z, transform] = procrustes(stim2_target_data_averged, stim1_2transform_data_averged); %#ok<ASGLU>
t_matrix = transform.T;
b = transform.b;
c = transform.c;
transformed_stim1 = b*stim1_test_data*t_matrix + c;

% partially transformed data - only 1 transform
scaling_only_transformed_stim1     = b*stim1_test_data;
rotation_only_transformed_stim1    = stim1_test_data*t_matrix;
translation_only_transformed_stim1 = stim1_test_data + c;
partially_transformed_stim1 = {scaling_only_transformed_stim1, ...
                               rotation_only_transformed_stim1, ...
                               translation_only_transformed_stim1};

% randomized data
shuffle_idx  = randperm(size(stim1_2transform_data_averged, 1));
toy_averaged = stim1_2transform_data_averged(shuffle_idx, :);
[d_fake, Z_fake, t_fake] = procrustes(stim2_target_data_averged, toy_averaged); %#ok<ASGLU>
t_fakematrix = t_fake.T;
toy_b        = t_fake.b;
toy_c        = t_fake.c;
toy_transd   = toy_b * stim1_test_data * t_fakematrix + toy_c;
end

function sample_trial_idx = sample_trial(trialnum_outof_ten)
% trial sampling
sample_trial_matrix = zeros(50, 10);
for i = 1:50
    sample_trial_matrix(i, :) = randperm(10);
end
sample_trial_4condition = sample_trial_matrix(:, 1:trialnum_outof_ten);
base = (0:10:490)';
sample_trial_idx = sample_trial_4condition + base;
sample_trial_idx = reshape(sample_trial_idx, [], 1);
sample_trial_idx = sort(sample_trial_idx);
end

function pair_wise_pca_data = pair_pcaloader(dca_data, stim1, stim2)
% use the input stimulus name to load pair_wise data
name2idx = struct('ac', 1, 'ec', 2, 'ex', 3);
stim1_data = dca_data{name2idx.(stim1)};
stim2_data = dca_data{name2idx.(stim2)};
pair_wise_pca_data = [stim1_data; stim2_data];
end

function [genAcc, accuracy_s1, accuracy_transformed_s1, accuracy_rand, ...
          accuracy_onlyscale, accuracy_onlyrotation, accuracy_onlytranslation, ...
          accuracy_non_transfer_control] = ...
    pro_decoding(training_data, training_label, stim1_test_data, stim1_test_label, ...
                 transformed_stim1, toy_transd, partially_transformed_stim1)
% get accs
stim2_model = fitcecoc(training_data, training_label);

CVMdl1    = crossval(stim2_model, 'KFold', 10);
genError1 = kfoldLoss(CVMdl1);
genAcc    = 1 - genError1;

% non-transformation transfer decoding
predicted_Label_s1 = predict(stim2_model, stim1_test_data);
accuracy_s1 = sum(stim1_test_label == predicted_Label_s1) / length(predicted_Label_s1);
% with-transformation transfer decoding
predicted_Label_transformed_s1 = predict(stim2_model, transformed_stim1);
accuracy_transformed_s1 = sum(stim1_test_label == predicted_Label_transformed_s1) / length(predicted_Label_transformed_s1);
% random-transformation transfer decoding
predicted_Label_rand = predict(stim2_model, toy_transd);
accuracy_rand = sum(stim1_test_label == predicted_Label_rand) / length(predicted_Label_rand);
% partial-transformation transfer decoding
predicted_Label_scale       = predict(stim2_model, partially_transformed_stim1{1,1});
accuracy_onlyscale          = sum(stim1_test_label == predicted_Label_scale)       / length(predicted_Label_scale);
predicted_Label_rotation    = predict(stim2_model, partially_transformed_stim1{1,2});
accuracy_onlyrotation       = sum(stim1_test_label == predicted_Label_rotation)    / length(predicted_Label_rotation);
predicted_Label_translation = predict(stim2_model, partially_transformed_stim1{1,3});
accuracy_onlytranslation    = sum(stim1_test_label == predicted_Label_translation) / length(predicted_Label_translation);
% non-transformation transfer decoding control
accuracy_non_transfer_control = sum(randperm(50)' == predicted_Label_s1) / length(predicted_Label_s1);
end
