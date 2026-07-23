%{
Filename: procrustes_decoding_basic_incremental_pr_trial_perturbed.m
Author: Zitong Wang
Date: 2026-07-23

Location:
    This script lives in code/noise_correlation_scripts/ so that it can call the
    standalone helpers apply_trial_perturbation.m and build_affine_shift_config.m
    without duplicating them. Relative data/results paths ('..','..',...) are the
    same depth as the other scripts in that folder, so they are unchanged.

Description:
    Incremental Procrustes-based transfer decoding with participation ratio (PR),
    run on trial-perturbed (noise-correlation-broken) pseudopopulation data.

    This is procrustes_decoding_basic_incremental_pr.m with two additions:

      (A) TRIAL PERTURBATION (from procrustes_decoding_trial_perturbed.m).
          A per-neuron, per-stimulus trial permutation is applied immediately
          after multiclass_svmloader_PT and before any downstream processing,
          so the decoding pipeline itself is unmodified. Modes:
            'none'    - pass-through baseline (identical to the unperturbed PR script)
            'shuffle' - independent random permutation of the 10 trials per neuron
                        per stimulus (matches noise_correlation_analysis.m)
            'affine'  - deterministic affine trial shift t -> (a*p + b) mod 10,
                        assigned per-neuron by session membership so that
                        simultaneously recorded neurons receive distinct shifts
                        (matches noise_correlation_analysis_affine.m)
          The marginal firing-rate distribution per neuron per stimulus is
          preserved exactly; only trial-aligned co-variability is destroyed.

          NOTE ON PR: because the perturbation only reorders trials within a
          (neuron, stimulus) block, the trial-averaged 50 x N signal manifold is
          numerically IDENTICAL to the unperturbed one. PR is therefore expected
          to be unchanged by the perturbation (up to neuron-sampling RNG), and
          serves here as a built-in sanity check rather than a second measurement.
          Any PR difference vs the unperturbed run at matched neuron counts
          indicates a bug, not a noise-correlation effect.

      (B) MULTI-SETTING BATCH. All four monkey x area combinations (FR/KO x V1/V2)
          can run in one invocation, while a single setting is still available by
          narrowing monkey_list / vp_list in the Configuration block. Each setting
          is self-contained: data are reloaded, the affine config is rebuilt for
          that monkey/area, and the RNG is reset the same way a standalone run
          would, so batch results are identical to running each setting alone.

    Everything else - the 10-unit neuron sweep, the PR computation, the struct
    output format, the correspondence-shuffled control, and all sampling - is
    carried over unchanged from procrustes_decoding_basic_incremental_pr.m.

    PR is computed on the z-scored, trial-averaged 50 x N response matrix (the
    signal manifold), centered across stimuli, as (sum lambda)^2 / sum(lambda^2)
    over the covariance spectrum. It is computed once per neuron-sampling repeat
    (independent of the trial split) and is invariant to the Procrustes transform
    (rotation / uniform scale / translation), so it characterizes the same
    geometry the decoder and PT operate on. No second renormalization is applied
    after trial-averaging: each unit is weighted by its stimulus-driven (reliable)
    variance, matching the decoding geometry.

    Column 4 of .acc does NOT hold the base script's random-transformation
    control: it has been replaced by a self-consistent, correspondence-shuffled
    control. One permutation pi shuffles BOTH the trial-averaged mean matrix and
    the held-out test trials; a rotation-only Procrustes transform is learned from
    the shuffled means and applied to the identically-shuffled test trials.
    Decoding against target labels 1:50 then measures whether PT can force a
    mapping even when the true stim1<->stim2 condition correspondence is destroyed
    (the shuffle-target accuracy). The shuffle-identity variant is not computed
    here. The matched true-correspondence baseline is the rotation-only column
    (col 6, accuracy_onlyrotation).

Output (per pair file): a 1 x numel(neuron_list) cell array; each cell is a
struct with fields:
    .neuron_num     scalar, number of sampled units
    .stim1, .stim2  condition labels ('ac'/'ec'/'ex'); stim1 transformed, stim2 trained
    .acc            [(neuron_sample_repeat*trial_sample_repeat) x 8] decoding accuracies.
                    Column order:
                      1 genAcc                        stim2 within-condition 10-fold CV accuracy (ceiling)
                      2 accuracy_s1                    raw stim1 test, no transform (vs 1:50)
                      3 accuracy_transformed_s1        TRUE corr., FULL transform b*X*T+c (vs 1:50)
                      4 accuracy_rand                  CORRESPONDENCE-SHUFFLED control, rotation-only,
                                                       decoded vs TARGET labels 1:50 (shuffle-target acc)
                      5 accuracy_onlyscale             TRUE corr., scaling-only b*X (vs 1:50)
                      6 accuracy_onlyrotation          TRUE corr., rotation-only X*T (vs 1:50) -- matched
                                                       baseline for column 4
                      7 accuracy_onlytranslation       TRUE corr., translation-only X+c (vs 1:50)
                      8 accuracy_non_transfer_control  random-permutation chance baseline
    .acc_repeat_id  [(...) x 1] neuron-repeat index for each acc row (for PR<->acc pairing)
    .pr_stim1       [neuron_sample_repeat x 1] PR of the stim1 manifold per repeat
    .pr_stim2       [neuron_sample_repeat x 1] PR of the stim2 manifold per repeat

    To pool a condition's PR across the analysis, gather pr_stim1 from the pairs
    where stim1 == that condition and pr_stim2 from the pairs where stim2 == that
    condition (each condition appears in 4 of the 6 ordered pairs).

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat
        Contains 'three_stim_array': cell array of [trial x neuron] matrices.
    - timewindow: two-element vector [start_ms end_ms] for spike counts.
    - perturb_mode: 'none' | 'shuffle' | 'affine'

Dependencies (same folder):
    multiclass_svmloader_PT.m, apply_trial_perturbation.m, build_affine_shift_config.m

Outputs:
    - MAT-files saved to:
      results/decoding_outputs/procrustes_decoding_basic_incremental_pr_trial_perturbed_results/<perturb_mode>/<monkey>/<area>/
      Files: acec_results.mat, ecex_results.mat, acex_results.mat,
             ecac_results.mat, exec_results.mat, exac_results.mat
%}

%% Initialize MATLAB environment
clc; clear;

%% Configuration -- which monkey / area settings to run
% Run all four combinations in one invocation:
monkey_list = {'FR','KO'};   % subset of {'FR','KO'}
vp_list     = {'V1','V2'};   % subset of {'V1','V2'}
%
% To run a SINGLE setting, narrow the lists, e.g.:
%   monkey_list = {'FR'};  vp_list = {'V2'};

%% Trial perturbation configuration
perturb_mode = 'affine';   % 'none' | 'shuffle' | 'affine'
perturb_seed = 1;          % RNG seed applied before perturbation (for 'shuffle')

%% Shared analysis parameters
cfg = struct();
cfg.timewindow           = [330 630];  % e.g., early: [340 410], late: [410 480]
cfg.perturb_mode         = perturb_mode;
cfg.perturb_seed         = perturb_seed;
cfg.neuron_sample_repeat = 20;         % fewer repeats for incremental test
cfg.trial_sample_repeat  = 15;
cfg.decode_seed          = 1;          % RNG seed applied before decoding

% If one setting errors (e.g. a missing session file for the affine config),
% keep going with the remaining settings instead of aborting the whole batch.
continue_on_error = true;

%% Build the run list (monkey x area)
run_list = cell(0,2);
for m = 1:numel(monkey_list)
    for v = 1:numel(vp_list)
        run_list(end+1,:) = {monkey_list{m}, vp_list{v}}; %#ok<SAGROW>
    end
end
n_runs = size(run_list,1);

fprintf('=== Incremental PR decoding | perturb_mode = %s | %d setting(s) ===\n', ...
    perturb_mode, n_runs);

%% Run every setting
run_status = cell(n_runs,1);
batch_timer = tic;
for r = 1:n_runs
    monkey = run_list{r,1};
    vp     = run_list{r,2};
    fprintf('\n--- [%d/%d] %s %s ---\n', r, n_runs, monkey, vp);

    if continue_on_error
        try
            run_one_setting(monkey, vp, cfg);
            run_status{r} = 'ok';
        catch ME
            run_status{r} = sprintf('FAILED (%s)', ME.message);
            fprintf(2, 'Setting %s %s failed: %s\n', monkey, vp, ME.message);
        end
    else
        run_one_setting(monkey, vp, cfg);
        run_status{r} = 'ok';
    end
end

%% Batch summary
fprintf('\n=== Batch summary (total %.1f min) ===\n', toc(batch_timer)/60);
for r = 1:n_runs
    fprintf('  %s %s : %s\n', run_list{r,1}, run_list{r,2}, run_status{r});
end

%% ==================== PER-SETTING DRIVER ====================
function run_one_setting(monkey, vp, cfg)
% Runs the full 6-pair incremental PR decoding for one monkey/area, with the
% configured trial perturbation applied before any downstream processing.

%% Load neuronal data for all stimuli
data_file = fullfile('..','..','neuronal_data', monkey, vp, ...
    sprintf('%s_%s_allstim.mat', monkey, vp));
tmp = load(data_file, 'three_stim_array');
spike_data = tmp.three_stim_array;

%% Prepare output directory (tagged by perturb_mode)
save_path = fullfile('..','..','results','decoding_outputs', ...
    'procrustes_decoding_basic_incremental_pr_trial_perturbed_results', ...
    cfg.perturb_mode, monkey, vp);
if ~exist(save_path, 'dir'), mkdir(save_path); end

fprintf('Monkey %s | %s | perturb_mode = %s\n', monkey, vp, cfg.perturb_mode);

%% Load trial-by-trial data
nconds = numel(spike_data);
data_trial  = cell(1,nconds);
label_trial = cell(1,nconds);
for i = 1:nconds
    [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, cfg.timewindow);
end

%% Build perturbation opts (affine only)
perturb_opts = struct();
if strcmp(cfg.perturb_mode, 'affine')
    N_neurons = size(data_trial{1}, 2);
    perturb_opts = build_affine_shift_config(monkey, vp, N_neurons);
    fprintf('Affine config: %d neurons across %d sessions.\n', ...
        N_neurons, numel(perturb_opts.session_files));
end

%% Apply trial perturbation independently to each rendering
if ~strcmp(cfg.perturb_mode, 'none')
    rng(cfg.perturb_seed);   % reproducible shuffles; no-op for deterministic affine
    for i = 1:nconds
        data_trial{i} = apply_trial_perturbation( ...
            data_trial{i}, label_trial{i}, cfg.perturb_mode, perturb_opts);
    end
    fprintf('Trial perturbation (%s) applied to all %d conditions.\n', ...
        cfg.perturb_mode, nconds);
end

%% Build incremental neuron-count list (uniform 10-unit steps)
% Base code maps monkey/area to maximum neurons (some areas list two reference
% sizes, e.g. a V1-matched size and the full population).
neuron_code = {[48]; [48 112]; [109]; [109 146]};
map = struct('FRV1',1, 'FRV2',2, 'KOV1',3, 'KOV2',4);
idx = map.(strcat(monkey, vp));
max_list = neuron_code{idx};

% Uniform 10-step grid up to the maximum, with the reference sizes forced in.
neuron_list = unique([10:10:max(max_list), max_list]);
labels = label_trial{1};

%% Decoding parameters
neuron_sample_repeat = cfg.neuron_sample_repeat;
trial_sample_repeat  = cfg.trial_sample_repeat;
rng(cfg.decode_seed);

tic;
% Run decoding for each transformation pair
acec_results = incre_decoding('ac','ec', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path,'acec_results.mat'),'acec_results');

ecex_results = incre_decoding('ec','ex', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path,'ecex_results.mat'),'ecex_results');

acex_results = incre_decoding('ac','ex', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path,'acex_results.mat'),'acex_results');

ecac_results = incre_decoding('ec','ac', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path,'ecac_results.mat'),'ecac_results');

exec_results = incre_decoding('ex','ec', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path,'exec_results.mat'),'exec_results');

exac_results = incre_decoding('ex','ac', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat);
save(fullfile(save_path,'exac_results.mat'),'exac_results');
toc;
end

%% ==================== DECODING (unchanged from the PR script) ====================
function results = incre_decoding(stim1,stim2,data_trial,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat)
% stim1 is transformed, stim2 is the training data

pair_wise_data_trial = pair_pcaloader(data_trial,stim1,stim2);

%load data
stim1_data = pair_wise_data_trial(1:500,:);
stim2_data = pair_wise_data_trial(501:1000,:);

% zscore (per unit, over the full set of trials)
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);

results = cell(1,length(neuron_num_list)); % each cell stores results for one neuron-number decoding

for neuron_squence = 1:length(neuron_num_list)
    disp([stim1,stim2,num2str(neuron_squence)])

    % sliced outputs for parfor: one acc block / one id block / one PR value per neuron-repeat
    acc_blocks = cell(neuron_sample_repeat,1);
    id_blocks  = cell(neuron_sample_repeat,1);
    pr1_vec    = zeros(neuron_sample_repeat,1);
    pr2_vec    = zeros(neuron_sample_repeat,1);

    parfor neuron_repeat = 1:neuron_sample_repeat
        % sample different neuron subset, for neuron_sample_repeat times
        sampled_neuron_idx = sample_neuron(neuron_num_list,neuron_squence);
        stim1_data_sample = stim1_data(:,sampled_neuron_idx);
        stim2_data_sample = stim2_data(:,sampled_neuron_idx);

        % participation ratio of each condition's signal manifold on this
        % subset (independent of the trial split -> computed once per repeat)
        pr1_vec(neuron_repeat) = compute_pr(stim1_data_sample);
        pr2_vec(neuron_repeat) = compute_pr(stim2_data_sample);

        % sample trials multiple times
        trial_result = zeros(trial_sample_repeat,8);
        for trial_repeat =1:trial_sample_repeat
            [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,shuffle_transformed,partially_transformed_stim1] = data_trial_sampler(stim1_data_sample,stim2_data_sample,labels);
            [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control] = pro_decoding(training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,shuffle_transformed,partially_transformed_stim1);

            trial_result(trial_repeat,:) = [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_rand,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control];

        end
        acc_blocks{neuron_repeat} = trial_result;
        id_blocks{neuron_repeat}  = repmat(neuron_repeat,trial_sample_repeat,1);
    end

    one = struct();
    one.neuron_num    = neuron_num_list(neuron_squence);
    one.stim1         = stim1;
    one.stim2         = stim2;
    one.acc           = cat(1, acc_blocks{:});   % (nrep*trep) x 8
    one.acc_repeat_id = cat(1, id_blocks{:});    % (nrep*trep) x 1
    one.pr_stim1      = pr1_vec;                 % nrep x 1
    one.pr_stim2      = pr2_vec;                 % nrep x 1
    results{neuron_squence} = one;
end

end

function pr = compute_pr(zscored_data)
% Participation ratio of the trial-averaged signal manifold.
% Input : zscored_data, a [500 x N] matrix (50 stimuli x 10 trials),
%         already z-scored per unit over the full set of trials.
% Method: trial-average to [50 x N], center across stimuli, then
%         PR = (sum lambda)^2 / sum(lambda^2) over the covariance spectrum.
%         Computed via SVD of the centered trial-averaged matrix, which is
%         scale-invariant and avoids rank-deficiency artifacts when N > 50.
avg = take_average(zscored_data, 10);   % 50 x N signal manifold
avg = avg - mean(avg, 1);               % center across the 50 stimuli
s = svd(avg, 'econ');                   % singular values
lambda = s.^2;                          % covariance eigenvalues (up to a constant)
lambda = lambda(lambda > 0);            % numerical safety
pr = (sum(lambda)^2) / sum(lambda.^2);
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

function [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,shuffle_transformed,partially_transformed_stim1] = data_trial_sampler(stim1_data_sample,stim2_data_sample,labels)
% sample trial
tial_out_of_ten = 1;
valid_trial_idx = sample_trial(tial_out_of_ten); % 10 fold
training_trial_idx = setdiff([1:500],valid_trial_idx);

% SVM training data, all stim2 data
training_data = stim2_data_sample;
training_label = labels;

% test data (one held-out trial per condition, ordered condition 1..50)
stim1_test_data = stim1_data_sample(valid_trial_idx,:);
stim1_test_label = labels(valid_trial_idx,:);

% transformed data
stim2_target_data = stim2_data_sample;
stim1_2transform_data = stim1_data_sample(training_trial_idx,:);
% take the mean
stim2_target_data_averged = take_average(stim2_target_data,10);
stim1_2transform_data_averged = take_average(stim1_2transform_data,10 - tial_out_of_ten);

% ---- TRUE correspondence transform (FULL transform) ----
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


% To use the FULL transform for the shuffle control instead of rotation-only:
%   shuffle_transformed = t_fake.b*stim1_test_data_shuffled*t_fakematrix + t_fake.c;


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
    = pro_decoding(training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,shuffle_transformed,partially_transformed_stim1)
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
% correspondence-shuffled control transfer decoding (vs target labels 1:50)
predicted_Label_rand = predict(stim2_model,shuffle_transformed);
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
