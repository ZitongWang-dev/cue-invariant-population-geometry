%{
Filename: procrustes_decoding_basic_incremental_pc.m
Author: Zitong Wang (PC-dimension variant)
Date: 2026-07-27

Description:
    Procrustes-based transfer decoding as a function of PRINCIPAL-COMPONENT
    dimension. Controlled companion to procrustes_decoding_basic_incremental_pr.m:
    instead of sweeping the number of sampled NEURONS (and measuring PR at each),
    it fixes the full recorded population and sweeps the number of PC DIMENSIONS
    retained, so dimension is controlled directly rather than measured.

    Structure is kept parallel to the pr script on purpose:
        incre_decoding  -> data_trial_sampler -> pro_decoding
    with take_average / sample_trial / pair_pcaloader unchanged. The differences
    are localized to (a) a per-cue mean-PCA basis + projection inside
    data_trial_sampler, (b) the new shuffle control, (c) the pdist output, and
    (d) an outer loop over all four monkey/area settings.

    Per ordered pair (stim1 transformed -> stim2 trained), per trial split, per
    PC-count k, data_trial_sampler:
      (1) z-score is already applied per cue over its own 500 trials (in
          incre_decoding, unchanged from the pr script);
      (2) builds a SEPARATE per-cue PCA basis from the 50 x N centered condition
          MEANS (the signal manifold PT/PR live on) -- stim2 from all trials,
          stim1 from TRAINING trials only, so held-out test trials never enter the
          basis (leakage-clean);
      (3) projects trials into the top-k of each cue's own basis;
      (4) fits Procrustes on the two 50 x k mean configs and returns the
          transformed test trials, the partial transforms, the shuffle control,
          and the residual shape distances. pro_decoding then decodes them with a
          stim2-trained ECOC SVM.

    Separate per-cue PCA (not joint): under the idealized claim X2 = X1*R, each
    cue's own PCA gives the same spectrum with eigenvectors related by R, so the
    top-k score configs coincide and PT reduces to a near-identity orthogonal map
    -- a clean, sensitive dimension sweep. The recovered T is therefore a RESIDUAL
    alignment after the change of basis, and the pdist columns are residual shape
    distances matched to these curves (not the full cross-cue rotation magnitude).

Output (per pair file): a 1 x numel(k_list) cell array; each cell is a struct:
    .pc_num   scalar, number of retained PCs k
    .stim1    stim1 condition label ('ac'/'ec'/'ex'), transformed
    .stim2    stim2 condition label, trained
    .acc      [trial_sample_repeat x 8] decoding accuracies, columns:
                1 genAcc                   stim2 within-cue 10-fold acc (per split)
                2 accuracy_s1              transfer, no transform
                3 accuracy_transformed_s1  transfer, full b*T + c
                4 accuracy_shuffle         rotation-only broken-correspondence null
                5 accuracy_onlyscale       transfer, b only
                6 accuracy_onlyrotation    transfer, T only  (matches col 4's null)
                7 accuracy_onlytranslation transfer, +c only
                8 accuracy_non_transfer_control  chance reference on the no-transform preds
    .pdist    [trial_sample_repeat x 3] residual shape distances on the 50 x k
              mean configs, all sharing MATLAB's denominator (scale of the stim2
              target), columns:
                1 d_before   no rotation      (T = I, b = 1)
                2 d_rot       rotation only    (b = 1, matches the test*T branch)
                3 d_full      MATLAB procrustes d (optimal b, T, c) -- reference floor

Note on cost: as in the pr script, pro_decoding retrains the stim2 ECOC per call,
so stim2's split-independent projection/model are recomputed each split. This is
redundant but kept for structural consistency; genAcc therefore carries the usual
per-split k-fold jitter rather than being constant.

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat
        Contains 'three_stim_array': cell array of [trial x neuron] matrices.
    - timewindow: two-element vector [start_ms end_ms] for spike counts.

Outputs:
    - MAT-files saved to:
      results/decoding_outputs/procrustes_decoding_basic_incremental_pc_results/<monkey>/<area>/
      Files: acec_results.mat, ecex_results.mat, acex_results.mat,
             ecac_results.mat, exec_results.mat, exac_results.mat
%}

%% Initialize MATLAB environment
clc; clear;
rng(1);

%% Configuration
settings   = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};  % all four in one run
timewindow = [330 630];       % e.g. early: [340 410], late: [410 480]
trial_sample_repeat = 100;    % trial-split repeats (absorbs the dropped neuron-repeat budget)

pairs      = {'ac','ec'; 'ec','ex'; 'ac','ex'; 'ec','ac'; 'ex','ec'; 'ex','ac'};
pair_names = {'acec','ecex','acex','ecac','exec','exac'};

%% Iterate over monkey / area settings
for si = 1:size(settings,1)
    monkey = settings{si,1};
    vp     = settings{si,2};

    data_file = fullfile('..','..','neuronal_data', monkey, vp, sprintf('%s_%s_allstim.mat', monkey, vp));
    tmp = load(data_file, 'three_stim_array');
    spike_data = tmp.three_stim_array;

    save_path = fullfile('..','..','results','decoding_outputs', ...
                         'procrustes_decoding_basic_incremental_pc_results', monkey, vp);
    if ~exist(save_path, 'dir'), mkdir(save_path); end

    %% Load trial-by-trial data
    nconds = numel(spike_data);
    data_trial  = cell(1,nconds);
    label_trial = cell(1,nconds);
    for i = 1:nconds
        [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, timewindow);
    end
    labels = label_trial{1};

    %% Build PC-dimension list (signal manifold rank <= 49 regardless of N)
    N = size(data_trial{1}, 2);
    k_list = 1:min(49, N);

    fprintf('=== %s %s (N = %d, k = 1..%d) ===\n', monkey, vp, N, k_list(end));
    tic;
    % Run decoding for each transformation pair
    for pi = 1:size(pairs,1)
        res   = incre_decoding(pairs{pi,1}, pairs{pi,2}, data_trial, labels, k_list, trial_sample_repeat);
        vname = [pair_names{pi} '_results'];
        out.(vname) = res;                                       
        save(fullfile(save_path,[vname '.mat']), '-struct', 'out', vname);
        clear out;
    end
    toc;
end

%%
function results = incre_decoding(stim1,stim2,data_trial,labels,k_list,trial_sample_repeat)
% stim1 is transformed, stim2 is the training data

pair_wise_data_trial = pair_pcaloader(data_trial,stim1,stim2);

%load data
stim1_data = pair_wise_data_trial(1:500,:);
stim2_data = pair_wise_data_trial(501:1000,:);

% zscore (per unit, over the full set of trials)
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);

results = cell(1,length(k_list)); % each cell stores results for one PC-dimension decoding

for k_sequence = 1:length(k_list)
    k = k_list(k_sequence);
    disp([stim1,stim2,' k=',num2str(k)])

    acc_block   = zeros(trial_sample_repeat,8);
    pdist_block = zeros(trial_sample_repeat,3);

    parfor trial_repeat = 1:trial_sample_repeat
        [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,shuffle_transformed,partially_transformed_stim1,pdist_row] = data_trial_sampler(stim1_data,stim2_data,labels,k);
        [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_shuffle,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control] = pro_decoding(training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,shuffle_transformed,partially_transformed_stim1);

        acc_block(trial_repeat,:)   = [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_shuffle,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control];
        pdist_block(trial_repeat,:) = pdist_row;
    end

    one = struct();
    one.pc_num = k;
    one.stim1  = stim1;
    one.stim2  = stim2;
    one.acc    = acc_block;    % trial_sample_repeat x 8
    one.pdist  = pdist_block;  % trial_sample_repeat x 3
    results{k_sequence} = one;
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

function [mu, V] = mean_basis(means_50xN)
% Separate per-cue PCA basis from the centered condition-mean manifold.
% mu : 1 x N grand mean across the 50 stimuli (centering used for projection)
% V  : N x r variance-ordered right singular vectors (principal directions).
% SVD of the centered means is scale-invariant and rank-safe when N > 50.
mu = mean(means_50xN, 1);
Xc = means_50xN - mu;
[~, ~, V] = svd(Xc, 'econ');
end

function [training_data,training_label,stim1_test_data,stim1_test_label,transformed_stim1,shuffle_transformed,partially_transformed_stim1,pdist_row] = data_trial_sampler(stim1_data_sample,stim2_data_sample,labels,k)
% sample trial
tial_out_of_ten = 1;
valid_trial_idx = sample_trial(tial_out_of_ten); % 10 fold
training_trial_idx = setdiff([1:500],valid_trial_idx);

% take the mean (stim2 from all trials; stim1 from training trials only)
stim2_target_data_averged     = take_average(stim2_data_sample, 10);
stim1_2transform_data_averged = take_average(stim1_data_sample(training_trial_idx,:), 10 - tial_out_of_ten);

% separate per-cue PCA basis on the centered means, top-k
[mu2, V2] = mean_basis(stim2_target_data_averged);  V2 = V2(:,1:k);
[mu1, V1] = mean_basis(stim1_2transform_data_averged);  V1 = V1(:,1:k);

% SVM training data, all stim2 data (projected into stim2's top-k)
training_data  = (stim2_data_sample - mu2) * V2;
training_label = labels;

% test data (projected into stim1's top-k)
stim1_test_data  = (stim1_data_sample(valid_trial_idx,:) - mu1) * V1;
stim1_test_label = labels(valid_trial_idx,:);

% mean configs in PC space (target = stim2 class centroids; source = stim1 train)
s2_means_k = (stim2_target_data_averged     - mu2) * V2;
s1_means_k = (stim1_2transform_data_averged - mu1) * V1;

[d_full,Z,transform] = procrustes(s2_means_k, s1_means_k);
t_matrix = transform.T;
b = transform.b;
c = transform.c;
transformed_stim1 = b*stim1_test_data*t_matrix+c;

% partially transformed data - only 1 transofrm
scaling_only_transformed_stim1     = b*stim1_test_data;
rotation_only_transformed_stim1    = stim1_test_data*t_matrix;
translation_only_transformed_stim1 = stim1_test_data + c;
partially_transformed_stim1 = {scaling_only_transformed_stim1,rotation_only_transformed_stim1,translation_only_transformed_stim1};

% shuffle control (correspondence destroyed for both mean matrix and test trials)
shuffle_idx = randperm(size(s1_means_k,1));
% shuffle the mean matrix and learn the (broken-correspondence) rotation
toy_averaged = s1_means_k(shuffle_idx,:);
[d_fake,Z_fake,t_fake] = procrustes(s2_means_k, toy_averaged);
t_fakematrix = t_fake.T;
% shuffle the held-out test trials the same way, then apply rotation only
stim1_test_data_shuffled = stim1_test_data(shuffle_idx,:);
shuffle_transformed = stim1_test_data_shuffled*t_fakematrix;   % rotation only (no toy_b, no toy_c)

% residual Procrustes shape distances (b = 1), shared stim2 denominator
A0  = s2_means_k - mean(s2_means_k,1);
B0  = s1_means_k - mean(s1_means_k,1);
ssA = sum(A0(:).^2);
d_before = sum((A0 - B0         ).^2,'all')/ssA;   % no rotation
d_rot    = sum((A0 - B0*t_matrix).^2,'all')/ssA;   % rotation only
pdist_row = [d_before, d_rot, d_full];

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

function pair_wise_pca_data = pair_pcaloader(dca_data,stim1,stim2)
% use the input stimulus name to load pair_wise dca data

name2idx = struct('ac',1,'ec',2,'ex',3);

stim1_data = dca_data{name2idx.(stim1)};
stim2_data = dca_data{name2idx.(stim2)};
pair_wise_pca_data = [stim1_data;stim2_data];
end

function [genAcc,accuracy_s1,accuracy_transformed_s1,accuracy_shuffle,accuracy_onlyscale,accuracy_onlyrotation,accuracy_onlytranslation,accuracy_non_transfer_control] ...
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
% rotation-only broken-correspondence shuffle control
predicted_Label_shuffle = predict(stim2_model,shuffle_transformed);
accuracy_shuffle = sum(stim1_test_label == predicted_Label_shuffle)/length(predicted_Label_shuffle);
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
