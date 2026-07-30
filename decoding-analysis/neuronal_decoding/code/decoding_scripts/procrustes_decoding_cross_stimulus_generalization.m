%{
Filename: procrustes_decoding_cross_stimulus_generalization.m
Author: Zitong Wang (CCGP / cross-stimulus generalization)
Date: 2026-07-30

Description:
    CCGP-style Procrustes transfer decoding with a STIMULUS hold-out design.
    Instead of fitting the cross-cue alignment on all 50 stimulus means and
    testing on held-out TRIALS of those same stimuli, the alignment (ROTATION
    ONLY) is estimated from a 40-stimulus training subset and applied to the 10
    held-out stimuli the rotation never saw. This tests whether the cue-1 ->
    cue-2 alignment is a GLOBAL property of the manifold (a rotation learned from
    part of it places the rest correctly) rather than something that requires
    seeing every stimulus.

    Everything is scored on the 10 held-out stimuli with a 10-way ECOC classifier
    trained on the TARGET cue (stim2). Only the rotation component T of the
    Procrustes solution is applied to the source test trials (X*T); the uniform
    scale b and translation c are NOT applied, matching the rotation-only claim.

    (Descends from procrustes_decoding_basic_incremental_pr.m; the participation-
    ratio machinery has been removed -- PR is not part of this analysis and lives
    in the _pr script, joinable by neuron_num if the gap-vs-PR view is wanted.)

Resampling structure (three independent knobs):
    n_partition           outer CCGP resampling: random 40/10 stimulus splits.
                          THIS is the dominant variance source (which 10 are held
                          out matters most), so it is its own explicit loop and
                          the one to keep large. parfor runs over this axis.
    neuron_sample_repeat  neuron subsets drawn per partition (subset overlap makes
                          this axis low-variance; small counts suffice).
    trial_sample_repeat   1-of-10 trial hold-outs per (partition, neuron subset);
                          smooths the single held-out source trial. Second-order.

    Efficiency: the target-H classifier, self_decode, and the fit-40 rotations
    (T_gen, T_rand) depend only on (partition, neuron subset) -- NOT on the trial
    hold-out -- so they are computed ONCE per (partition, neuron subset) and
    reused across the trial loop. Only T_ceil / T_floor and the source test/anchor
    trials are redrawn per trial. Consequence: self_decode is CONSTANT across the
    trial_sample_repeat rows of a (partition, neuron) cell -- collapse on
    partition_id x neuron_rep_id before treating it as independent samples.

Regimes and controls (all rotation-only, all scored on the 10 held-out stimuli):
    self_decode   within-target-cue 10-way k-fold CV. Absolute (no-transfer)
                  ceiling for these 10 stimuli. Chance = 1/10.
    no_transform  target classifier applied to UNtransformed source test trials.
    pt_gen        MAIN. Rotation fit on the 40 training-stimulus means (true
                  correspondence), applied to the 10 held-out source test trials.
                  The generalization / CCGP measure. Read RELATIVE to the bracket
                  below: pt_gen ~ pt_ceiling => the rotation generalizes.
    pt_ceiling    Rotation fit on the 10 held-out stimulus means (true
                  correspondence, from the 9 train trials), applied to the same 10.
                  In-sample PT ceiling.
    pt_floor      Rotation fit on the 10 held-out stimuli under a SELF-CONSISTENT
                  correspondence shuffle: one permutation pi shuffles BOTH the
                  source mean matrix (used to fit the rotation) AND the source test
                  trials, while labels stay in true order. The correct null for
                  "can the transform force a mapping when true correspondence is
                  destroyed?".
                  NOTE: OVERFITTING-MATCHED null, not a chance floor. With N large
                  relative to 10 anchor points a rotation has enough d.o.f. to
                  align an arbitrary bijection, so pt_floor rides toward pt_ceiling
                  at large N and only separates in the constrained (small-N)
                  regime. The ceiling<->floor gap vs neuron count is the
                  overfitting diagnostic.
    rand_rot_40   Rotation fit on the 40 with a NON-self-consistent shuffle (means
                  shuffled, test unshuffled) applied to the 10. Genuine chance
                  floor, kept as a trivial sanity reference. (A self-consistent
                  null is undefined at the 40-fit scale because fit set != test
                  set -- which is why the floor/ceiling bracket is built on the 10.)
    chance        Random-label baseline (~1/10).

    z-scoring: per unit, once over all 500 trials (50-stimulus calibration frame).
    Treated as fixed sensor gain, not a learned alignment parameter, so estimating
    it on all stimuli is a calibration choice, not correspondence/label leakage.
    To switch to a strict train-only (40-stimulus) frame, move the zscore() calls
    into score_one_trial and fit on F; this would make the frame partition-
    dependent (deliberately kept global here).

Output (per pair file): a 1 x numel(neuron_list) cell array; each cell a struct:
    .neuron_num     scalar, number of sampled units
    .stim1, .stim2  condition labels; stim1 = source (transformed), stim2 = target
    .acc            [(n_partition*neuron_sample_repeat*trial_sample_repeat) x 7]
                    accuracies, columns in this FIXED order:
                      1 self_decode   2 no_transform   3 pt_gen   4 pt_ceiling
                      5 pt_floor      6 rand_rot_40     7 chance
    .partition_id   [(...) x 1] which 40/10 stimulus draw each acc row came from
    .neuron_rep_id  [(...) x 1] which neuron subset (within its partition) per row

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat
        Contains 'three_stim_array': cell array of [trial x neuron] matrices.
    - timewindow: two-element vector [start_ms end_ms] for spike counts.

Execution:
    One click runs all four monkey x area combinations (FR/KO x V1/V2) in a
    single batch loop; each combination loads its own data and neuron-count list
    and writes its own results folder.

Outputs:
    - MAT-files saved to:
      results/decoding_outputs/procrustes_decoding_cross_stimulus_generalization_results/<monkey>/<area>/
      Files: acec_results.mat, ecex_results.mat, acex_results.mat,
             ecac_results.mat, exec_results.mat, exac_results.mat
%}

%% Initialize MATLAB environment
clc; clear;

%% Batch configuration -- one click over every monkey x area combination
combos = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};  % {monkey, area} per row
timewindow = [330 630];  % e.g., early: [340 410], late: [410 480]

% max units available per monkey/area (indexes into neuron_code via map)
neuron_code = {[48]; [48 112]; [109]; [109 146]};
map = struct('FRV1',1, 'FRV2',2, 'KOV1',3, 'KOV2',4);

%% Resampling parameters (three independent knobs -- see header)
n_partition          = 100;  % random 40/10 stimulus splits (dominant variance axis)
neuron_sample_repeat = 10;   % neuron subsets per partition
trial_sample_repeat  = 10;   % 1-of-10 trial hold-outs per (partition, neuron subset)
rng(1);

total_timer = tic;
for c = 1:size(combos,1)
    monkey = combos{c,1};   % 'FR' or 'KO'
    vp     = combos{c,2};   % 'V1' or 'V2'
    fprintf('=== %s %s ===\n', monkey, vp);
    combo_timer = tic;

    % ---- load trial-by-trial data for this monkey/area ----
    data_file = fullfile('..','..','neuronal_data', monkey, vp, sprintf('%s_%s_allstim.mat', monkey, vp));
    tmp = load(data_file, 'three_stim_array');
    spike_data = tmp.three_stim_array;

    save_path = fullfile('..','..','results','decoding_outputs','procrustes_decoding_cross_stimulus_generalization_results', monkey, vp);
    if ~exist(save_path, 'dir'), mkdir(save_path); end

    nconds = numel(spike_data);
    data_trial  = cell(1,nconds);
    label_trial = cell(1,nconds);
    for i = 1:nconds
        [data_trial{i}, label_trial{i}] = multiclass_svmloader_PT(spike_data{i}, timewindow);
    end

    % ---- neuron-count list (uniform 10-unit steps) for this monkey/area ----
    max_list    = neuron_code{map.(strcat(monkey, vp))};
    neuron_list = unique([10:10:max(max_list), max_list]);
    labels      = label_trial{1};

    % ---- run all six ordered pairs (stim1 = source, stim2 = target) ----
    acec_results = cross_stim_decoding('ac','ec', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat, n_partition);
    save(fullfile(save_path,'acec_results.mat'),'acec_results');

    ecex_results = cross_stim_decoding('ec','ex', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat, n_partition);
    save(fullfile(save_path,'ecex_results.mat'),'ecex_results');

    acex_results = cross_stim_decoding('ac','ex', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat, n_partition);
    save(fullfile(save_path,'acex_results.mat'),'acex_results');

    ecac_results = cross_stim_decoding('ec','ac', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat, n_partition);
    save(fullfile(save_path,'ecac_results.mat'),'ecac_results');

    exec_results = cross_stim_decoding('ex','ec', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat, n_partition);
    save(fullfile(save_path,'exec_results.mat'),'exec_results');

    exac_results = cross_stim_decoding('ex','ac', data_trial, labels, neuron_list, neuron_sample_repeat, trial_sample_repeat, n_partition);
    save(fullfile(save_path,'exac_results.mat'),'exac_results');

    fprintf('    %s %s done in %.1f min\n', monkey, vp, toc(combo_timer)/60);
end
fprintf('all combinations done in %.1f min\n', toc(total_timer)/60);

%%
function results = cross_stim_decoding(stim1,stim2,data_trial,labels,neuron_num_list,neuron_sample_repeat,trial_sample_repeat,n_partition)
% stim1 = source (transformed), stim2 = target (classifier trained here).
% labels is unused: held-out stimulus ids are derived from the 40/10 partition
% inside draw_partition. Kept in the signature for call-site parity.

pair_wise_data_trial = pair_pcaloader(data_trial,stim1,stim2);

% load data
stim1_data = pair_wise_data_trial(1:500,:);     % source
stim2_data = pair_wise_data_trial(501:1000,:);  % target

% zscore (per unit, over the full set of trials -- calibration frame, see header)
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);

results = cell(1,length(neuron_num_list)); % each cell stores results for one neuron-number decoding

for neuron_squence = 1:length(neuron_num_list)
    disp([stim1,stim2,num2str(neuron_squence)])

    % one sliced block per partition (parfor over the CCGP resampling axis)
    acc_blocks  = cell(n_partition,1);
    part_blocks = cell(n_partition,1);
    neu_blocks  = cell(n_partition,1);
    rows_per_part = neuron_sample_repeat*trial_sample_repeat;

    parfor part = 1:n_partition
        % draw one 40/10 stimulus partition, fixed across the neuron & trial loops
        [H,F] = draw_partition();

        acc_p = zeros(rows_per_part,7);
        neu_p = zeros(rows_per_part,1);

        for neuron_repeat = 1:neuron_sample_repeat
            sampled_neuron_idx = sample_neuron(neuron_num_list,neuron_squence);
            src = stim1_data(:,sampled_neuron_idx);
            tgt = stim2_data(:,sampled_neuron_idx);

            % ---- computed ONCE per (partition, neuron subset): trial-independent ----
            [model, self_decode] = train_target_classifier(tgt, H);
            source_F_mean = cond_mean_all(src, F);
            target_F_mean = cond_mean_all(tgt, F);
            target_H_mean = cond_mean_all(tgt, H);
            T_gen  = rot_only(target_F_mean, source_F_mean);                       % fit-40 true
            T_rand = rot_only(target_F_mean, source_F_mean(randperm(numel(F)),:)); % fit-40 chance

            % ---- trial hold-out loop (redraws only T_ceil/T_floor + source trials) ----
            for trial_repeat = 1:trial_sample_repeat
                row = (neuron_repeat-1)*trial_sample_repeat + trial_repeat;
                acc_p(row,:) = score_one_trial(model, self_decode, src, H, target_H_mean, T_gen, T_rand);
                neu_p(row)   = neuron_repeat;
            end
        end

        acc_blocks{part}  = acc_p;
        part_blocks{part} = repmat(part, rows_per_part, 1);
        neu_blocks{part}  = neu_p;
    end

    one = struct();
    one.neuron_num    = neuron_num_list(neuron_squence);
    one.stim1         = stim1;
    one.stim2         = stim2;
    one.acc           = cat(1, acc_blocks{:});    % (n_partition*nrep*trep) x 7
    one.partition_id  = cat(1, part_blocks{:});   % which 40/10 draw
    one.neuron_rep_id = cat(1, neu_blocks{:});    % which neuron subset within the partition
    results{neuron_squence} = one;
end

end

%% ---------- one 40/10 stimulus partition ----------
function [H,F] = draw_partition()
perm = randperm(50);
H = sort(perm(1:10));      % 10 held-out (test) stimuli
F = sort(perm(11:50));     % 40 fit stimuli
end

%% ---------- target classifier (trained once per partition x neuron subset) ----------
function [model, self_decode] = train_target_classifier(tgt, H)
% 10-way ECOC on all target-cue trials of the held-out stimuli; self_decode is
% the within-cue 10-fold CV accuracy (no-transfer ceiling).
[Xtr, ytr] = stim_trials(tgt, H);      % 100 x n, labels = stim id
model       = fitcecoc(Xtr, ytr);
cv          = crossval(model,'KFold',10);
self_decode = 1 - kfoldLoss(cv);
end

%% ---------- one trial hold-out -> 7 rotation-only measures ----------
function acc = score_one_trial(model, self_decode, src, H, target_H_mean, T_gen, T_rand)
% T_gen / T_rand are pre-fit (trial-independent); T_ceil / T_floor are fit here
% (they depend on the 9-trial source mean). All applied transforms are X*T.
testIdx        = randi(10, 10, 1);                 % held-out trial per H stim
source_H_test  = pick_trial(src, H, testIdx);      % 10 x n (transfer test set)
source_H_train = cond_mean_excl(src, H, testIdx);  % 10 x n (rotation anchor)
test_label     = H(:);

T_ceil = rot_only(target_H_mean, source_H_train);

% self-consistent shuffle on the 10 (same pi on fit means AND test trials)
pi10               = randperm(10);
T_floor            = rot_only(target_H_mean, source_H_train(pi10,:));
source_H_test_shuf = source_H_test(pi10,:);        % labels stay = test_label

no_transform = mean(predict(model, source_H_test)                == test_label);
pt_gen       = mean(predict(model, source_H_test      * T_gen)   == test_label);
pt_ceiling   = mean(predict(model, source_H_test      * T_ceil)  == test_label);
pt_floor     = mean(predict(model, source_H_test_shuf * T_floor) == test_label);
rand_rot_40  = mean(predict(model, source_H_test      * T_rand)  == test_label);
chance       = mean(test_label(randperm(10))                     == test_label);

% column order: 1 self_decode 2 no_transform 3 pt_gen 4 pt_ceiling
%               5 pt_floor    6 rand_rot_40  7 chance
acc = [self_decode, no_transform, pt_gen, pt_ceiling, pt_floor, rand_rot_40, chance];
end

%% ---------- rotation-only Procrustes ----------
function T = rot_only(target_means, source_means)
% Procrustes maps source_means -> target_means; return the orthogonal
% (rotation/reflection) component T only. Scale b and translation c are dropped,
% matching the rotation-only convention (X*T) used elsewhere in the pipeline.
[~,~,tf] = procrustes(target_means, source_means);
T = tf.T;
end

%% ---------- per-stimulus indexing helpers (500 rows = 50 stim x 10 trials) ----------
function [X, lab] = stim_trials(data, stim_ids)
% all 10 trials of each stimulus in stim_ids, labelled by stimulus id
k = numel(stim_ids);
X   = zeros(k*10, size(data,2));
lab = zeros(k*10,1);
for i = 1:k
    s    = stim_ids(i);
    rows = (s-1)*10 + (1:10);
    X((i-1)*10 + (1:10), :) = data(rows,:);
    lab((i-1)*10 + (1:10))  = s;
end
end

function X = pick_trial(data, stim_ids, trialIdx)
% one specified trial (trialIdx(i) in 1..10) per stimulus
k = numel(stim_ids);
X = zeros(k, size(data,2));
for i = 1:k
    s = stim_ids(i);
    X(i,:) = data((s-1)*10 + trialIdx(i), :);
end
end

function M = cond_mean_excl(data, stim_ids, exclIdx)
% per-stimulus mean over the 9 trials EXCLUDING exclIdx(i)
k = numel(stim_ids);
M = zeros(k, size(data,2));
for i = 1:k
    s          = stim_ids(i);
    rows       = (s-1)*10 + (1:10);
    rows(exclIdx(i)) = [];
    M(i,:)     = mean(data(rows,:), 1);
end
end

function M = cond_mean_all(data, stim_ids)
% per-stimulus mean over all 10 trials
k = numel(stim_ids);
M = zeros(k, size(data,2));
for i = 1:k
    s      = stim_ids(i);
    rows   = (s-1)*10 + (1:10);
    M(i,:) = mean(data(rows,:), 1);
end
end

%% ---------- carried over from procrustes_decoding_basic_incremental_pr.m ----------
function sampled_neuron_idx = sample_neuron(neuron_num_list,sequence)
% neuron sampling. sequence: the idx of how many neurons to be sampled.
origin_neuron_num = max(neuron_num_list);
sampled_neuron_num = neuron_num_list(sequence);
sampled_neuron_idx = randperm(origin_neuron_num,sampled_neuron_num);
% disp([origin_neuron_num,sampled_neuron_num])  % silenced: fires n_partition*nrep times under parfor
end

function pair_wise_pca_data = pair_pcaloader(dca_data,stim1,stim2)
% use the input stimulus name to load pair_wise dca data
name2idx = struct('ac',1,'ec',2,'ex',3);
stim1_data = dca_data{name2idx.(stim1)};
stim2_data = dca_data{name2idx.(stim2)};
pair_wise_pca_data = [stim1_data;stim2_data];
end
