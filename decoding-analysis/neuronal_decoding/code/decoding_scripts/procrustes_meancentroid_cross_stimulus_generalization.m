%{
Filename: procrustes_meancentroid_cross_stimulus_generalization.m
Author: Zitong Wang (CCGP / cross-stimulus generalization -- mean level)
Date: 2026-08-01

Description:
    MEAN-LEVEL companion to procrustes_decoding_cross_stimulus_generalization.m.
    Same generalization question -- does a rotation-only cross-cue alignment fit
    on part of the manifold place the rest correctly -- but softened: instead of
    learning the rotation from means and then decoding NOISY SINGLE TRIALS through
    an SVM, we decode the held-out MEANS themselves by nearest centroid. This
    strips trial noise and asks the structural question directly: after applying
    the rotation, does each held-out source mean land on its own target mean?

    Rotation is fit on the (50 - n_stim_hold) fit stimuli (F) and applied,
    ROTATION ONLY (M1_H * T; no scale b, no translation c), to the n_stim_hold
    held-out source means (H). Because the per-unit z-score centers each cue's
    full 50-stimulus cloud at the origin, the dropped Procrustes translation is
    ~0, so nearest-centroid (which is translation-sensitive) is not thrown off,
    and the convention stays consistent with the trial-level script.

Two reference sets (both scored on the same transformed held-out source means):
    GLOBAL  nearest among ALL 50 target means. The fit-stimulus targets act as
            distractors, so this tests GLOBAL placement. Chance = 1/50.
    WITHIN  nearest among only the n_stim_hold held-out target means. The easier
            rank-within-the-held-out test. Chance = 1/n_stim_hold.

Three variants under each reference set:
    gen           rotation fit on F (true correspondence) -- the MAIN measure.
    no_transform  raw source means, no rotation. Baseline: the SAME neurons are
                  recorded across cues, so this is how much cue-invariance already
                  sits in the raw responses. gen > no_transform => the rotation
                  adds alignment; gen < no_transform => a mis-generalizing rotation
                  actively hurts.
    rand          rotation fit on F with a shuffled correspondence -- chance floor.

Deliberately omitted (vs the trial-level bracket): in-sample pt_ceiling and the
    self-consistent pt_floor. At the mean level with rotation-only, a rotation fit
    ON the held-out means aligns them near-perfectly whenever N >= n_stim_hold
    (always true here), so both saturate to ~1 on a neuron axis and carry no
    information. They only become graded under a PC-dimension sweep (future work);
    the `gen` rotation, fit on the large F set, is genuinely constrained and its
    neuron-sweep curve is real signal.

Resampling (two knobs -- no trial loop, means use all 10 trials):
    n_partition           random splits holding out n_stim_hold stimuli (dominant
                          variance axis; parfor runs over it).
    neuron_sample_repeat  neuron subsets drawn per partition.

Output (per pair file): a 1 x numel(neuron_list) cell array; each cell a struct:
    .neuron_num     scalar, number of sampled units
    .stim1, .stim2  condition labels; stim1 = source (transformed), stim2 = target
    .acc            [(n_partition*neuron_sample_repeat) x 6] accuracies, columns:
                      1 gen_global  2 notransform_global  3 rand_global
                      4 gen_within  5 notransform_within  6 rand_within
    .partition_id   [(...) x 1] which stimulus draw (partition) each acc row is
    .neuron_rep_id  [(...) x 1] which neuron subset (within its partition) per row
    .stim_hold_num  n_stim_hold (chance = 1/50 for global cols, 1/stim_hold_num
                    for within cols)

Inputs:
    - neuronal_data/<monkey>/<area>/<monkey>_<area>_allstim.mat  ('three_stim_array')
    - timewindow: two-element vector [start_ms end_ms] for spike counts.

Execution:
    One click over all four monkey x area combinations (FR/KO x V1/V2).

Outputs:
    - MAT-files saved to:
      results/decoding_outputs/procrustes_meancentroid_cross_stimulus_generalization_results/HoldStim<n>/<monkey>/<area>/
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

%% Resampling parameters (two knobs -- no trial loop; see header)
n_partition          = 100;  % random stimulus splits (dominant variance axis)
neuron_sample_repeat = 10;   % neuron subsets per partition
n_stim_hold          = 2;    % number of held-out stimuli (fit set = 50 - n_stim_hold)
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

    save_path = fullfile('..','..','results','decoding_outputs', ...
        'procrustes_meancentroid_cross_stimulus_generalization_results', ...
        ['HoldStim',num2str(n_stim_hold)], monkey, vp);
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
    acec_results = meancentroid_decoding('ac','ec', data_trial, labels, neuron_list, neuron_sample_repeat, n_partition, n_stim_hold);
    save(fullfile(save_path,'acec_results.mat'),'acec_results');

    ecex_results = meancentroid_decoding('ec','ex', data_trial, labels, neuron_list, neuron_sample_repeat, n_partition, n_stim_hold);
    save(fullfile(save_path,'ecex_results.mat'),'ecex_results');

    acex_results = meancentroid_decoding('ac','ex', data_trial, labels, neuron_list, neuron_sample_repeat, n_partition, n_stim_hold);
    save(fullfile(save_path,'acex_results.mat'),'acex_results');

    ecac_results = meancentroid_decoding('ec','ac', data_trial, labels, neuron_list, neuron_sample_repeat, n_partition, n_stim_hold);
    save(fullfile(save_path,'ecac_results.mat'),'ecac_results');

    exec_results = meancentroid_decoding('ex','ec', data_trial, labels, neuron_list, neuron_sample_repeat, n_partition, n_stim_hold);
    save(fullfile(save_path,'exec_results.mat'),'exec_results');

    exac_results = meancentroid_decoding('ex','ac', data_trial, labels, neuron_list, neuron_sample_repeat, n_partition, n_stim_hold);
    save(fullfile(save_path,'exac_results.mat'),'exac_results');

    fprintf('    %s %s done in %.1f min\n', monkey, vp, toc(combo_timer)/60);
end
fprintf('all combinations done in %.1f min\n', toc(total_timer)/60);

%%
function results = meancentroid_decoding(stim1,stim2,data_trial,labels,neuron_num_list,neuron_sample_repeat,n_partition,n_stim_hold) %#ok<INUSD>
% stim1 = source (transformed), stim2 = target (centroids live here).
% labels is unused: held-out stimulus ids come from the partition. Kept for parity.

pair_wise_data_trial = pair_pcaloader(data_trial,stim1,stim2);

% load data
stim1_data = pair_wise_data_trial(1:500,:);     % source
stim2_data = pair_wise_data_trial(501:1000,:);  % target

% zscore (per unit, over all 500 trials -- centers each cue's cloud at origin,
% which is what makes the dropped Procrustes translation ~0; see header)
stim1_data = zscore(stim1_data);
stim2_data = zscore(stim2_data);

results = cell(1,length(neuron_num_list));

for neuron_squence = 1:length(neuron_num_list)
    disp([stim1,stim2,num2str(neuron_squence)])

    acc_blocks  = cell(n_partition,1);
    part_blocks = cell(n_partition,1);
    neu_blocks  = cell(n_partition,1);

    parfor part = 1:n_partition
        % draw one stimulus partition (n_stim_hold held out), fixed across neuron loop
        [H,F] = draw_partition(n_stim_hold);

        acc_p = zeros(neuron_sample_repeat,6);
        for neuron_repeat = 1:neuron_sample_repeat
            sampled_neuron_idx = sample_neuron(neuron_num_list,neuron_squence);
            src = stim1_data(:,sampled_neuron_idx);
            tgt = stim2_data(:,sampled_neuron_idx);

            % ---- condition means (all 10 trials; no trial hold-out at mean level) ----
            M1_F   = cond_mean_all(src, F);        % (50-h) x n  source fit means
            M2_F   = cond_mean_all(tgt, F);        % (50-h) x n  target fit means
            M1_H   = cond_mean_all(src, H);        % h x n       source held-out means (queries)
            M2_H   = cond_mean_all(tgt, H);        % h x n       target held-out means (within refs)
            M2_all = cond_mean_all(tgt, 1:50);     % 50 x n      all target means (global refs)

            % ---- rotation-only Procrustes fit on F ----
            T_gen  = rot_only(M2_F, M1_F);                        % true correspondence
            T_rand = rot_only(M2_F, M1_F(randperm(numel(F)),:));  % shuffled correspondence (floor)

            Q_gen  = M1_H * T_gen;    % transformed held-out source means
            Q_rand = M1_H * T_rand;
            Q_raw  = M1_H;            % no transform

            % ---- GLOBAL reference (all 50 target means; chance 1/50) ----
            gen_g = nearest_centroid(Q_gen,  H, M2_all, (1:50)');
            nt_g  = nearest_centroid(Q_raw,  H, M2_all, (1:50)');
            rnd_g = nearest_centroid(Q_rand, H, M2_all, (1:50)');

            % ---- WITHIN reference (held-out target means; chance 1/n_stim_hold) ----
            gen_w = nearest_centroid(Q_gen,  H, M2_H, H(:));
            nt_w  = nearest_centroid(Q_raw,  H, M2_H, H(:));
            rnd_w = nearest_centroid(Q_rand, H, M2_H, H(:));

            % column order: 1 gen_global 2 notransform_global 3 rand_global
            %               4 gen_within 5 notransform_within 6 rand_within
            acc_p(neuron_repeat,:) = [gen_g, nt_g, rnd_g, gen_w, nt_w, rnd_w];
        end

        acc_blocks{part}  = acc_p;
        part_blocks{part} = repmat(part, neuron_sample_repeat, 1);
        neu_blocks{part}  = (1:neuron_sample_repeat)';
    end

    one = struct();
    one.neuron_num    = neuron_num_list(neuron_squence);
    one.stim1         = stim1;
    one.stim2         = stim2;
    one.acc           = cat(1, acc_blocks{:});    % (n_partition*nrep) x 6
    one.partition_id  = cat(1, part_blocks{:});   % which draw
    one.neuron_rep_id = cat(1, neu_blocks{:});    % which neuron subset within the partition
    one.stim_hold_num = n_stim_hold;
    results{neuron_squence} = one;
end

end

%% ---------- one stimulus partition (n_stim_hold held out) ----------
function [H,F] = draw_partition(n)
perm = randperm(50);
H = sort(perm(1:n));      % n held-out (test) stimuli
F = sort(perm(n+1:50));   % 50-n fit stimuli
end

%% ---------- nearest-centroid assignment ----------
function acc = nearest_centroid(Q, q_labels, C, c_labels)
% Assign each query row in Q to the nearest centroid row in C (Euclidean).
% Accuracy = fraction of queries whose nearest centroid carries the query's label.
D = pdist2(Q, C);          % nQ x nC
[~, j] = min(D, [], 2);
acc = mean(c_labels(j) == q_labels(:));
end

%% ---------- rotation-only Procrustes ----------
function T = rot_only(target_means, source_means)
% Procrustes maps source_means -> target_means; return the orthogonal
% (rotation/reflection) component T only. Scale b and translation c are dropped.
[~,~,tf] = procrustes(target_means, source_means);
T = tf.T;
end

%% ---------- per-stimulus mean helper (500 rows = 50 stim x 10 trials) ----------
function M = cond_mean_all(data, stim_ids)
% per-stimulus mean over all 10 trials, one row per stimulus (in stim_ids order)
k = numel(stim_ids);
M = zeros(k, size(data,2));
for i = 1:k
    s      = stim_ids(i);
    rows   = (s-1)*10 + (1:10);
    M(i,:) = mean(data(rows,:), 1);
end
end

%% ---------- carried over from procrustes_decoding_cross_stimulus_generalization.m ----------
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
