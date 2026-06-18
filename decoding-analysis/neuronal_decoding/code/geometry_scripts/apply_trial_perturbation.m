function data_trial_perturbed = apply_trial_perturbation(data_trial, labels, mode, opts)
%APPLY_TRIAL_PERTURBATION  Break noise correlations in pseudopopulation trial data.
%
%   Reorders trials within each (neuron, stimulus) block so that neurons no
%   longer share trial-aligned co-variability, while the marginal firing-rate
%   distribution per neuron per stimulus is preserved exactly. Intended to be
%   dropped in right after multiclass_svmloader_PT in decoding scripts so the
%   downstream pipeline (zscore, Procrustes, SVM, sampling) is unmodified.
%
%   data_trial_perturbed = APPLY_TRIAL_PERTURBATION(data_trial, labels, mode, opts)
%
% INPUTS
%   data_trial : [N_trials_total x N_neurons] firing rates grouped by
%                stimulus label (as returned by multiclass_svmloader_PT).
%   labels     : [N_trials_total x 1] stimulus labels (e.g. 1..50).
%   mode       : 'none'    - pass-through (no perturbation).
%                'shuffle' - independent random permutation of the 10 trials
%                            per neuron per stimulus. Matches the random-
%                            shuffle baseline in noise_correlation_analysis.m.
%                'affine'  - deterministic affine trial shift per neuron
%                            based on session membership. Matches the
%                            cyclic/affine shift in noise_correlation_analysis_affine.m.
%   opts       : (optional struct, required for 'affine')
%                .neuron_pattern_id : [N_neurons x 1] pattern index in 1..24.
%                .affine_patterns   : [24 x N_trials_per_stim] permutation table.
%                See build_affine_shift_config.m.
%
% OUTPUT
%   data_trial_perturbed : same size as data_trial. Row i still carries
%                          label labels(i); only the response at row i has
%                          been swapped with another trial of the same
%                          stimulus (per neuron, independently).
%
% NOTES
%   * RNG: this function does NOT reset the RNG. Call rng(...) in the caller
%     before invoking with mode='shuffle' if reproducibility is needed.
%   * mode='none' returns data_trial unchanged so the call site can sit
%     unconditionally inside the loader loop.

    if nargin < 3 || isempty(mode), mode = 'none'; end
    if nargin < 4,                    opts = struct(); end

    if strcmpi(mode, 'none')
        data_trial_perturbed = data_trial;
        return;
    end

    %% Reorganize into [N_stimuli x N_trials_per_stim x N_neurons]
    stim_ids  = unique(labels);
    N_stimuli = numel(stim_ids);
    N_neurons = size(data_trial, 2);

    % Require a uniform number of trials per stimulus (matches the
    % assumptions of the NC analysis and the decoding sampler).
    counts = arrayfun(@(s) sum(labels == s), stim_ids);
    assert(all(counts == counts(1)), ...
        'apply_trial_perturbation: unequal trial counts per stimulus: [%s]', ...
        num2str(counts(:)'));
    N_trials_per_stim = counts(1);

    data_by_stim = zeros(N_stimuli, N_trials_per_stim, N_neurons);
    for s = 1:N_stimuli
        data_by_stim(s, :, :) = data_trial(labels == stim_ids(s), :);
    end

    %% Apply the requested permutation
    switch lower(mode)
        case 'shuffle'
            % Independent random permutation of trials, per neuron per stim.
            for n = 1:N_neurons
                for s = 1:N_stimuli
                    data_by_stim(s, :, n) = ...
                        data_by_stim(s, randperm(N_trials_per_stim), n);
                end
            end

        case 'affine'
            assert(isfield(opts, 'neuron_pattern_id') && ...
                   isfield(opts, 'affine_patterns'), ...
                'apply_trial_perturbation: affine mode requires opts.neuron_pattern_id and opts.affine_patterns.');
            neuron_pattern_id = opts.neuron_pattern_id;
            affine_patterns   = opts.affine_patterns;
            assert(size(affine_patterns, 2) == N_trials_per_stim, ...
                'affine_patterns has %d trial columns; expected %d.', ...
                size(affine_patterns, 2), N_trials_per_stim);
            assert(numel(neuron_pattern_id) == N_neurons, ...
                'neuron_pattern_id length (%d) != N_neurons (%d).', ...
                numel(neuron_pattern_id), N_neurons);
            assert(all(neuron_pattern_id >= 1 & neuron_pattern_id <= size(affine_patterns, 1)), ...
                'neuron_pattern_id values out of range [1, %d].', size(affine_patterns, 1));
            for n = 1:N_neurons
                pat = affine_patterns(neuron_pattern_id(n), :);
                for s = 1:N_stimuli
                    data_by_stim(s, :, n) = data_by_stim(s, pat, n);
                end
            end

        otherwise
            error('apply_trial_perturbation: unknown mode "%s" (use ''none'', ''shuffle'', or ''affine'').', mode);
    end

    %% Reshape back to [N_trials_total x N_neurons], preserving row/label alignment
    data_trial_perturbed = zeros(size(data_trial));
    for s = 1:N_stimuli
        data_trial_perturbed(labels == stim_ids(s), :) = ...
            reshape(data_by_stim(s, :, :), N_trials_per_stim, N_neurons);
    end
end
