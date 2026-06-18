function opts = build_affine_shift_config(monkey, vp, N_neurons, data_root)
%BUILD_AFFINE_SHIFT_CONFIG  Build opts struct for apply_trial_perturbation('affine', ...).
%
%   opts = BUILD_AFFINE_SHIFT_CONFIG(monkey, vp, N_neurons) constructs the
%   24 affine trial-permutation patterns t -> (a*p + b) mod 10 and assigns
%   one pattern to each neuron based on its session membership, matching
%   the scheme used in noise_correlation_analysis_affine.m.
%
% INPUTS
%   monkey    : 'FR' or 'KO'
%   vp        : 'V1' or 'V2'
%   N_neurons : total neurons in the pseudopopulation (for sanity check)
%   data_root : (optional) path to the neuronal_data root directory.
%               Defaults to fullfile('..','..','neuronal_data').
%
% OUTPUT
%   opts.affine_patterns    : [24 x 10] permutation table.
%   opts.neuron_pattern_id  : [N_neurons x 1] pattern index per neuron.
%   opts.neuron_session_id  : [N_neurons x 1] session index per neuron.
%   opts.session_files      : cell array of session .mat filenames used.
%
% NOTE
%   Within each recording session, neurons are assigned patterns 1, 2, 3, ...
%   sequentially so simultaneously recorded pairs always receive distinct
%   shifts. The 24 patterns are split into forward (a=1, b=0..9), reverse
%   (a=9, b=0..9), and stride-3 (a=3, b=0..3) groups.

    if nargin < 4 || isempty(data_root)
        data_root = fullfile('..','..','neuronal_data');
    end

    N_trials_per_stim = 10;

    %% Build the 24 affine patterns
    affine_patterns = zeros(24, N_trials_per_stim);
    k = 0;
    % Group 1: forward shifts (a = 1, b = 0..9)
    for b = 0:9
        k = k + 1;
        for p = 0:(N_trials_per_stim - 1)
            affine_patterns(k, p+1) = mod(1*p + b, 10) + 1;
        end
    end
    % Group 2: reverse shifts (a = 9, b = 0..9)
    for b = 0:9
        k = k + 1;
        for p = 0:(N_trials_per_stim - 1)
            affine_patterns(k, p+1) = mod(9*p + b, 10) + 1;
        end
    end
    % Group 3: stride-3 (a = 3, b = 0..3)
    for b = 0:3
        k = k + 1;
        for p = 0:(N_trials_per_stim - 1)
            affine_patterns(k, p+1) = mod(3*p + b, 10) + 1;
        end
    end
    assert(k == 24, 'Expected 24 affine patterns, got %d', k);

    %% Read session files from the 'ac' folder; bad_channel is the same across renderings
    session_dir = fullfile(data_root, monkey, vp, 'ac');
    all_files = dir(fullfile(session_dir, '*.mat'));

    % Filter out the combined population file (e.g. FR_V1_ac.mat)
    combined_pattern = sprintf('%s_%s_ac.mat', monkey, vp);
    session_files = {};
    for f = 1:numel(all_files)
        if ~strcmpi(all_files(f).name, combined_pattern)
            session_files{end+1} = all_files(f).name; %#ok<AGROW>
        end
    end
    session_files = sort(session_files);
    N_sessions = numel(session_files);
    assert(N_sessions > 0, ...
        'build_affine_shift_config: no session .mat files found in %s', session_dir);

    %% Build neuron -> session mapping and assign patterns
    neuron_session_id = [];
    neuron_pattern_id = [];
    for sess = 1:N_sessions
        sess_data = load(fullfile(session_dir, session_files{sess}), 'bad_channel');
        good_channels = setdiff(1:24, sess_data.bad_channel);
        n_good = numel(good_channels);
        assert(n_good <= 24, ...
            'Session %d has %d neurons, exceeding 24 available patterns', sess, n_good);
        neuron_session_id = [neuron_session_id; repmat(sess, n_good, 1)]; %#ok<AGROW>
        neuron_pattern_id = [neuron_pattern_id; (1:n_good)']; %#ok<AGROW>
    end

    assert(numel(neuron_session_id) == N_neurons, ...
        'Session neuron count (%d) does not match population size (%d). Check data_root/session files.', ...
        numel(neuron_session_id), N_neurons);

    %% Pack output
    opts = struct();
    opts.affine_patterns    = affine_patterns;
    opts.neuron_pattern_id  = neuron_pattern_id;
    opts.neuron_session_id  = neuron_session_id;
    opts.session_files      = session_files;
end
