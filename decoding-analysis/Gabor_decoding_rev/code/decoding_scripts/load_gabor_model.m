function pop = load_gabor_model(model, opts)
%{
Filename: load_gabor_model.m
Author:   Zitong Wang
Date:     2026-10-15

Description:
    Loads one Gabor-filter cell model for the revision decoding pipeline and
    returns it in the same form as neuron_decoding_rev/load_population.m, so the
    shared decoder (pt_transfer_core.m) treats model cells exactly like neurons.

      'simple'   192 model simple cells: even and odd Gabor filters, 3 scales x
                 16 orientations, each half-wave rectified into an ON and an OFF
                 channel (Methods; Fig 15). Stored as 'even and odd combined'.
      'complex'  48 model complex cells: the energy (sum of squared even and odd
                 responses) of the same filters.

    The trials are the ones the published analyses used (Fig 15, rebuttal point
    9): 10 Poisson draws per cell and stimulus around the filter response, made
    once by Gabor_decoding/code/decoding_scripts/generate_gabor_trial_data.m
    (rng(1)) and read here in place from
    Gabor_decoding/Gabor_filters_response_data/sampled_trials (not copied).
    Draws are independent across cells, so there are no noise correlations, and
    the counts are large (hundreds to thousands per trial), so the trials are
    nearly noiseless.

    Cells that never vary within one rendering (some rectified simple-cell
    channels are silent for a whole rendering) are kept, as in the published
    runs: z-scoring turns them into all-zero columns for that rendering. Their
    number per rendering is reported in pop.n_silent.

Inputs:
    model  : 'simple' or 'complex'
    opts   : optional struct
               .data_root  folder holding <stored name>/<stored name>.mat
                           (default: ../../../Gabor_decoding/Gabor_filters_response_data/sampled_trials,
                           relative to this file)

Output: struct pop with fields
    model, name          e.g. 'simple', 'Gabor simple'
    cue_names            {'ac','ec','ex'} (order of response_array)
    data                 1 x 3 cell of [n_stim*n_rep x N] Poisson counts
    labels               [n_stim*n_rep x 1] stimulus labels, blocked 1..n_stim
    N, n_stim, n_rep     cells, stimuli (50), trials per stimulus (10)
    n_silent             [1 x 3] cells constant within each rendering
    file
%}

if nargin < 2 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'data_root') || isempty(opts.data_root)
    here = fileparts(mfilename('fullpath'));
    opts.data_root = fullfile(here, '..', '..', '..', 'Gabor_decoding', 'Gabor_filters_response_data', 'sampled_trials');
end
model = lower(model);
switch model
    case 'simple',  stored = 'even and odd combined';
    case 'complex', stored = 'complex';
    otherwise, error('load_gabor_model: model must be ''simple'' or ''complex''.');
end

file = fullfile(opts.data_root, stored, [stored '.mat']);
assert(exist(file, 'file') == 2, 'load_gabor_model: missing %s', file);
L = load(file, 'response_array');
ra = L.response_array;                       % 1 x 3 cell (ac, ec, ex) of {n_stim x N} cells of [n_rep x 1]
assert(numel(ra) == 3, 'load_gabor_model: expected 3 renderings in %s', file);

data = cell(1, 3);
for c = 1:3
    data{c} = cell2mat(ra{c});               % stimulus s occupies rows (s-1)*n_rep+1 .. s*n_rep
end
[n_stim, N] = size(ra{1});
n_rep = size(ra{1}{1}, 1);
assert(all(cellfun(@(x) isequal(size(x), [n_stim * n_rep, N]), data)), 'load_gabor_model: renderings differ in size.');

pop = struct();
pop.model     = model;
pop.name      = ['Gabor ' model];
pop.cue_names = {'ac', 'ec', 'ex'};
pop.data      = data;
pop.labels    = kron((1:n_stim)', ones(n_rep, 1));
pop.N         = N;
pop.n_stim    = n_stim;
pop.n_rep     = n_rep;
pop.n_silent  = cellfun(@(x) nnz(std(x, 0, 1) == 0), data);
pop.file      = file;
end
