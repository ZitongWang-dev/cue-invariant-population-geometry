function pop = load_ann_layer(net, layer, opts)
%{
Filename: load_ann_layer.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Loads one CNN pooling layer for the revision decoding pipeline and returns
    it in the same form as neuron_decoding_rev/load_population.m, so the shared
    decoder (pt_transfer_core.m) treats CNN units exactly like neurons.

    The trials are the ones the published analyses used (Figs 3d, 4d): for each
    unit and stimulus, 10 Poisson draws around the unit's activation at the
    centre of the feature map, made once by
    ANN_decoding/code/decoding_scripts/generate_ann_trial_data.m (rng(1)) and
    read here in place from ANN_decoding/ANN_response_data (not copied). Draws
    are independent across units, so there are no noise correlations.

    Units that never vary within one rendering (silent for all its stimuli) are
    kept, as in the published runs: z-scoring turns them into all-zero columns
    for that rendering (MATLAB's zscore leaves constant columns at 0). Their
    number per rendering is reported in pop.n_silent.

Inputs:
    net    : 'alexnet' or 'vgg' (VGG-19)
    layer  : pooling layer number (alexnet: 1, 2, 5; vgg: 1-5); MATLAB's AlexNet
             has three pooling layers, named pool1, pool2 and pool5
    opts   : optional struct
               .data_root  folder holding <net>_pooling_resps/<net>_pooling_resps_layer<L>.mat
                           (default: ../../../ANN_decoding/ANN_response_data,
                           relative to this file)

Output: struct pop with fields
    net, layer, name     e.g. 'alexnet', 5, 'AlexNet pool5'
    cue_names            {'ac','ec','ex'} (order of response_array)
    data                 1 x 3 cell of [n_stim*n_rep x N] Poisson counts
    labels               [n_stim*n_rep x 1] stimulus labels, blocked 1..n_stim
    N, n_stim, n_rep     units, stimuli (50), trials per stimulus (10)
    n_silent             [1 x 3] units constant within each rendering
    file
%}

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'data_root') || isempty(opts.data_root)
    here = fileparts(mfilename('fullpath'));
    opts.data_root = fullfile(here, '..', '..', '..', 'ANN_decoding', 'ANN_response_data');
end
net = lower(net);
assert(any(strcmp(net, {'alexnet', 'vgg'})), 'load_ann_layer: net must be ''alexnet'' or ''vgg''.');

file = fullfile(opts.data_root, [net '_pooling_resps'], sprintf('%s_pooling_resps_layer%d.mat', net, layer));
assert(exist(file, 'file') == 2, 'load_ann_layer: missing %s', file);
L = load(file, 'response_array');
ra = L.response_array;                       % 1 x 3 cell (ac, ec, ex) of {n_stim x N} cells of [n_rep x 1]
assert(numel(ra) == 3, 'load_ann_layer: expected 3 renderings in %s', file);

data = cell(1, 3);
for c = 1:3
    data{c} = cell2mat(ra{c});               % stimulus s occupies rows (s-1)*n_rep+1 .. s*n_rep
end
[n_stim, N] = size(ra{1});
n_rep = size(ra{1}{1}, 1);
assert(all(cellfun(@(x) isequal(size(x), [n_stim * n_rep, N]), data)), 'load_ann_layer: renderings differ in size.');

pop = struct();
pop.net       = net;
pop.layer     = layer;
pop.name      = sprintf('%s pool%d', net_label(net), layer);
pop.cue_names = {'ac', 'ec', 'ex'};
pop.data      = data;
pop.labels    = kron((1:n_stim)', ones(n_rep, 1));
pop.N         = N;
pop.n_stim    = n_stim;
pop.n_rep     = n_rep;
pop.n_silent  = cellfun(@(x) nnz(std(x, 0, 1) == 0), data);
pop.file      = file;
end

function s = net_label(net)
if strcmp(net, 'alexnet'), s = 'AlexNet'; else, s = 'VGG'; end
end
