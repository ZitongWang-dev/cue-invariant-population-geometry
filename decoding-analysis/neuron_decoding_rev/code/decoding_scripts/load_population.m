function pop = load_population(monkey, area, opts)
%{
Filename: load_population.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Loads one recorded population (monkey x area) for the revision decoding
    pipeline and returns the spike counts in the decoding window as one
    [trials x neurons] matrix per rendering (cue). The data are read in place
    from neuronal_decoding/neuronal_data (not copied into this folder).

    Row order is the one multiclass_svmloader_PT produces: 50 stimulus blocks of
    10 repeats, and row r of a block is the SAME stored repeat for every neuron.
    Simultaneously recorded (same-session) neurons therefore keep their real
    trial-to-trial noise correlations; neurons from different sessions are
    paired arbitrarily. Nothing in this pipeline reorders trials per neuron.

Inputs:
    monkey : 'FR' or 'KO'
    area   : 'V1' or 'V2'
    opts   : optional struct
               .timewindow  [start end] samples of the stored 1 ms trace
                            (default [330 630] = 30-330 ms after stimulus onset,
                            as in every published decoding script)
               .data_root   folder holding <monkey>/<area>/<monkey>_<area>_allstim.mat
                            (default: ../../../neuronal_decoding/neuronal_data,
                            relative to this file)

Output: struct pop with fields
    monkey, area, name   e.g. 'KO', 'V1', 'KO V1'
    cue_names            {'ac','ec','ex'} (order of three_stim_array)
    data                 1 x 3 cell of [n_stim*n_rep x N] spike counts
    labels               [n_stim*n_rep x 1] stimulus labels, blocked 1..n_stim
    N, n_stim, n_rep     neurons, stimuli (50), repeats per stimulus (10)
    timewindow, file
%}

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'timewindow') || isempty(opts.timewindow), opts.timewindow = [330 630]; end
if ~isfield(opts, 'data_root') || isempty(opts.data_root)
    here = fileparts(mfilename('fullpath'));
    opts.data_root = fullfile(here, '..', '..', '..', 'neuronal_decoding', 'neuronal_data');
end

file = fullfile(opts.data_root, monkey, area, sprintf('%s_%s_allstim.mat', monkey, area));
assert(exist(file, 'file') == 2, 'load_population: data file not found: %s', file);
tmp = load(file, 'three_stim_array');
arr = tmp.three_stim_array;

cue_names = {'ac', 'ec', 'ex'};
assert(numel(arr) == numel(cue_names), 'load_population: expected 3 renderings in %s', file);

data = cell(1, numel(cue_names));
labels = [];
for c = 1:numel(cue_names)
    [data{c}, lab] = multiclass_svmloader_PT(arr{c}, opts.timewindow);
    if c == 1
        labels = lab;
    else
        assert(isequal(lab, labels), 'load_population: label order differs between renderings.');
        assert(size(data{c}, 2) == size(data{1}, 2), 'load_population: neuron count differs between renderings.');
    end
end

n_stim = max(labels);
n_rep  = numel(labels) / n_stim;
assert(n_rep == round(n_rep) && isequal(labels, kron((1:n_stim)', ones(n_rep, 1))), ...
    'load_population: labels must come in equal, consecutive stimulus blocks.');

pop = struct('monkey', monkey, 'area', area, 'name', sprintf('%s %s', monkey, area), ...
    'cue_names', {cue_names}, 'data', {data}, 'labels', labels, ...
    'N', size(data{1}, 2), 'n_stim', n_stim, 'n_rep', n_rep, ...
    'timewindow', opts.timewindow, 'file', file);
end
