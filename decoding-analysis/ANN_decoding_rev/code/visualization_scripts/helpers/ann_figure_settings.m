function S = ann_figure_settings()
%{
Filename: ann_figure_settings.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Shared settings of the ANN_decoding_rev figures: paths, the eight CNN layers,
    their colors (as the published Figs 3d and 4d: pool1 blue, pool2 orange,
    then yellow, purple, green; AlexNet's third pooling layer, pool5, is yellow
    as in the published figure), the bar-edge styles of the published figures
    (dashed = non-PT, thick solid = PT, dotted = control) and the bootstrap size.

    Also puts on the path:
      neuron_decoding_rev/code/decoding_scripts         transfer_metrics.m (shared core)
      neuron_decoding_rev/code/visualization_scripts/helpers
          plot_grouped_bars, bootstrap_metrics, to_resampling_units, mean_sem,
          save_figure, save_table (generic; they read the fields set here), and
          figure_settings / load_all_results of the neuron figures (used by
          plot_cnn_vs_neurons.m).
%}

here     = fileparts(mfilename('fullpath'));                  % .../ANN_decoding_rev/code/visualization_scripts/helpers
rev_root = fullfile(here, '..', '..', '..');                  % ANN_decoding_rev
da_root  = fullfile(rev_root, '..');                          % decoding-analysis
S.core_dir       = fullfile(da_root, 'neuron_decoding_rev', 'code', 'decoding_scripts');
S.neuron_helpers = fullfile(da_root, 'neuron_decoding_rev', 'code', 'visualization_scripts', 'helpers');
addpath(S.core_dir);
addpath(S.neuron_helpers);
S.results_root = fullfile(rev_root, 'results', 'decoding_outputs');
S.fig_root     = fullfile(rev_root, 'results', 'figures');
S.old_root     = fullfile(da_root, 'ANN_decoding', 'results', 'decoding_outputs');

% layers: {net, pooling layer}; keys name the result files
S.layers = {'alexnet',1; 'alexnet',2; 'alexnet',5; 'vgg',1; 'vgg',2; 'vgg',3; 'vgg',4; 'vgg',5};
S.keys   = arrayfun(@(i) sprintf('%s_pool%d', S.layers{i, 1}, S.layers{i, 2}), 1:size(S.layers, 1), 'UniformOutput', false);
S.nets      = {'alexnet', 'vgg'};
S.net_label = struct('alexnet', 'AlexNet', 'vgg', 'VGG-19');
S.fig3d_N   = struct('alexnet', 96, 'vgg', 64);               % published Fig 3d size (each net's pool1 size)

% colors of the published CNN figures
blue = [0 0.4470 0.7410];  orange = [0.8500 0.3250 0.0980];  yellow = [0.9290 0.6940 0.1250];
purple = [0.4940 0.1840 0.5560];  green = [0.4660 0.6740 0.1880];
S.c.layer = struct('alexnet_pool1', blue, 'alexnet_pool2', orange, 'alexnet_pool5', yellow, ...
    'vgg_pool1', blue, 'vgg_pool2', orange, 'vgg_pool3', yellow, 'vgg_pool4', purple, 'vgg_pool5', green);
S.ls.net  = struct('alexnet', '-', 'vgg', '--');             % network told apart by line style when both share an axis
S.mk.net  = struct('alexnet', 'o', 'vgg', 's');
S.c.dir   = [blue; orange];                                  % hierarchy: pool2-to-pool5, pool5-to-pool2
S.c.V1    = blue;  S.c.V2 = orange;                          % neurons (plot_cnn_vs_neurons)
S.c.err   = [1 0 0];
S.edge.self = {'-',  0.5};
S.edge.nopt = {'--', 1.5};
S.edge.pt   = {'-',  1.5};
S.edge.ctrl = {':',  1.5};

S.k_main = 5;
S.chance = 0.02;
S.n_boot = 2000;
end
