function S = gabor_figure_settings()
%{
Filename: gabor_figure_settings.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Shared settings of the Gabor_decoding_rev figures: paths, the two models
    (simple, complex), their colors as the published Fig 15 (simple blue,
    complex orange), the bar-edge styles of the published figures (dashed =
    non-PT, thick solid = PT, dotted = control) and the bootstrap size.

    Also puts on the path:
      neuron_decoding_rev/code/decoding_scripts          transfer_metrics.m (shared core)
      neuron_decoding_rev/code/visualization_scripts/helpers
          plot_grouped_bars, bootstrap_metrics, to_resampling_units, mean_sem,
          save_figure, save_table (generic; they read the fields set here), and
          figure_settings / load_all_results of the neuron figures
      ANN_decoding_rev/code/visualization_scripts/helpers
          ann_figure_settings / load_ann_results (used by plot_models_vs_neurons.m)
%}

here     = fileparts(mfilename('fullpath'));                  % .../Gabor_decoding_rev/code/visualization_scripts/helpers
rev_root = fullfile(here, '..', '..', '..');                  % Gabor_decoding_rev
da_root  = fullfile(rev_root, '..');                          % decoding-analysis
S.core_dir       = fullfile(da_root, 'neuron_decoding_rev', 'code', 'decoding_scripts');
S.neuron_helpers = fullfile(da_root, 'neuron_decoding_rev', 'code', 'visualization_scripts', 'helpers');
S.ann_helpers    = fullfile(da_root, 'ANN_decoding_rev', 'code', 'visualization_scripts', 'helpers');
addpath(S.core_dir);
addpath(S.neuron_helpers);
addpath(S.ann_helpers);
S.results_root = fullfile(rev_root, 'results', 'decoding_outputs');
S.fig_root     = fullfile(rev_root, 'results', 'figures');
S.old_root     = fullfile(da_root, 'Gabor_decoding', 'results', 'decoding_outputs');

S.models    = {'simple', 'complex'};
S.label     = struct('simple', 'Gabor simple', 'complex', 'Gabor complex');
S.fig15a_N  = 48;                                             % published Fig 15a size (complex-cell count)

% colors and styles of the published figures
S.c.model = struct('simple', [0 0.4470 0.7410], 'complex', [0.8500 0.3250 0.0980]);
S.c.err   = [1 0 0];
S.edge.self = {'-',  0.5};
S.edge.nopt = {'--', 1.5};
S.edge.pt   = {'-',  1.5};
S.edge.ctrl = {':',  1.5};

S.k_main = 5;
S.chance = 0.02;
S.n_boot = 2000;
end
