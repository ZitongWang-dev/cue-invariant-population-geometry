function S = figure_settings()
%{
Filename: figure_settings.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Shared settings of the neuron_decoding_rev figures: paths, the colors and
    bar-edge styles of the published figures (V1 blue, V2 orange; condition told
    apart by the bar edge: dashed = non-PT, thick solid = PT, dotted = control),
    the anchor sizes of the main figures, and the bootstrap size. Also puts
    code/decoding_scripts on the path (transfer_metrics.m).
%}

here = fileparts(mfilename('fullpath'));                  % .../code/visualization_scripts/helpers
addpath(fullfile(here, '..', '..', 'decoding_scripts'));
S.results_root = fullfile(here, '..', '..', '..', 'results', 'decoding_outputs');
S.fig_root     = fullfile(here, '..', '..', '..', 'results', 'figures');
S.old_root     = fullfile(here, '..', '..', '..', '..', 'neuronal_decoding', 'results', 'decoding_outputs');

% colors and styles of the published figures
S.c.V1    = [0 0.4470 0.7410];
S.c.V2    = [0.8500 0.3250 0.0980];
S.c.err   = [1 0 0];
S.c.pairs = [0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560; 0.9290 0.6940 0.1250];   % EC-AC, EC-EX, AC-EX
S.c.cues  = [0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560; 0.3010 0.7450 0.9330];   % AC, EC, EX
S.edge.self = {'-',  0.5};
S.edge.nopt = {'--', 1.5};
S.edge.pt   = {'-',  1.5};
S.edge.ctrl = {':',  1.5};

S.fig1 = [100 100 1080 540];      % one panel, as the published figures
S.fig2 = [100 100 1250 900];      % KO (top) and FR (bottom), as the paper panels

S.monkeys  = {'KO', 'FR'};        % panel order of the paper: KO top, FR bottom
S.pops     = {'KO_V1', 'KO_V2', 'FR_V1', 'FR_V2'};
S.anchor_N = struct('FR', 48, 'KO', 109);    % Fig 3b: V1 size; V2 subsampled to it
S.anchor_h = struct('FR', 24, 'KO', 54);     % Fig 4a: half size of the paper
S.k_main   = 5;
S.chance   = 0.02;
S.n_boot   = 2000;
end
