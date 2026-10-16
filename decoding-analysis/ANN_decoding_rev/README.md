# ANN_decoding_rev — CNN decoding with the corrected control

Revision reruns of the CNN (AlexNet, VGG-19) decoding of the J Neurosci resubmission:
Fig 3d (cue transfer per pooling layer), Fig 4d (cue transfer vs number of units), the
AlexNet pool2 ↔ pool5 "along the hierarchy" analysis cited in the rebuttal (point 1),
and a PC-k version of both for comparing CNN layers with neurons at the same rotation
freedom. Every analysis scores Procrustes transfer (PT) against the **corrected
control** (one permutation of the stimuli applied to both the source class means and
the held-out source trials). The old scripts in `ANN_decoding/` used the control that
stays at chance by construction.

**Shared core.** The drivers here call the decoder of `neuron_decoding_rev`
(`pt_transfer_core.m`, plus `population_pr.m`, `transfer_metrics.m`,
`progress_queue.m`) through `addpath`, so CNN units and neurons go through exactly the
same fold design, control, PC projection and participation ratio. See
`../neuron_decoding_rev/README.md` for how one comparison is decoded.

**Data** are read in place from `../ANN_decoding/ANN_response_data` (not copied): the
saved trials of the published analyses, i.e. 10 Poisson draws per unit and stimulus
around the unit's activation at the centre of the feature map
(`ANN_decoding/code/decoding_scripts/generate_ann_trial_data.m`, `rng(1)`).

## Layout

```
ANN_decoding_rev/
  code/
    decoding_scripts/
      load_ann_layer.m                  one CNN layer: [500 x N] Poisson counts per rendering + labels
      run_ann_cue_transfer_incremental.m Fig 3d (one point) + Fig 4d: cue transfer vs number of units
      run_ann_hierarchy.m               AlexNet pool2 <-> pool5, both directions, 3 x 3 renderings
      run_ann_pc_sweep.m                k = 1..49 sweep for cue transfer (8 layers) and the hierarchy
      run_all_ann.m                     MASTER: runs the drivers back to back, unattended, with a log
      run_ann_driver.m                  runs one driver in its own workspace (used by run_all_ann)
    visualization_scripts/              one script per figure series; each writes results/figures/<NN_series>/
      plot_sanity_old_vs_new.m              01  rerun vs published pipeline (Figs 3d, 4d, hierarchy)
      plot_paper_style_bars.m               02  Fig 3d and hierarchy bars, published style + corrected control (full, k = 5)
      plot_paper_style_curves.m             03  Fig 4d, published style + corrected control
      plot_normalized_metric_bars.m         04  PT/control, PT - control, headroom and its denominator
      plot_normalized_metrics_vs_size.m     05  the same metrics vs number of units
      plot_pc_sweep_accuracy.m              06  self, PT, control, non-PT vs number of PCs
      plot_pc_sweep_metrics.m               07  normalized metrics vs number of PCs
      plot_dimensionality.m                 08  participation ratio and variance spectra
      plot_cnn_vs_neurons.m                 09  CNN layers and recorded neurons on the same axes
      helpers/                              ann_figure_settings, load_ann_results; the generic helpers
                                            (plot_grouped_bars, bootstrap_metrics, ...) are shared from
                                            neuron_decoding_rev/code/visualization_scripts/helpers
  results/
    decoding_outputs/<analysis>/<net>_pool<L>.mat   (hierarchy: <net>_pool<lo>_pool<hi>.mat)
    logs/run_all_ann_<set>_<date>.log
    figures/<NN_series>/                figures (.png, .fig) + a .csv of every plotted number;
                                        index in results/figures/README.md, comments in FIGURE_COMMENTS.md
```

Each result file holds one variable `res` with the same fields as the matching neuron
driver (`cue_transfer_incremental`, `across_area`, `pc_sweep_cue`, `pc_sweep_across`),
plus `net`, `layer` and `n_silent`; the driver headers list every field.

## Layers and design

| Network | Layer | Units | Unit grid (cue transfer vs units) |
| --- | --- | --- | --- |
| AlexNet | pool1 | 96 | 10, 20, 40, 64, 80, 96 |
| AlexNet | pool2 | 250 | 10 … 96, 160, 250 |
| AlexNet | pool5 | 249 | 10 … 96, 160, 249 |
| VGG-19 | pool1 | 64 | 10, 20, 40, 64 |
| VGG-19 | pool2 | 128 | 10 … 96, 128 |
| VGG-19 | pool3 | 256 | 10 … 96, 160, 256 |
| VGG-19 | pool4 | 511 | 10 … 96, 160, 320, 511 |
| VGG-19 | pool5 | 507 | 10 … 96, 160, 320, 507 |

- **Fig 3d** is the 96-unit point of every AlexNet layer and the 64-unit point of every
  VGG layer (each network's pool1 size, as published); **Fig 4d** is the whole curve.
  MATLAB's AlexNet has three pooling layers, named pool1, pool2 and pool5; the published
  Fig 4d legend calls the last one "Pool 3".
- 20 random unit subsets per size (1 for the whole layer), shared by the six ordered
  rendering pairs; 10 folds; 10 control permutations per fold; participation ratio
  (raw and noise-corrected, 20 half-splits) at every point.
- **Hierarchy**: pool2 (250 units) is subsampled to the size of pool5 (249), 5 draws.
  "pool2-to-pool5" = decoder trained on pool2, tested on pool5 (as "V1-to-V2" across
  areas). Non-PT maps different units onto each other and is only a sanity check.
- **PC sweep**: whole layers, k = 1–10, 15–45, 49, each side projected onto its own top-k
  class-mean PCs; non-PT (cue transfer) projects the source onto the target basis.

**How CNN populations differ from the neurons.** Poisson draws are independent across
units, so there are no noise correlations, and the noise level is set by the activation
scale, which differs by more than 100× across layers (VGG pool3 averages about 600,
AlexNet pool5 about 2). Self-decoding is 0.90–1.00 at the Fig 3d sizes, against
0.2–0.5 for the neurons. Many deep-layer units are silent within one rendering (e.g.
127 of 249 in AlexNet pool5, 278 of 507 in VGG pool5); they are kept, as published,
and become all-zero columns after z-scoring (`res.n_silent` gives the counts).

## Running

**Unattended (recommended)**, from MATLAB (any current folder):

```matlab
run_all_ann('test')      % ~6 min, tiny configs in a temp folder: run once per machine first
run_all_ann('all')       % every driver, ~11 h on one 8-worker machine
run_all_ann('desktop')   % or split across two machines: ~6.5 h here ...
run_all_ann('laptop')    % ... and ~5 h here
```

A job that fails is logged and skipped, and calling `run_all_ann` again with the same set
resumes every driver where it stopped. **One driver at a time:** every driver runs all
its layers by default; to change a setting, set an override first, e.g.

```matlab
cfg_override = struct('layers', {{'vgg',4; 'vgg',5}});
run_ann_cue_transfer_incremental
clear cfg_override    % before running another driver
```

Seeds are fixed per (layer, size, draw, item), so results do not depend on which machine
runs which layer, and the two machines never write the same file.

## Run time (measured 2026-10-12, R2025b, single-threaded worker)

One `pt_transfer_core` call takes about 28 s at 96 units and 38 s at 511 units in full
space, and about 24 s in PC mode, on one single-threaded worker of the laptop; with 8
workers in parallel expect roughly 600–700 calls per hour.

| Driver | Calls | Time, one machine |
| --- | --- | --- |
| run_ann_cue_transfer_incremental | 6,168 (AlexNet 606 / 846 / 846; VGG 366 / 726 / 846 / 966 / 966) | ~9.5 h |
| run_ann_hierarchy | 90 (5 draws × 9 combinations × 2 directions) | ~0.2 h |
| run_ann_pc_sweep, cue | 864 (18 k values × 6 pairs per layer) | ~1.2 h |
| run_ann_pc_sweep, hierarchy | 324 (18 × 9 × 2) | ~0.5 h |
| **Total** | **7,446** | **~11.4 h** |

The `run_all_ann` split gives the faster machine ('desktop') the four large layers
(AlexNet pool5, VGG pool3–5), the hierarchy and its PC sweep, and the other ('laptop')
the four small layers and the cue PC sweep of all eight layers.
