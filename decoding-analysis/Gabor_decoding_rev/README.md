# Gabor_decoding_rev — Gabor-model decoding with the corrected control

Revision reruns of the Gabor-filter cell-model decoding of the J Neurosci resubmission:
Fig 15 (simple vs complex model cells, at 48 cells and with the whole models) and the
Gabor side of rebuttal point 9 (neurons vs filter bank vs population size and
participation ratio), plus a PC-k version for comparing the models with neurons and CNN
layers at the same rotation freedom. Every analysis scores Procrustes transfer (PT)
against the **corrected control** (one permutation of the stimuli applied to both the
source class means and the held-out source trials). The old scripts in
`Gabor_decoding/` used the control that stays at chance by construction, including
the Gabor sweep behind rebuttal point 9.

**Shared core.** The drivers here call the decoder of `neuron_decoding_rev`
(`pt_transfer_core.m`, plus `population_pr.m`, `transfer_metrics.m`,
`progress_queue.m`) through `addpath`, so model cells, CNN units and neurons go
through exactly the same fold design, control, PC projection and participation ratio.
See `../neuron_decoding_rev/README.md` for how one comparison is decoded.

**Data** are read in place from `../Gabor_decoding/Gabor_filters_response_data/sampled_trials`
(not copied): the saved trials of the published analyses, i.e. 10 Poisson draws per
cell and stimulus around the filter response (`Gabor_decoding/code/decoding_scripts/generate_gabor_trial_data.m`,
`rng(1)`).

## Layout

```
Gabor_decoding_rev/
  code/
    decoding_scripts/
      load_gabor_model.m                    one model ('simple' or 'complex'): [500 x N] counts per rendering + labels
      run_gabor_cue_transfer_incremental.m  Fig 15a/b + rebuttal point 9: cue transfer vs number of cells
      run_gabor_pc_sweep.m                  k = 1..49 sweep of cue transfer, whole models
      run_all_gabor.m                       MASTER: runs both drivers back to back, unattended, with a log
      run_gabor_driver.m                    runs one driver in its own workspace (used by run_all_gabor)
    visualization_scripts/                  one script per figure series; each writes results/figures/<NN_series>/
      plot_sanity_old_vs_new.m                  01  rerun vs published pipeline (Fig 15, point 9 sweep)
      plot_paper_style_bars.m                   02  Fig 15, published style + corrected control (full, k = 5)
      plot_paper_style_curves.m                 03  self, PT, control, non-PT vs number of cells
      plot_normalized_metric_bars.m             04  PT/control, PT - control, headroom and its denominator
      plot_normalized_metrics_vs_size.m         05  the same metrics vs number of cells
      plot_pc_sweep_accuracy.m                  06  accuracies vs number of PCs
      plot_pc_sweep_metrics.m                   07  normalized metrics vs number of PCs
      plot_dimensionality.m                     08  participation ratio and variance spectra
      plot_models_vs_neurons.m                  09  rebuttal point 9 redone; neurons vs Gabor vs CNN layers
      helpers/                                  gabor_figure_settings, load_gabor_results; the generic helpers
                                                are shared from neuron_decoding_rev (and the CNN loader from
                                                ANN_decoding_rev)
  results/
    decoding_outputs/<analysis>/<model>.mat
    logs/run_all_gabor_<set>_<date>.log
    figures/<NN_series>/                    figures (.png, .fig) + a .csv of every plotted number;
                                            index in results/figures/README.md, comments in FIGURE_COMMENTS.md
```

Each result file holds one variable `res` with the same fields as the matching neuron
driver (`cue_transfer_incremental`, `pc_sweep_cue`), plus `model` and `n_silent`; the
driver headers list every field.

## Models and design

| Model | Cells | Cell grid (cue transfer vs cells) |
| --- | --- | --- |
| simple (`even and odd combined`) | 192: 3 scales × 16 orientations × even/odd × ON/OFF (half-wave rectified) | 10, 20, …, 190, 192, and 48 |
| complex | 48: energy of the even and odd filters | 10, 20, 30, 40, 48 |

- **Fig 15a** is the 48-cell point of both models (the simple model down-sampled to
  the complex-cell count, as published); **Fig 15b** is the whole models; the 10-cell
  grid matches the neural sweep and the published rebuttal sweep exactly.
- 20 random cell subsets per size (1 for the whole model), shared by the six ordered
  rendering pairs; 10 folds; 10 control permutations per fold; participation ratio
  (raw and noise-corrected, 20 half-splits) at every point.
- **PC sweep**: whole models, k = 1–10, 15–45 and 49 (48 for the complex model), each side
  projected onto its own top-k class-mean PCs; non-PT projects the source onto the
  target basis.

**How the model cells differ from the neurons.** Poisson draws are independent across
cells, so there are no noise correlations, and the counts are large (means of about
300–4,000 per trial), so the trials are nearly noiseless: self-decoding is 1.000 at
every size. Some rectified simple-cell channels are silent for a whole rendering (16 in
AC, 9 in EX); they are kept, as published, and become all-zero columns after z-scoring
(`res.n_silent`). A noise-matched variant (rates rescaled so self-decoding matches the
neurons) is deferred until the Poisson results are in, as for the CNNs.

## Running

From MATLAB (any current folder):

```matlab
run_all_gabor('test')    % ~6 min, tiny configs in a temp folder: run once on the machine first
run_all_gabor('all')     % both drivers, ~5 h on one 8-worker machine
```

A job that fails is logged and skipped, and calling `run_all_gabor` again resumes every
driver where it stopped. **One driver at a time:** both drivers run both models by
default; to change a setting, set an override first, e.g.

```matlab
cfg_override = struct('models', {{'complex'}});
run_gabor_cue_transfer_incremental
clear cfg_override    % before running another driver
```

Seeds are fixed per (model, size, draw, item), so results do not depend on which
machine runs which model.

## Run time

One `pt_transfer_core` call takes about 45–55 s on a worker of the laptop's 8-worker
pool at these sizes (2026-10-15 test), about 600–650 calls per hour.

| Driver | Calls | Time, one machine |
| --- | --- | --- |
| run_gabor_cue_transfer_incremental | 2,892 (simple 2,406, complex 486) | ~4.5 h |
| run_gabor_pc_sweep | 216 (18 k values × 6 pairs per model) | ~0.3 h |
| **Total** | **3,108** | **~5 h** |

The test run reproduces the published Fig 15b values exactly (simple, 192 cells: self
1.000, non-PT 0.062, PT 0.397; complex, 48 cells: self 1.000, non-PT 0.038, PT 0.139).
