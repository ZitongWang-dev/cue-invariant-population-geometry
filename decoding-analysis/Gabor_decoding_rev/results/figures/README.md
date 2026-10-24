# Gabor_decoding_rev figures

Each folder is written by one script in `code/visualization_scripts/` (same name). Every
script reads `results/decoding_outputs/` (and, for `09`, the neuron and CNN results in
`neuron_decoding_rev/` and `ANN_decoding_rev/`), runs on its own in 2–11 s, and saves
each figure as `.png` (150 dpi) and `.fig`, plus a `.csv` of every number it plots.

Interpretation and suggested use of every figure: [FIGURE_COMMENTS.md](FIGURE_COMMENTS.md).

## Folders

| Folder | Script | What it is for |
| --- | --- | --- |
| `01_sanity_old_vs_new/` | `plot_sanity_old_vs_new.m` | Rerun vs the published pipeline (Fig 15, point 9 sweep), point by point |
| `02_paper_style_bars/` | `plot_paper_style_bars.m` | Fig 15 in the published style, corrected control added; full space and k = 5 |
| `03_paper_style_curves/` | `plot_paper_style_curves.m` | Self, PT, control, non-PT vs number of model cells |
| `04_normalized_metric_bars/` | `plot_normalized_metric_bars.m` | PT/control, PT − control, headroom and its denominator, with bootstrap CIs |
| `05_normalized_metrics_vs_size/` | `plot_normalized_metrics_vs_size.m` | The same metrics vs number of cells |
| `06_pc_sweep_accuracy/` | `plot_pc_sweep_accuracy.m` | Accuracies vs the number of PCs, k = 1–49 |
| `07_pc_sweep_metrics/` | `plot_pc_sweep_metrics.m` | Normalized metrics vs the number of PCs |
| `08_dimensionality/` | `plot_dimensionality.m` | Participation ratio and variance spectra of the models |
| `09_models_vs_neurons/` | `plot_models_vs_neurons.m` | Rebuttal point 9 redone, and neurons vs Gabor models vs CNN layers |

## Figure list

| File | Shows | Data |
| --- | --- | --- |
| **01_sanity_old_vs_new** | | |
| `old_vs_new_all_points` | Old (x) vs rerun (y) for self, non-PT, PT and control: Fig 15 and the point 9 sweep | old `Gabor_decoding` outputs + rerun |
| **02_paper_style_bars** | | |
| `fig15_gabor_cue_transfer_full` | (a) both models at 48 cells; (b) whole models (simple 192, complex 48) | gabor_cue_transfer_incremental |
| `fig15_gabor_cue_transfer_k5` | Whole models at k = 5 | gabor_pc_sweep |
| **03_paper_style_curves** | | |
| `gabor_vs_cells` | Self, PT, control, non-PT vs number of cells, one panel per model | gabor_cue_transfer_incremental |
| **04_normalized_metric_bars** | | |
| `gabor_metrics_by_design` | 48 cells, whole models, whole models at k = 5 | both |
| `gabor_metrics_per_pair_48` | Per rendering pair at 48 cells | gabor_cue_transfer_incremental |
| **05_normalized_metrics_vs_size** | | |
| `gabor_metrics_vs_cells` | PT/control, PT − control, headroom, self − control vs cells | gabor_cue_transfer_incremental |
| **06_pc_sweep_accuracy** | | |
| `pc_gabor_accuracy_vs_k` | Self, PT, control, non-PT vs k (+ full-space markers) | gabor_pc_sweep |
| **07_pc_sweep_metrics** | | |
| `pc_gabor_metrics_vs_k` | PT/control, PT − control, headroom vs k (+ full-space markers) | gabor_pc_sweep |
| **08_dimensionality** | | |
| `pr_vs_cells` | PR vs cells (raw solid, noise-corrected dashed; they coincide) | gabor_cue_transfer_incremental |
| `cumvar_vs_k` | Variance captured by the top k PCs, whole models | gabor_pc_sweep |
| `metrics_vs_pr` | PT − control and headroom vs noise-corrected PR | gabor_cue_transfer_incremental |
| **09_models_vs_neurons** | | |
| `point9_vs_cells` | Neurons vs Gabor: self, PT and control, PT − control, headroom vs population size | neuron + Gabor cue_transfer_incremental |
| `point9_vs_pr` | Neurons vs Gabor: PT, PT − control, headroom vs noise-corrected PR | |
| `k5_models_vs_neurons` | Bars at k = 5: 4 neural populations, 2 Gabor models, 8 CNN layers | all three pc_sweep folders |
| `pc_sweep_models_vs_neurons` | Self, PT − control, headroom vs k, all three families | |

## Conventions

- **Colors** follow the published Fig 15: simple model blue, complex orange; bar edges:
  dashed = non-PT, thick solid = PT (rotation only), dotted = corrected control. In `09`,
  neurons are V1 blue / V2 orange (KO solid, FR dashed), Gabor simple black and complex
  gray (diamonds), CNN layers greens (AlexNet) and purples (VGG), darker = deeper.
- **Paper-style error bars** (02, 03, 06) are the SEM over all repeat estimates, as in
  the old scripts; they understate uncertainty.
- **Normalized metrics** are ratios of pooled means (`transfer_metrics.m`) with
  percentile-bootstrap 95% CIs over cell draws, or over folds when a whole model is used once.
- **Headroom is well defined here.** The Gabor control stays at or below 0.34, so
  self − control is 0.66 or more (unlike the large CNN layers).
