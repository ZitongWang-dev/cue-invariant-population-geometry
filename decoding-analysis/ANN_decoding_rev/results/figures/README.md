# ANN_decoding_rev figures

Each folder is written by one script in `code/visualization_scripts/` (same name). Every
script reads `results/decoding_outputs/` (and, for `09`, the neuron results in
`neuron_decoding_rev/results/decoding_outputs/`), runs on its own in 5–20 s, and saves
each figure as `.png` (150 dpi) and `.fig`, plus a `.csv` of every number it plots.

Interpretation and suggested use of every figure: [FIGURE_COMMENTS.md](FIGURE_COMMENTS.md).

## Folders

| Folder | Script | What it is for |
| --- | --- | --- |
| `01_sanity_old_vs_new/` | `plot_sanity_old_vs_new.m` | Rerun vs the published pipeline (Figs 3d, 4d, hierarchy), point by point |
| `02_paper_style_bars/` | `plot_paper_style_bars.m` | Fig 3d and the hierarchy bars in the published style, corrected control added; full space and k = 5 |
| `03_paper_style_curves/` | `plot_paper_style_curves.m` | Fig 4d in the published style, corrected control added |
| `04_normalized_metric_bars/` | `plot_normalized_metric_bars.m` | PT/control, PT − control, headroom and its denominator, with bootstrap CIs |
| `05_normalized_metrics_vs_size/` | `plot_normalized_metrics_vs_size.m` | The same metrics vs number of units (the "plateau" test) |
| `06_pc_sweep_accuracy/` | `plot_pc_sweep_accuracy.m` | Accuracies vs the number of PCs, k = 1–49 |
| `07_pc_sweep_metrics/` | `plot_pc_sweep_metrics.m` | Normalized metrics vs the number of PCs |
| `08_dimensionality/` | `plot_dimensionality.m` | Participation ratio and variance spectra of the CNN layers |
| `09_cnn_vs_neurons/` | `plot_cnn_vs_neurons.m` | CNN layers and recorded neurons on the same axes |

## Figure list

| File | Shows | Data |
| --- | --- | --- |
| **01_sanity_old_vs_new** | | |
| `old_vs_new_all_points` | Old (x) vs rerun (y) for self, non-PT, PT and control: Fig 3d, Fig 4d at shared unit counts, hierarchy | old `ANN_decoding` outputs + rerun |
| **02_paper_style_bars** | | |
| `fig3d_cnn_cue_transfer_full` / `_k5` | Fig 3d: AlexNet (top), VGG-19 (bottom), per rendering pair; 96 / 64 units, or k = 5 on whole layers | ann_cue_transfer_incremental, ann_pc_sweep_cue |
| `hierarchy_by_training_full` / `_k5` | AlexNet pool2 ↔ pool5 per training rendering (Fig 4b layout) | ann_hierarchy, ann_pc_sweep_hierarchy |
| `hierarchy_all_combos_full` / `_k5` | The same, all 9 train × test combinations (Fig 12 layout) | |
| **03_paper_style_curves** | | |
| `fig4d_cnn_vs_units` | Fig 4d: self, PT and control vs units, per layer | ann_cue_transfer_incremental |
| `fig4d_cnn_vs_units_per_layer` | One panel per layer, with non-PT | |
| **04_normalized_metric_bars** | | |
| `cue_metrics_full_vs_k5` | Whole layers: full space vs k = 5 | ann_cue_transfer_incremental (whole layer), ann_pc_sweep_cue |
| `cue_metrics_fig3d_size_per_pair` | Per rendering pair at the Fig 3d size | |
| `hierarchy_metrics_full_vs_k5` | Per direction and pooled | |
| **05_normalized_metrics_vs_size** | | |
| `cue_metrics_vs_units` | PT/control, PT − control, headroom, self − control vs units | ann_cue_transfer_incremental |
| **06_pc_sweep_accuracy** | | |
| `pc_cue_accuracy_vs_k`, `pc_hierarchy_accuracy_vs_k` | Self, PT, control, non-PT vs k (+ full-space markers) | ann_pc_sweep_* |
| **07_pc_sweep_metrics** | | |
| `pc_cue_metrics_vs_k`, `pc_hierarchy_metrics_vs_k` | PT/control, PT − control, headroom vs k (+ full-space markers) | ann_pc_sweep_* |
| **08_dimensionality** | | |
| `pr_vs_units` | PR vs units (raw solid, noise-corrected dashed) | ann_cue_transfer_incremental |
| `cumvar_vs_k` | Variance captured by the top k PCs, whole layers | ann_pc_sweep_cue |
| `metrics_vs_pr` | PT − control and headroom vs noise-corrected PR along the unit sweep | ann_cue_transfer_incremental |
| **09_cnn_vs_neurons** | | |
| `k5_cnn_vs_neurons` | Self, PT − control, headroom, PT/control at k = 5, 4 neural populations + 8 layers | both pc_sweep_cue folders |
| `pc_sweep_cnn_vs_neurons` | Self, PT − control, headroom vs k | |
| `vs_size_cnn_vs_neurons` | The same vs number of neurons / units, full space | both cue_transfer_incremental folders |
| `vs_pr_cnn_vs_neurons` | PT − control and headroom vs noise-corrected PR, full space | |

## Conventions

- **Colors** follow the published CNN figures: pool1 blue, pool2 orange, then yellow,
  purple, green (AlexNet's third pooling layer, pool5, is yellow). Bar edges: dashed =
  non-PT, thick solid = PT (rotation only), dotted = corrected control. Where AlexNet
  and VGG share an axis: AlexNet solid / circles, VGG dashed / squares. In `09`,
  neurons are V1 blue / V2 orange and the CNN layers greens (AlexNet) and purples (VGG).
- **Paper-style error bars** (02, 03, 06) are the SEM over all repeat estimates, as in
  the old scripts; they understate uncertainty.
- **Normalized metrics** are ratios of pooled means (`transfer_metrics.m`) with
  percentile-bootstrap 95% CIs over unit draws, or over folds when a whole layer is used
  once (and for the hierarchy, whose five draws differ by one unit).
- **Headroom near ceiling.** CNN self-decoding is close to 1 and the corrected control
  is high in full space, so headroom's denominator (self − control) can fall below 0.1
  and headroom then swings with small changes in PT or control. Figures 04 and 05 show
  the denominator; PT − control is the more stable full-space measure.
