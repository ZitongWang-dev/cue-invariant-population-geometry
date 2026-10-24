# neuron_decoding_rev figures

Each folder is written by one script in `code/visualization_scripts/` (same number, same
name). Every script reads `results/decoding_outputs/`, runs on its own in a few seconds to
half a minute (no decoding), and saves each figure as `.png` (150 dpi) and `.fig`, plus a
`.csv` of every number it plots.

Interpretation and suggested use of every figure: [FIGURE_COMMENTS.md](FIGURE_COMMENTS.md).

## Folders

| Folder | Script | What it is for |
| --- | --- | --- |
| `01_sanity_old_vs_new/` | `plot_sanity_old_vs_new.m` | Rerun vs the published pipeline, bar by bar |
| `02_paper_style_bars/` | `plot_paper_style_bars.m` | Bar panels of the main figures in the published style, corrected control added; full space and k = 5 |
| `03_paper_style_curves/` | `plot_paper_style_curves.m` | Incremental curves (Figs 4c, 13) in the published style, corrected control added |
| `04_normalized_metric_bars/` | `plot_normalized_metric_bars.m` | PT/control, headroom and PT − control as bars, with bootstrap CIs |
| `05_normalized_metrics_vs_size/` | `plot_normalized_metrics_vs_size.m` | The same metrics vs population / half size |
| `06_pc_sweep_accuracy/` | `plot_pc_sweep_accuracy.m` | Accuracies vs the number of PCs, k = 1–49 |
| `07_pc_sweep_metrics/` | `plot_pc_sweep_metrics.m` | Normalized metrics vs the number of PCs |
| `08_dimensionality/` | `plot_dimensionality.m` | Participation ratio (raw and noise-corrected) and variance spectra |
| `09_nc_injection/` | `plot_nc_injection.m` | Noise correlations injected into the pseudo-populations (rebuttal Point 8) |

## Figure list

| File | Shows | Data |
| --- | --- | --- |
| **01_sanity_old_vs_new** | | |
| `old_vs_new_all_bars` | Every bar of old Figs 3b, 4a, 4b/12: old (x) vs rerun (y) for self, non-PT, PT and control | old `neuronal_decoding` outputs + rerun |
| `old_vs_new_bars.csv` | Every plotted point | |
| **02_paper_style_bars** | | |
| `fig3b_cue_transfer_full` | Cue transfer per pair, V1 vs V2 (V1 whole, V2 subsampled to the V1 count) | cue_transfer_incremental at N = 48 / 109 |
| `fig4a_halves_full` | Within-area halves per rendering, halves of 24 (FR) / 54 (KO) | halves_incremental |
| `fig4b_across_area_full` | Across areas per training rendering, averaged over test renderings | across_area |
| `fig12_across_area_all_combos_full` | Across areas, all 9 train × test renderings | across_area |
| `*_k5` (same four names) | The same panels at k = 5, whole populations (halves of floor(N/2)); no non-PT where the neurons differ | pc_sweep_* |
| `paper_style_bars.csv` | Every bar (mean, SEM) | |
| **03_paper_style_curves** | | |
| `fig4c_cue_transfer_vs_N` | Self, PT, control, non-PT vs number of neurons | cue_transfer_incremental |
| `fig13_halves_vs_size` | The same vs half size | halves_incremental |
| `paper_style_curves.csv` | Every point (mean, SEM) | |
| **04_normalized_metric_bars** | | |
| `cue_metrics_full_vs_k5` | Cue transfer, whole populations: full space vs k = 5 | cue (N = all), pc_sweep_cue |
| `cue_metrics_per_pair_full` / `_k5` | The same per rendering pair (full space at the Fig 3b size) | |
| `halves_metrics_full_vs_k5` | Halves of floor(N/2): full space vs k = 5 | |
| `across_metrics_full_vs_k5` | Across areas per direction and pooled | |
| `cue_vs_halves_matched_size` | Cue transfer vs halves at the same neuron count (FR 24, KO 54) | |
| `normalized_metric_bars.csv` | Every bar with its CI | |
| **05_normalized_metrics_vs_size** | | |
| `cue_metrics_vs_N`, `halves_metrics_vs_size` | Metrics vs population / half size, full space | cue_transfer_incremental, halves_incremental |
| `normalized_metrics_vs_size.csv` | Every point with its CI | |
| **06_pc_sweep_accuracy** | | |
| `pc_cue_accuracy_vs_k`, `pc_halves_accuracy_vs_k`, `pc_across_accuracy_vs_k` | Self, PT, control, non-PT vs k | pc_sweep_* |
| `pc_sweep_accuracy.csv` | Every point (mean, SEM) | |
| **07_pc_sweep_metrics** | | |
| `pc_cue_metrics_vs_k`, `pc_halves_metrics_vs_k`, `pc_across_metrics_vs_k` | Metrics vs k ("full" markers = full space, same neurons) | pc_sweep_*, cue / halves at full size |
| `pc_cue_metrics_vs_k_per_pair` | Cue transfer metrics per rendering pair | |
| `pc_headroom_all_comparisons` | Cue transfer vs halves vs across areas, per monkey | |
| `pc_sweep_metrics.csv` | Every point with its CI | |
| **08_dimensionality** | | |
| `pr_vs_N`, `pr_per_cue_vs_N` | PR vs number of neurons (raw solid, noise-corrected dashed) | cue_transfer_incremental |
| `cumvar_vs_k` | Variance captured by the top k PCs, whole populations | pc_sweep_cue |
| `metrics_vs_pr` | PT/control and headroom vs PR along the N sweep | cue_transfer_incremental |
| `pr_vs_N.csv`, `cumvar_vs_k.csv` | PR per population, size and rendering; cumulative variance per k | |
| **09_nc_injection** | | |
| `nc_injection_metrics` | Self, PT and control, headroom, PT/control per condition (original, NC removed, injected r = 0–0.15) | nc_injection |
| `nc_injection_change` | Paired change from the NC-removed data vs achieved mean NC | nc_injection |
| `nc_injection_check` | Achieved NC (all / same- / cross-session pairs), spread, noise PR, per-pair NC distributions | nc_injection |
| `nc_injection.csv` | Every plotted number, with bootstrap CIs over seeds | |

## Conventions

- **Colors and bar edges** follow the published figures: V1 blue, V2 orange (across
  areas: blue = V1-to-V2, decoder trained on V1; orange = V2-to-V1). Edge style: dashed
  = non-PT, thick solid = PT (rotation only), dotted = corrected control. Dashed line = 2% chance.
- **Paper-style error bars** (02, 03, 06) are the SEM over all repeat estimates, as in the
  old scripts. They understate uncertainty, and points that use a whole population once
  (V1 in Fig 3b, all PC cue / across points) only reflect fold variability.
- **Normalized metrics** (04, 05, 07) are ratios of pooled means (`transfer_metrics.m`);
  error bars and bands are percentile-bootstrap 95% CIs over neuron draws (folds when there
  is a single draw).
- **Headroom** = (PT − control) / (self − control), with the TARGET's self-decoding.
  Pooled over both directions it is symmetric. Per direction it can exceed 1 when the
  tested population is more decodable than the target (KO V2-to-V1 across areas: 0.63
  in full space).
- **Halves in the old-vs-new figure** sit slightly below the identity line (up to 0.056):
  the rerun no longer trains the target half on the trial the source half is tested on.
