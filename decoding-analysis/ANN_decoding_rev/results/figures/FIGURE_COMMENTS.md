# Figure comments: ANN_decoding_rev (2026-10-13)

What each CNN figure in this folder shows and how we might use it. Click a name to open the PNG.

**Suggested use**
- **Replace (main)**: replaces the published main-figure panel of the same number
- **Add (main)**: a new main-figure panel
- **Supplementary**: Extended Data / appendix
- **Reference**: for our own checks or the rebuttal, not for the paper

**Terms.** "Full space" is the published design: 96 units per AlexNet layer and 64 per VGG
layer for Fig 3d, or the given unit count. "Whole layer" is every unit of the layer (64–511).
"Control" is the corrected control. Headroom = (PT − control) / (self − control).

## Why headroom

**What it measures.** Headroom places PT on a scale whose two ends are measured on the same population, with the same decoder and the same rotation freedom:

- **0 is the corrected control:** what a rotation with the same degrees of freedom reaches once the true stimulus correspondence is destroyed. Everything up to this level could be produced by rotation freedom alone.
- **1 is the target's self-decoding:** what the decoder reaches when the rendering does not change at all.

Headroom is the fraction of that gap which the true correspondence closes. In words: how much of the decodable structure survives the change of rendering, beyond what any rotation could produce.

**Why not the other summaries.**

- **Raw PT** grows with how decodable a population is.
- **PT − control** is in accuracy units, so it scales with decodability in the same way.
- **PT/control** moves with the level of the control, which tracks rotation freedom rather than shared geometry. It peaks at k ≈ 5 only because the control falls to chance there. In the partially shared simulations it rises from 1.82 to 3.18 as the number of stimuli grows from 20 to 50, at a fixed shared fraction (α = 0.8).

**Evidence that it tracks shared geometry** (Simulations_v3, `partial_sharing_numbers.txt`):

- **No shared structure (α = 0):** headroom is −0.07 to 0.05 at every population size from 20 neurons and every number of stimuli.
- **A shared fraction α of 5 latent dimensions per rendering:** headroom follows α, at 0.09, 0.22, 0.38, 0.62 and 0.94 for α = 0.2–1.0 (50 stimuli). It changes by at most 0.03 between 20 and 50 stimuli.
- **Recorded populations:** headroom stays flat over random subsets of 20–50 stimuli, while PT/control rises with the number of stimuli.
- **Reference points:** two halves of one area (same geometry, different neurons) reach 0.48–0.70, and neural cue transfer reaches 0.14–0.25.

**For the CNN layers specifically, headroom only works where its denominator is large.**

- **Where it fails:** self-decoding is about 1, and in full space the control climbs to 0.7–0.95, so self − control falls below 0.1 for the large layers. Headroom then swings with small changes in PT or control (e.g. VGG pool 4: 0.43 on a denominator of 0.08). For full-space CNN statements, use PT − control.
- **Where it is reliable:** at k = 5 and along most of the PC sweep, the control stays near chance and self − control is about 0.95.
- **The noise-regime caveat:** CNN self-decoding (about 1) is far above the neurons' (0.24–0.41). Comparing CNN with neural headroom therefore still carries this difference, which is why the noise-matched variant is open (end of this document).

## At a glance

| Paper panel | Use | Figure |
| --- | --- | --- |
| Fig 3d | Replace (main) | `02/fig3d_cnn_cue_transfer_full` |
| Fig 4d | Replace (main) | `03/fig4d_cnn_vs_units` |
| New model-comparison panel | Add (main) | `09/k5_cnn_vs_neurons` (self and headroom panels) |
| Extended Data | Supplementary | Fig 3d at k = 5, metrics vs units, PC sweep, PR, CNN vs neurons along k |
| Rebuttal point 1 (hierarchy) | Reference | `02/hierarchy_*`, `04/hierarchy_*`, `07/pc_hierarchy_metrics_vs_k` |

**Message after the correction.** For CNN layers the corrected control is high: 0.17–0.86 at
the Fig 3d sizes and up to 0.95 for whole layers. PT stays only slightly above it.

- **Fig 3d sizes:** PT/control is 1.03–1.34 and PT − control is 0.01–0.08.
- **Published claims that do not hold:** "PT closely approached self-decoding (23–90% of
  self)", "plateau at 250–300 units" and "almost perfectly preserved". The fitted rotation's
  freedom produces most of that accuracy, and the published curves flatten together with the
  control.
- **At the same rotation freedom (k = 5):** neural headroom is 0.09–0.18 and CNN headroom is
  0.01–0.08.
- **Caveat on that gap:** CNN self-decoding is about 1 and neural self-decoding is 0.24–0.41.
  In absolute PT − control the two overlap (neurons 0.02–0.07, CNN 0.01–0.08), so the
  difference comes from the denominator. A noise-matched CNN variant is needed before the
  model figure is final; see "Open decision" at the end.

## 01_sanity_old_vs_new

| Figure | Use | Comment |
| --- | --- | --- |
| [old_vs_new_all_points](01_sanity_old_vs_new/old_vs_new_all_points.png) | Reference (could back one rebuttal sentence) | Self, non-PT and PT reproduce the published numbers: mean \|new − old\| is 0.004 for Fig 3d, 0.006 for Fig 4d at the shared unit counts, and 0.003 for the hierarchy. Only the control moves, from 0.02 to 0.04–0.99 (0.68–0.98 in the hierarchy). |

## 02_paper_style_bars (SEM bars, as published)

| Figure | Use | Comment |
| --- | --- | --- |
| [fig3d_cnn_cue_transfer_full](02_paper_style_bars/fig3d_cnn_cue_transfer_full.png) | **Replace (main Fig 3d)** | Self is 0.87–1.00. The control (dotted) sits just below PT in every layer and pair: AlexNet PT 0.64–0.90 vs control 0.61–0.86; VGG pools 2–5 PT 0.47–0.65 vs 0.45–0.60; VGG pool 1 PT 0.22 vs 0.17. The large PT bars of the published panel are mostly rotation freedom. |
| [fig3d_cnn_cue_transfer_k5](02_paper_style_bars/fig3d_cnn_cue_transfer_k5.png) | Supplementary (Extended Data to Fig 3d) | At 5 PCs per rendering, self stays at 0.95–1.00 but PT falls to 0.05–0.12, against a near-chance control (0.04). The CNN geometry shared across renderings does not lie in the top PCs. The clearest case is AlexNet pool 1 EC-AC (PT 0.21). |
| [hierarchy_by_training_full](02_paper_style_bars/hierarchy_by_training_full.png) | Reference (rebuttal point 1) | AlexNet pool 2 ↔ pool 5: PT 0.86–0.90 vs control 0.83–0.88 (PT/control 1.02–1.04). In full space the layer-to-layer alignment cited in the rebuttal is almost all rotation freedom. |
| [hierarchy_by_training_k5](02_paper_style_bars/hierarchy_by_training_k5.png) | Reference | At k = 5, PT is 0.06–0.21 vs control 0.04. The largest value is pool 5 → pool 2 with EC as the training rendering. |
| [hierarchy_all_combos_full](02_paper_style_bars/hierarchy_all_combos_full.png) | Reference | Same-rendering combinations (AC→AC, EX→EX, EC→EC) carry the transfer: PT − control is 0.04–0.06 on the diagonal vs 0.01–0.02 off it, as for V1 ↔ V2. In some combinations the control exceeds PT (e.g. EC→EX). |
| [hierarchy_all_combos_k5](02_paper_style_bars/hierarchy_all_combos_k5.png) | Reference | Same pattern, sharper: 0.07–0.16 on the diagonal vs 0.01–0.03 off it. |

## 03_paper_style_curves (SEM bars, as published)

| Figure | Use | Comment |
| --- | --- | --- |
| [fig4d_cnn_vs_units](03_paper_style_curves/fig4d_cnn_vs_units.png) | **Replace (main Fig 4d)** | The published curves with the control added. The control (dotted) climbs alongside PT (dashed) for every layer and reaches 0.65–0.95 at full size (VGG pool 1: 0.17), so the "plateau" in the published figure is shared by the control. Self saturates near 1 by 20–160 units. |
| [fig4d_cnn_vs_units_per_layer](03_paper_style_curves/fig4d_cnn_vs_units_per_layer.png) | Supplementary | One panel per layer, with non-PT (0.03–0.16). The gap between PT and control is visible only in AlexNet pool 1–2 and VGG pool 1–3, and only at moderate sizes. |

## 04_normalized_metric_bars (bootstrap 95% CI)

| Figure | Use | Comment |
| --- | --- | --- |
| [cue_metrics_full_vs_k5](04_normalized_metric_bars/cue_metrics_full_vs_k5.png) | Supplementary | Whole layers. In full space PT/control is 1.01–1.34 and PT − control 0.01–0.08. Headroom (0.06–0.43) is unreliable here because its denominator is 0.05–0.36 in 7 of 8 layers; VGG pool 4's 0.43 rests on a denominator of 0.08. At k = 5, PT/control rises to 1.2–2.9 only because the control falls to chance; headroom drops to 0.01–0.08. |
| [cue_metrics_fig3d_size_per_pair](04_normalized_metric_bars/cue_metrics_fig3d_size_per_pair.png) | Supplementary | EC-AC has the largest PT − control in 6 of 8 layers (up to 0.15 in AlexNet pool 1); AlexNet pool 5 and VGG pool 4 are near-ties. AC-EX is the weakest pair: AlexNet pool 1 AC-EX is indistinguishable from the control. This is the same pair ranking as the neurons. |
| [hierarchy_metrics_full_vs_k5](04_normalized_metric_bars/hierarchy_metrics_full_vs_k5.png) | Reference (rebuttal point 1) | In full space PT/control is 1.02–1.04 and PT − control 0.02–0.04 (headroom 0.20 on a denominator of 0.10–0.17). At k = 5, pool 5 → pool 2 shows PT/control 2.8 and PT − control 0.07, larger than pool 2 → pool 5 (0.03). |

## 05_normalized_metrics_vs_size (bootstrap 95% CI bands)

| Figure | Use | Comment |
| --- | --- | --- |
| [cue_metrics_vs_units](05_normalized_metrics_vs_size/cue_metrics_vs_units.png) | Supplementary (Extended Data to Fig 4d; supports the text change) | The test of the plateau claim. PT/control falls toward 1 as units are added (1.01–1.08 for the layers of 249–511 units). PT − control peaks at 40 units in AlexNet pool 1–2 and VGG pool 3 (0.08–0.12) and then declines (VGG pool 2 rises to 0.07 at its full 128 units); the deep layers (AlexNet pool 5, VGG pool 4–5) stay at or below 0.06. Headroom keeps rising in several layers only because self − control shrinks to 0.05–0.08 (bottom row). |

## 06_pc_sweep_accuracy (SEM bars)

| Figure | Use | Comment |
| --- | --- | --- |
| [pc_cue_accuracy_vs_k](06_pc_sweep_accuracy/pc_cue_accuracy_vs_k.png) | Supplementary | Self reaches 0.95–1.00 by k = 5 in every layer. PT and the control stay at or below 0.12 up to k = 5 and then rise together, reaching the full-space values at k = 49. Non-PT on the target's PCs stays at 0.02–0.15. |
| [pc_hierarchy_accuracy_vs_k](06_pc_sweep_accuracy/pc_hierarchy_accuracy_vs_k.png) | Reference | Pool 5 → pool 2 shows a clear PT-over-control gap at k = 4–30 (e.g. 0.42 vs 0.25 at k = 20). Pool 2 → pool 5 shows a small one. |

## 07_pc_sweep_metrics (bootstrap 95% CI bands)

| Figure | Use | Comment |
| --- | --- | --- |
| [pc_cue_metrics_vs_k](07_pc_sweep_metrics/pc_cue_metrics_vs_k.png) | Supplementary | PT − control peaks at intermediate k, not at k = 5 or in full space. Peaks: AlexNet pool 1 0.165 and pool 2 0.128 (k = 15); VGG pool 3 0.115 (k = 30); VGG pool 2 and 4 0.09–0.10 (k = 15); AlexNet pool 5, VGG pool 1 and 5 at or below 0.065. The cue-shared CNN geometry sits in a 10–30-dimensional subspace, which full space dilutes and k = 5 misses. |
| [pc_hierarchy_metrics_vs_k](07_pc_sweep_metrics/pc_hierarchy_metrics_vs_k.png) | Reference (best support for rebuttal point 1) | Pool 5 → pool 2 PT − control peaks at 0.19 (k = 15), and the pooled value is 0.13 at k ≈ 20. This is real layer-to-layer alignment that full space (0.02–0.03) hides. If the rebuttal keeps the hierarchy argument, this is the figure to cite. |

## 08_dimensionality

| Figure | Use | Comment |
| --- | --- | --- |
| [pr_vs_units](08_dimensionality/pr_vs_units.png) | Supplementary (with the neural PR figure) | Raw and noise-corrected PR coincide, because the Poisson CNN trials are nearly noiseless. Whole-layer PR is 4.6 (VGG pool 1), 9.2–20.6 (AlexNet) and 11.9–29.9 (VGG pool 2–5), against 4.3–5.2 for the neurons. Except in VGG pool 1, PR keeps rising with the number of units. |
| [cumvar_vs_k](08_dimensionality/cumvar_vs_k.png) | Reference | The top 5 PCs hold 40–64% of the signal variance in AlexNet and 29–84% in VGG, against 64–90% (noise-corrected) for the neurons. k = 5 is a much stronger reduction for the CNN layers. |
| [metrics_vs_pr](08_dimensionality/metrics_vs_pr.png) | Reference | No consistent relation between PR and transfer. Within a layer, PT − control first rises with PR (more units) and then falls in the early layers. |

## 09_cnn_vs_neurons

| Figure | Use | Comment |
| --- | --- | --- |
| [k5_cnn_vs_neurons](09_cnn_vs_neurons/k5_cnn_vs_neurons.png) | **Add (main)**, as the model comparison that replaces "CNNs show similar invariance" | At the same rotation freedom, neural headroom is 0.09–0.18 and CNN headroom 0.01–0.08. But self-decoding is 0.24–0.41 for neurons and 0.95–1.00 for CNN layers, and PT − control overlaps (0.02–0.07 vs 0.01–0.08). Show the self panel next to the headroom panel, and settle the noise-matched variant first. |
| [pc_sweep_cnn_vs_neurons](09_cnn_vs_neurons/pc_sweep_cnn_vs_neurons.png) | Supplementary | Along k, neural PT − control plateaus at 0.03–0.09, while the early CNN layers peak higher (0.10–0.17 at k = 10–30) and fall back by k = 49. Headroom curves overlap from k ≈ 15 on. The CNN advantage in absolute PT − control appears only where CNN self-decoding is already at ceiling. |
| [vs_size_cnn_vs_neurons](09_cnn_vs_neurons/vs_size_cnn_vs_neurons.png) | Supplementary | Full space vs population size (log axis). Neural self is 0.1–0.5 and CNN self about 1. Neural PT − control grows steadily with N (to 0.08 in KO V1). CNN PT − control peaks early and falls in AlexNet pool 1–2 and VGG pool 3. Neural headroom sits in the CNN range at the same N. |
| [vs_pr_cnn_vs_neurons](09_cnn_vs_neurons/vs_pr_cnn_vs_neurons.png) | Supplementary (rebuttal point 9, matched-PR style) | Neurons reach headroom 0.10–0.25 at PR 3–5. CNN layers reach the same headroom only at PR 7–30, except VGG pool 1 (PR 3–5, headroom 0.01–0.07). Per effective dimension, the neural code carries more cue-transferable structure, in line with the Gabor argument of rebuttal point 9. |

## Caveats for whichever figures we keep

- **Noise regime.** CNN trials are Poisson draws around activations whose scale differs by
  more than 100× across layers, so self-decoding is about 1 and there are no noise
  correlations. Every comparison with the neurons should show self-decoding next to the
  normalized metric.
- **Headroom near ceiling.** For CNN layers in full space, self − control is often below 0.1.
  Use PT − control (and PT/control) for full-space CNN statements; headroom is reliable only
  where its denominator is large, e.g. at k = 5.
- **Error bars.** Paper-style SEM bars (02, 03, 06) treat draws × pairs × folds as
  independent. For the paper, use bootstrap CIs over unit draws or folds (as in 04, 05, 07, 09).
- **Layer names.** MATLAB's AlexNet has three pooling layers named pool1, pool2 and pool5;
  the published Fig 4d legend calls the last one "Pool 3". Use "pool5" or "pool 3 (of 3)"
  consistently.

## Open decision

**Noise-matched CNN variant.** We agreed to decide this after the Poisson run, and the
model comparison now depends on it. Rescaling each layer's activations so that its
self-decoding at a reference size matches the neurons' (about 0.3–0.5) would put CNN and
neural headroom on the same footing. This would be a rerun of the k = 5 / PC sweep, which
takes about 1.5 h, plus the unit sweep if Fig 4d should also use it.
