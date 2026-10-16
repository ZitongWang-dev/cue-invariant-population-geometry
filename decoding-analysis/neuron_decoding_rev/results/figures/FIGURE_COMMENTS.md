# Figure comments: neuron_decoding_rev (2026-10-12)

What each figure in this folder shows and how we might use it. Click a name to open the PNG.

**Suggested use**
- **Replace (main)**: replaces the published main-figure panel of the same number
- **Replace (appendix)**: replaces a published appendix figure (Figs 12, 13)
- **Add (main)**: a new main-figure panel
- **Supplementary**: Extended Data / appendix
- **Reference**: for our own checks, not for the paper

**Terms.** "Full space" is the published design: V1 whole, V2 subsampled to the V1 count (FR 48, KO 109), unless a figure says otherwise. "Control" is the corrected control. Headroom = (PT − control) / (self − control).

## Why headroom

**What it measures.** Headroom places PT on a scale whose two ends are measured on the same population, with the same decoder and the same rotation freedom:

- **0 is the corrected control:** what a rotation with the same degrees of freedom reaches once the true stimulus correspondence is destroyed. Everything up to this level could be produced by rotation freedom alone.
- **1 is the target's self-decoding:** what the decoder reaches when the rendering does not change at all.

Headroom is the fraction of that gap which the true correspondence closes. In words: how much of the decodable structure survives the change of rendering, beyond what any rotation could produce.

**Why not the other summaries.**

- **Raw PT** grows with how decodable a population is. KO V1 has a higher PT than FR V1 partly because its self-decoding is twice as high (0.50 vs 0.26).
- **PT − control** is in accuracy units, so it scales with decodability in the same way.
- **PT/control** moves with the level of the control, which tracks rotation freedom rather than shared geometry. In the PC sweep it peaks at k ≈ 5 only because the control falls to chance there. In the partially shared simulations it rises from 1.82 to 3.18 as the number of stimuli grows from 20 to 50, at a fixed shared fraction (α = 0.8).

**Evidence that it tracks shared geometry** (Simulations_v3, `partial_sharing_numbers.txt`):

- **No shared structure (α = 0):** headroom is −0.07 to 0.05 at every population size from 20 neurons and every number of stimuli.
- **A shared fraction α of 5 latent dimensions per rendering:** headroom follows α, at 0.09, 0.22, 0.38, 0.62 and 0.94 for α = 0.2–1.0 (50 stimuli). It changes by at most 0.03 between 20 and 50 stimuli, while PT/control rises by 16% (α = 0.2) to 97% (α = 1.0).
- **Recorded populations behave the same way.** On random subsets of 20–50 stimuli, headroom stays flat (FR V1 0.14–0.15, FR V2 0.15–0.18, KO V1 0.25–0.28, KO V2 0.21–0.23), while PT/control rises from 1.09–1.18 to 1.24–1.47.
- **It has interpretable reference points.** Two halves of one area (same geometry, different neurons) reach 0.48–0.70. Cue transfer at 0.14–0.25 corresponds to α ≈ 0.3–0.46 in the simulations.

**Limits.**

- **Small denominator:** headroom becomes unstable when self − control is small. That happens for small populations, and for the near-noiseless CNN and Gabor models; report PT − control there.
- **Per direction:** a single direction is normalized by its own target's self-decoding, so across areas it can exceed 1. Pool both directions instead.
- **How it is computed:** always as a ratio of pooled means (`transfer_metrics.m`), never as a mean of per-fold ratios.

## At a glance

| Paper panel | Use | Figure |
| --- | --- | --- |
| Fig 3b | Replace (main) | `02/fig3b_cue_transfer_full` |
| Fig 4a | Replace (main) | `02/fig4a_halves_full` |
| Fig 4b | Replace (main) | `02/fig4b_across_area_full` |
| Fig 4c | Replace (main) | `03/fig4c_cue_transfer_vs_N` |
| Figs 12, 13 | Replace (appendix) | `02/fig12_across_area_all_combos_full`, `03/fig13_halves_vs_size` |
| New summary panel | Add (main), one of two | `07/pc_headroom_all_comparisons` (headroom row) or `04/cue_vs_halves_matched_size` |
| Extended Data | Supplementary | k = 5 bar figures, metrics vs N, PC sweep, PR, variance spectra |

**Message after the correction.** The control is no longer at chance: it is 0.04–0.32 instead of 0.02, and it sits above non-PT everywhere. The claim therefore rests on PT vs the control, not PT vs non-PT.

- Pooled over pairs, PT beats the control in every population and comparison.
- Cue transfer closes only 14–25% of the control-to-self gap. Two halves of the same area close 48–70%.
- Across-area transfer (24–26%, both directions pooled) is at the cue-transfer level.
- The PC sweep gives the same ordering at every k.

## 01_sanity_old_vs_new

| Figure | Use | Comment |
| --- | --- | --- |
| [old_vs_new_all_bars](01_sanity_old_vs_new/old_vs_new_all_bars.png) | Reference (could back one rebuttal sentence) | Self, non-PT and PT reproduce the published bars: mean \|new − old\| is 0.006 for cue transfer, 0.013 for halves and 0.003 across areas. Halves PT drops by up to 0.056 because the target half is no longer trained on the tested trial. Only the control moves, from 0.02 to 0.04–0.32. |

## 02_paper_style_bars (SEM bars, as published)

| Figure | Use | Comment |
| --- | --- | --- |
| [fig3b_cue_transfer_full](02_paper_style_bars/fig3b_cue_transfer_full.png) | **Replace (main Fig 3b)** | The control (dotted) is 0.06–0.22, above non-PT in all 12 cells, so non-PT is no longer the relevant lower bound. PT/control is 1.11–1.65 per pair. PT exceeds the control in 11 of 12 pairs; the exception is FR V1 EC-EX, whose CI includes 1. The text's "PT 2–4× non-PT" framing must be restated against the control. |
| [fig3b_cue_transfer_k5](02_paper_style_bars/fig3b_cue_transfer_k5.png) | Supplementary (Extended Data to Fig 3b; answers the low-dimensional request) | With 5 PCs per cue the control falls near chance (0.035–0.04) and PT stays above it. In V1 (same neurons as full space), self drops modestly (KO 0.50 → 0.41) while PT drops by more than half (0.25 → 0.11). Non-PT (source projected onto the target's PCs) is close to PT in V1 (KO 0.092 vs 0.106; FR 0.048 vs 0.069), so in 5-D part of the transfer needs no fitted rotation. |
| [fig4a_halves_full](02_paper_style_bars/fig4a_halves_full.png) | **Replace (main Fig 4a)** | The same-geometry reference: PT is 1.9–2.5× the control (headroom 0.48–0.68), and non-PT is at chance because the halves are different neurons. |
| [fig4a_halves_k5](02_paper_style_bars/fig4a_halves_k5.png) | Supplementary (with fig3b_k5) | Halves keep their headroom at k = 5 (0.49–0.67), so the geometry two halves share lives in the top PCs. Cue transfer's headroom drops at k = 5. |
| [fig4b_across_area_full](02_paper_style_bars/fig4b_across_area_full.png) | **Replace (main Fig 4b)** | PT is above the control in both directions and both monkeys (PT/control 1.35–1.52); non-PT is at chance. Pooled headroom is 0.24–0.26, the same level as cue transfer. Per direction it is lopsided (KO V2-to-V1 0.63 vs V1-to-V2 0.13) because each direction is normalized by its training area's self-decoding, so report the pooled value. The published panel used green/magenta for the two directions; match that before submission. |
| [fig4b_across_area_k5](02_paper_style_bars/fig4b_across_area_k5.png) | Supplementary | PT is small (0.07–0.13) but 2.3–2.7× a near-chance control. Pooled headroom is 0.18 for KO and 0.24 for FR. |
| [fig12_across_area_all_combos_full](02_paper_style_bars/fig12_across_area_all_combos_full.png) | **Replace (appendix Fig 12)** | Across-area transfer is strongest when both areas see the same rendering. PT − control is 0.07–0.11 on the diagonal (AC→AC, EX→EX, EC→EC) vs 0.02–0.05 off it, so Fig 4b's averages mix a same-rendering effect with weak cross-rendering transfer. |
| [fig12_across_area_all_combos_k5](02_paper_style_bars/fig12_across_area_all_combos_k5.png) | Supplementary (optional) | Same diagonal pattern, sharper at k = 5: 0.08–0.14 on the diagonal vs 0.02–0.04 off it. |

## 03_paper_style_curves (SEM bars, as published)

| Figure | Use | Comment |
| --- | --- | --- |
| [fig4c_cue_transfer_vs_N](03_paper_style_curves/fig4c_cue_transfer_vs_N.png) | **Replace (main Fig 4c)** | Self, PT and the control all rise with N. The control rises almost as fast as PT (more neurons give the rotation more freedom), so a PT curve alone overstates how transfer grows with N. Non-PT stays at 0.04–0.09. Suggest showing self, PT and control only; in KO the V1 control overlaps the V2 PT curve. |
| [fig13_halves_vs_size](03_paper_style_curves/fig13_halves_vs_size.png) | **Replace (appendix Fig 13)** | For halves, PT stays well above the control at every size; the gap reaches 0.06–0.16 at the largest halves. |

## 04_normalized_metric_bars (bootstrap 95% CI)

| Figure | Use | Comment |
| --- | --- | --- |
| [cue_vs_halves_matched_size](04_normalized_metric_bars/cue_vs_halves_matched_size.png) | **Add (main)**, option A | At the same neuron count (FR 24, KO 54), cue transfer recovers 13–22% of the control-to-self gap and halves 48–68%. This is the simplest full-space statement of the revised result: well above the control, far from the same-geometry reference. |
| [cue_metrics_full_vs_k5](04_normalized_metric_bars/cue_metrics_full_vs_k5.png) | Supplementary | PT/control is 1.25–1.46 in full space vs 1.64–2.66 at k = 5, yet headroom is lower at k = 5 in 3 of 4 populations (0.09–0.18 vs 0.14–0.25; FR V1 is the exception). The ratio rises only because the control falls toward chance, so compare spaces by headroom, not PT/control. |
| [cue_metrics_per_pair_full](04_normalized_metric_bars/cue_metrics_per_pair_full.png) | Supplementary | EC-AC has the highest headroom in 3 of 4 populations (0.24–0.42). EC-EX and AC-EX sit at 0.06–0.13, except KO V1 EC-EX at 0.33. FR V1 EC-EX is the only pair whose CI includes 0. This matches the decoder-free geometry, where EC-AC is the most preserved. |
| [cue_metrics_per_pair_k5](04_normalized_metric_bars/cue_metrics_per_pair_k5.png) | Reference | Same ranking at k = 5: EC-AC is best in all four populations (tied with EC-EX in KO V1). |
| [halves_metrics_full_vs_k5](04_normalized_metric_bars/halves_metrics_full_vs_k5.png) | Reference | Halves headroom is unchanged by the PC reduction: 0.48–0.70 in full space vs 0.49–0.67 at k = 5. |
| [across_metrics_full_vs_k5](04_normalized_metric_bars/across_metrics_full_vs_k5.png) | Supplementary (if Fig 4b reports headroom) | Shows the per-direction asymmetry and that pooled headroom is stable: 0.24–0.26 in full space, 0.18–0.24 at k = 5. |

## 05_normalized_metrics_vs_size (bootstrap 95% CI bands)

| Figure | Use | Comment |
| --- | --- | --- |
| [cue_metrics_vs_N](05_normalized_metrics_vs_size/cue_metrics_vs_N.png) | Supplementary (Extended Data to Fig 4c) | PT/control is flat with N (KO V1 1.45–1.57, KO V2 1.16–1.26, FR 1.2–1.33). Headroom is flat (FR, about 0.10–0.17) or rises slowly (KO V1 0.17 → 0.25, V2 0.10 → 0.23), so the effect is not a population-size artifact. PT − control grows with N. |
| [halves_metrics_vs_size](05_normalized_metrics_vs_size/halves_metrics_vs_size.png) | Reference | Halves headroom rises from 0.3–0.4 at 5 neurons per half to 0.48–0.70. |

## 06_pc_sweep_accuracy (SEM bars)

| Figure | Use | Comment |
| --- | --- | --- |
| [pc_cue_accuracy_vs_k](06_pc_sweep_accuracy/pc_cue_accuracy_vs_k.png) | Supplementary (with pc_cue_metrics_vs_k) | Self saturates by k ≈ 10 at the full-space value; the old 0.68–0.78 was the PCA basis leak. PT and the control rise steadily with k. Non-PT (target basis) is close to PT for small k (up to 4–5) in V1 and flat after. |
| [pc_halves_accuracy_vs_k](06_pc_sweep_accuracy/pc_halves_accuracy_vs_k.png) | Reference | PT stays within about 0.05–0.12 of self at every k; the control rises slowly from chance. |
| [pc_across_accuracy_vs_k](06_pc_sweep_accuracy/pc_across_accuracy_vs_k.png) | Reference | Same shape as cue transfer: PT and the control rise together with k. |

## 07_pc_sweep_metrics (bootstrap 95% CI bands)

| Figure | Use | Comment |
| --- | --- | --- |
| [pc_headroom_all_comparisons](07_pc_sweep_metrics/pc_headroom_all_comparisons.png) | **Add (main)**, option B; otherwise Supplementary | Headroom row: halves are at 0.5–0.7 for every k. For k ≥ 2, cue transfer (0.1–0.25) and across areas (0.15–0.3) stay far below. This covers the PC request and the same-geometry reference in one panel. Use the headroom row only; the PT/control row has the denominator effect. |
| [pc_cue_metrics_vs_k](07_pc_sweep_metrics/pc_cue_metrics_vs_k.png) | Supplementary (rebuttal on the PC request) | PT/control peaks at k ≈ 4–6 only because the control falls to chance. Headroom rises with k and reaches the full-space value at k = 49, so no low-dimensional subspace carries relatively more cue-invariant structure than the full space. |
| [pc_halves_metrics_vs_k](07_pc_sweep_metrics/pc_halves_metrics_vs_k.png) | Reference | Halves headroom is about 0.5–0.7 at every k and matches full space. |
| [pc_across_metrics_vs_k](07_pc_sweep_metrics/pc_across_metrics_vs_k.png) | Reference | Pooled headroom rises slowly with k (about 0.1 → 0.27). Per-direction curves cross because of the normalization. |
| [pc_cue_metrics_vs_k_per_pair](07_pc_sweep_metrics/pc_cue_metrics_vs_k_per_pair.png) | Reference | EC-AC has the highest headroom at most k in KO V2, FR V1 and FR V2; in KO V1, EC-EX overtakes it for k ≥ 20. |

## 08_dimensionality

| Figure | Use | Comment |
| --- | --- | --- |
| [pr_vs_N](08_dimensionality/pr_vs_N.png) | Supplementary (answers the PR request) | Noise-corrected PR saturates at 4.3–5.2 by about 40–60 neurons. Raw PR (7.0–9.5 at full size) keeps rising because trial noise inflates it, so report the noise-corrected value. |
| [cumvar_vs_k](08_dimensionality/cumvar_vs_k.png) | Supplementary (justifies k = 5) | The top 5 PCs hold 64–90% of the noise-corrected signal variance (50–72% raw). |
| [pr_per_cue_vs_N](08_dimensionality/pr_per_cue_vs_N.png) | Reference | EC has the lowest noise-corrected PR in 3 of 4 populations (2.5–3.5); EX is the highest in V2. |
| [metrics_vs_pr](08_dimensionality/metrics_vs_pr.png) | Reference | No consistent relation between PR and PT/control. Headroom follows N rather than PR: it keeps rising after PR has saturated (KO). This is not evidence for a PR–transfer link. |

## Caveats for whichever figures we keep

- **Error bars.** The paper-style SEM bars (02, 03, 06) treat draws × pairs × folds as independent and understate uncertainty. For the paper, consider bootstrap CIs over neuron draws (as in 04, 05, 07) or the SEM over draws only.
- **PT/control across spaces.** PT/control is not comparable across spaces or k, because it rises whenever the control falls toward chance. Use headroom for those comparisons.
- **V2 size.** Full-space cue and across-area bars subsample V2 to the V1 count, while the k = 5 bars use whole populations. Compare V2 between the two spaces with that in mind.
