# Figure comments: Gabor_decoding_rev (2026-10-21)

What each Gabor-model figure in this folder shows and how we might use it. Click a name to open the PNG.

**Suggested use**
- **Replace (appendix)**: replaces a published appendix figure (Fig 15)
- **Add (main)**: a new main-figure panel
- **Supplementary**: Extended Data / appendix
- **Reference**: for our own checks or the rebuttal, not for the paper

**Terms.** "Simple" is the 192-cell simple-cell model, "complex" the 48-cell complex-cell model.
"48 cells" is the published Fig 15a design (simple down-sampled to the complex-cell count),
"whole model" the Fig 15b design. "Control" is the corrected control.
Headroom = (PT − control) / (self − control).

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

**For the Gabor models, headroom is well defined.** Self-decoding is 1.000 and the control stays at or below 0.34, so the denominator self − control is 0.66 or more, unlike for the large CNN layers. The noise regime still differs from the neurons'. The Gabor trials are nearly noiseless (self 1.000 vs 0.2–0.5 for neurons), so comparisons with the neurons should show self-decoding alongside headroom, and the noise-matched variant remains the decisive check (end of this document).

## At a glance

| Paper item | Use | Figure |
| --- | --- | --- |
| Fig 15 | Replace (appendix) | `02/fig15_gabor_cue_transfer_full` |
| Model comparison (Gabor + CNN vs neurons) | Add (main) | `09/k5_models_vs_neurons` (self and headroom panels) |
| Rebuttal point 9, Figs 1–3 | Reference (rebuttal) | `09/point9_vs_cells`, `09/point9_vs_pr` |
| Extended Data | Supplementary | Fig 15 at k = 5, metrics by design, the three-family PC sweep |

**Message after the correction.** The Gabor models carry little cue-transferable structure beyond rotation freedom.

- **Against the control:** PT exceeds the corrected control by only 0.04–0.06 in full space (PT/control 1.15–1.35; headroom 0.04–0.09), and by 0.01 at k = 5.
- **Against the neurons at the same size:** at 48–50 cells or neurons, the four neural populations have headroom 0.12–0.21 against 0.04–0.05 for the Gabor models. At about 110, it is 0.18–0.25 against 0.07.
- **Against the neurons at the same dimensionality:** the Gabor complex model (PR 3–5) reaches at most 0.04, where the neurons at PR 3–5 reach 0.09–0.25.
- **Text that has to change:** the paper's statements that Gabor PT is "significantly above" the floor, that "simple cells approximately doubled PT cue-transfer decoding" with population size, and that "classical Gabor models are sufficient to produce a cue-invariant geometric structure". The doubling of PT from 48 to 192 cells (0.28 to 0.40) is matched by the control (0.24 to 0.34).
- **The rebuttal's point-9 dissociation now holds against the correct control on both sides, and is larger.** The old figures compared a corrected neural control with the old Gabor control.
- **Caveat:** in absolute PT − control, only KO V1 clearly exceeds the Gabor bank at matched size (0.053 vs 0.036 at about 50; 0.080 vs 0.047 at about 110). The other populations are at or below it, so the headroom advantage depends on the neurons reaching it from self-decoding of 0.2–0.5. That is why the noise-matched variant matters.

## 01_sanity_old_vs_new

| Figure | Use | Comment |
| --- | --- | --- |
| [old_vs_new_all_points](01_sanity_old_vs_new/old_vs_new_all_points.png) | Reference (could back one rebuttal sentence) | Self, non-PT and PT reproduce the published numbers: mean \|new − old\| is 0.002 for Fig 15 and 0.004 for the point-9 sweep. Only the control moves, from 0.02 to 0.06–0.47. |

## 02_paper_style_bars (SEM bars, as published)

| Figure | Use | Comment |
| --- | --- | --- |
| [fig15_gabor_cue_transfer_full](02_paper_style_bars/fig15_gabor_cue_transfer_full.png) | **Replace (appendix Fig 15)** | Self is 1.000 for both models. The control (dotted) sits just below PT in every pair. (a) At 48 cells: simple PT 0.28 vs control 0.24, complex 0.14 vs 0.10; simple EC-EX is 0.33 vs 0.32. (b) Whole models: simple 0.40 vs 0.34. The published "PT significantly above chance" is mostly rotation freedom. |
| [fig15_gabor_cue_transfer_k5](02_paper_style_bars/fig15_gabor_cue_transfer_k5.png) | Supplementary (Extended Data to Fig 15) | At 5 PCs per rendering, self stays at 1.000 while PT falls to 0.05 against a control of 0.04. There is essentially no cue-shared structure in the models' top PCs; non-PT (target PCs) equals PT for the simple model. |

## 03_paper_style_curves (SEM bars, as published)

| Figure | Use | Comment |
| --- | --- | --- |
| [gabor_vs_cells](03_paper_style_curves/gabor_vs_cells.png) | Reference (superseded by `09/point9_vs_cells`) | PT and the control rise together with the number of cells. For the simple model the gap widens slowly, from 0.01 at 10 cells to 0.06 at 192. Self is 1.000 from 20 cells on, and non-PT stays at 0.04–0.06. |

## 04_normalized_metric_bars (bootstrap 95% CI)

| Figure | Use | Comment |
| --- | --- | --- |
| [gabor_metrics_by_design](04_normalized_metric_bars/gabor_metrics_by_design.png) | Supplementary (Extended Data to Fig 15) | At 48 cells, simple and complex are equal in PT − control (0.036), with headroom 0.048 and 0.040. The whole simple model reaches headroom 0.087 (PT/control 1.17). At k = 5 both models drop to headroom 0.01. The denominator stays at 0.66–0.96, so these values are stable. |
| [gabor_metrics_per_pair_48](04_normalized_metric_bars/gabor_metrics_per_pair_48.png) | Reference | Simple model: EC-AC transfers best (PT − control 0.057, headroom 0.074) and EC-EX barely (0.011, CI reaching 0). Complex: the three pairs are similar (0.029–0.043). The neural ranking (EC-AC best) appears only in the simple model. |

## 05_normalized_metrics_vs_size (bootstrap 95% CI bands)

| Figure | Use | Comment |
| --- | --- | --- |
| [gabor_metrics_vs_cells](05_normalized_metrics_vs_size/gabor_metrics_vs_cells.png) | Reference | Headroom rises slowly with the number of cells: simple 0.01 → 0.09, complex 0.01 → 0.04. PT/control is flat at 1.1–1.2 (simple) and 1.3 (complex). The denominator shrinks only modestly (simple 0.92 → 0.66). |

## 06_pc_sweep_accuracy (SEM bars)

| Figure | Use | Comment |
| --- | --- | --- |
| [pc_gabor_accuracy_vs_k](06_pc_sweep_accuracy/pc_gabor_accuracy_vs_k.png) | Reference | Self reaches 1.000 by k = 2–3. PT and the control stay together near chance up to k ≈ 7 (simple) and rise in parallel after. PT clears the control only at large k for the simple model (0.40 vs 0.34 at k = 49, the full-space value). For the complex model it does so from k ≈ 6, staying 0.03–0.05 above the control. |

## 07_pc_sweep_metrics (bootstrap 95% CI bands)

| Figure | Use | Comment |
| --- | --- | --- |
| [pc_gabor_metrics_vs_k](07_pc_sweep_metrics/pc_gabor_metrics_vs_k.png) | Reference | Headroom stays at or below 0.05 up to k = 30 for both models. The simple model climbs to 0.10 only at k = 40–49; the complex model plateaus at 0.045 from k ≈ 15. Unlike the CNN layers (PT − control peaks of 0.10–0.17 at k = 10–30), the Gabor models have no intermediate-k subspace with a clear transfer advantage. |

## 08_dimensionality

| Figure | Use | Comment |
| --- | --- | --- |
| [pr_vs_cells](08_dimensionality/pr_vs_cells.png) | Reference (rebuttal point 9) | Raw and noise-corrected PR coincide (nearly noiseless trials). The simple model's PR rises 5.2 → 10.6 and the complex model's 3.2 → 4.8, against 4.3–5.2 for the neurons. This agrees with the rebuttal's "Gabor PR∞ ≈ 11" for the simple model. |
| [cumvar_vs_k](08_dimensionality/cumvar_vs_k.png) | Reference | The top 5 PCs hold 62% (simple) and 83% (complex) of the signal variance. The complex model is as low-dimensional as the neurons, yet transfers far less. |
| [metrics_vs_pr](08_dimensionality/metrics_vs_pr.png) | Reference | For both models, transfer rises with PR as cells are added. Its level stays far below the neurons' at every PR (see `09/point9_vs_pr`). |

## 09_models_vs_neurons

| Figure | Use | Comment |
| --- | --- | --- |
| [k5_models_vs_neurons](09_models_vs_neurons/k5_models_vs_neurons.png) | **Add (main)**, as the single model-comparison panel (supersedes `ANN_decoding_rev/09/k5_cnn_vs_neurons`) | At the same rotation freedom, neural headroom is 0.09–0.18, CNN layers 0.01–0.08 and Gabor models 0.01. In PT − control the neurons (0.02–0.07) are above the Gabor models (0.01) and overlap the CNN layers (0.01–0.08). Self-decoding: neurons 0.24–0.41, all models 0.95–1.00. |
| [point9_vs_cells](09_models_vs_neurons/point9_vs_cells.png) | Reference (replaces rebuttal point 9, Fig 1) | At matched size the Gabor bank still wins on self (1.000) and on raw PT (0.36 vs KO V1 0.25 at about 110), as the rebuttal said. But its control is high too (0.31), so in headroom every neural population sits above both models at every matched size (neurons 0.09–0.25 over their whole range, Gabor 0.01–0.09). |
| [point9_vs_pr](09_models_vs_neurons/point9_vs_pr.png) | Reference (replaces rebuttal point 9, Figs 2–3) | At matched noise-corrected PR the neurons sit above both Gabor models in headroom. The complex model shares the neural PR range (3–5) but reaches only 0.01–0.04, against 0.09–0.25. The simple model reaches 0.04–0.09 only at PR 9–10.6. This is the reviewer's "comparable invariance in fewer effective dimensions", now against the correct control. |
| [pc_sweep_models_vs_neurons](09_models_vs_neurons/pc_sweep_models_vs_neurons.png) | Supplementary | Along k, neural headroom (0.08–0.25 from k ≈ 3) stays above both Gabor models (at or below 0.05 up to k = 30) at every k. The early CNN layers catch up with the neurons only from k ≈ 15, where their self-decoding has been at ceiling since k = 3. |

## Caveats for whichever figures we keep

- **Noise regime.** Gabor trials are Poisson draws with means of 300–4,000 per trial, so
  self-decoding is 1.000 and there are no noise correlations. Show self-decoding next to
  every neuron-vs-model comparison.
- **Absolute vs normalized transfer.** In PT − control, KO V1 exceeds the Gabor bank but
  KO V2, FR V1 and FR V2 do not. The neural advantage is in headroom, i.e. relative to
  how decodable each population is.
- **Error bars.** Paper-style SEM bars (02, 03, 06) treat draws × pairs × folds as
  independent. For the paper, use bootstrap CIs over cell draws or folds (as in 04, 05, 07, 09).

## Open decision

**Noise-matched variant (CNN and Gabor together).** Regenerating the model trials at a lower
rate scale, so that self-decoding at a reference size matches the neurons' (about 0.3–0.5),
would test whether the neural headroom advantage survives when the models are as noisy as
the neurons. This is the remaining step before the model figure is final. A Gabor rerun of
the k = 5 / PC sweep takes about 15 min; the CNN one about 1.5 h.
