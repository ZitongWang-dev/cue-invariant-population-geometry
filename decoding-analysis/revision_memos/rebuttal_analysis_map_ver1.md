# Rebuttal analysis map — reviewer points 1–10 → scripts, outputs, figures

Compiled 2026-08-06 from the `revising` branch. All paths are relative to
`decoding-analysis/`. Shorthand: **ND** = `neuronal_decoding/`, **NDC** =
`neuronal_decoding_control/`, **SIM** = `Simulations/`, **GAB** = `Gabor_decoding/`.

Shared conventions worth remembering while writing:
- 8-column `.acc` order everywhere: `1 self(genAcc) · 2 no-transform · 3 full PT ·
  4 rand/shuffle null · 5 scale-only · 6 rotation-only · 7 translation-only · 8 chance ctrl`.
  **Column 6 (rotation-only) is the headline PT metric** in the paper. Chance = 0.02.
- Six ordered pair files per population: `acec ecac | ecex exec | acex exac`;
  complementary directions are combined into three unordered pairs (EC-AC, EC-EX, AC-EX).
- Populations: FR V1 = 48, FR V2 = 112, KO V1 = 109, KO V2 = 146 neurons.
  Window `[330 630]` ms, 50 stimuli × 10 trials.
- `saveas` is commented out during development in several scripts — that is why
  some analyses have `.mat` outputs but no `.png` yet (flagged per point below).

---

## Point 1 — Terminology: "cue-invariant" vs "rendering-equivariant"

**No computation.** This is a framing/wording change, plus the reviewer's aside that
early visual cortex *cannot* be truly rendering-invariant (otherwise renderings would be
indistinguishable) should be added to the motivation.

Optional supporting numbers you already have, if you want to say "the representation is
equivariant, not invariant, and here is how far from invariant it is":
- **`no-transform` transfer (col 2)** vs **rotation-only PT (col 6)** in
  `ND/results/decoding_outputs/procrustes_decoding_basic_results/<M>/<A>/` — the gap is
  literally the amount of transformation required, i.e. the equivariance.
- `ND/results/figures/{FR,KO} with zscore only rotation.png` already show this contrast.

---

## Point 2 — Test the hypothesis directly with Procrustes-invariant quantities (distances, angles), no decoder

This is the **decoder-free geometry** work stream, in `ND/code/geometry_scripts/`.

### Primary
| | |
|---|---|
| Compute | `ND/code/geometry_scripts/geometric_metrics_basic.m` |
| Outputs | `ND/results/geometric_metrics_outputs/<M>/<A>/geom_{raw,zscore,pt}_results.mat` (all 4 populations; `pt` mode only FR V1 / KO V1) |
| Figures | `ND/results/figures/simple_geometrics/mean_of_trials/geometric_metrics_<M>_<A>/` → `{cos,dist}_distributions_{raw,zscore}.png`, `{cos,dist}_heatmap_{avg,rs}_*.png`, `{cos,dist}_scatter_{avg,rs}_*.png` (all 4 populations) |

**What it does.** For each rendering it builds the 50×50 pairwise **Euclidean distance**
matrix `D` and the 50×50 **centroid-anchored cosine** matrix `Cos` (angle at the population
centroid between stimulus *i* and *j*), from trial-averaged responses. Both are
Procrustes-invariant — exactly the quantities the reviewer names. It then correlates these
matrices **across renderings** (1225 unique upper-triangle pairs, diagonal excluded) and
puts that correlation between two references, both computed inside the same bootstrap so
the noise structure matches:
- **Null** — stimulus-label shuffle nested in the bootstrap (≈ 0).
- **Ceiling** — two independent bootstrap resamples *within* one rendering, i.e. the
  split-half reliability of that rendering's own geometry (~0.88–0.95).
- **Observed** — the real cross-rendering correlation, plus a black-triangle trial-averaged
  point estimate.

**The answer to give the reviewer:** cross-rendering geometry correlation is clearly above
null but well below ceiling — e.g. KO V1 z-scored: AC-EC ≈ 0.38, EC-EX ≈ 0.32, AC-EX ≈ 0.19,
against ceilings 0.88–0.91 and null ≈ 0. Partial, not perfect, preservation of the shape —
which is the honest decoder-free version of the claim. Note `raw` distance is nearly flat
(0.06–0.11) and z-scoring is what reveals the structure; the paper's decoding pipeline
z-scores too, so this is consistent, but you should state it.

### Supporting / documented negatives (use to show the search was thorough)
| Script | Figures | What it shows |
|---|---|---|
| `geometric_metrics_pc1_diagnostic.m` | **none saved** (`saveas` commented, lines 268–269) | Per-stimulus PC1 variance fraction + full scree vs a within-neuron trial-shuffle null. **Observed scree ≈ null scree** → the pseudopopulation's within-cloud covariance is close to independent per-neuron variability. Doubles as a Reviewer-3/7 point. |
| `geometric_metrics_pc1_angle.m` | `simple_geometrics/pc1_angle/<M>_<A>/angle_{within,across}_{raw,zscore}_{full,split}.png` (all 4 populations, 8 figs each) | Angle between each stimulus cloud's PC1 and the signal directions to other stimuli. Clean negative: z-score+split medians ~83–84° vs an 85.6° random-orientation null. The `full→split` (~2°) and `raw→zscore` (~6–7°) progressions are the two artifact controls — keep both panels, they *are* the argument. |
| `geometric_metrics_characteristic_trials.m` | `simple_geometrics/characteristic_trials/<M>_<A>/{cos,dist}_ksweep_{raw,zscore}_{global,own}.png` | Tests whether high-magnitude ("characteristic") trials preserve cue-invariant geometry better. Refuted in every panel — **far is worst**, ordering near ≥ random > far. Cue-invariant geometry lives in typical, not extreme, responses. |

**Gap:** `geometric_metrics_pc1_diagnostic.m` has no saved figure. Uncomment its `saveas`
if you want to cite the scree comparison.

---

## Point 3 — Does successful transfer decoding actually imply special geometric structure? (the rotation has ~N² parameters)

**This is the most important point and the one where you have the strongest new result:
the PT overfitting artifact.** The rotation has ~N²/2 free parameters against 50 stimulus
anchors; when N ≥ 50 it can align two *structureless* clouds using whatever correspondence
it is given. Three independent lines of evidence, all already computed:

### 3a. Synthetic null (Model 1) — the reviewer's own "simplest case"
| | |
|---|---|
| Generator | `SIM/code/generation_scripts/generate_model1_trial_data.m` |
| Driver | `SIM/code/decoding_scripts/procrustes_decoding_simulation.m` |
| Outputs | `SIM/results/decoding_outputs/model1/rho0_N{10,50,100}_ws*_poisson/simulation_results.mat` |

Means drawn i.i.d. per (neuron, stimulus, rendering) — *no shared structure whatsoever*.
Result: **PT ≈ 0.385 at N=100, but ≈ 0.029 (chance) at N=10.** Rotation-only ≈ full PT
(0.398 vs 0.385), confirming the ~N²/2-parameter rotation carries the artifact. The
existing shuffled-correspondence null misses this because shuffling maps to *wrong* labels
→ chance regardless. **The 1/50 chance baseline in the paper is too lenient.**

### 3b. Real-data null floors (structureless surrogates at the real N)
| Script (both in `NDC/code/decoding_scripts/`) | Outputs | Coverage |
|---|---|---|
| `procrustes_decoding_basic_null_floor_mmtmatched.m` | `NDC/results/decoding_outputs/procrustes_decoding_basic_null_floor_mmtmatched_results/KO/V1/` | **KO V1 only** |
| `procrustes_decoding_basic_null_floor_prmatched.m` | `..._null_floor_prmatched_results/KO/V2/` | **KO V2 only** |

- **Moment-matched**: replaces the 50 stimulus means with per-neuron Gaussian draws
  (destroys structure, preserves per-neuron signal scale), re-adds real residuals.
  KO V1 result: **PT ≈ 0.91–0.95, no-transform and rand at chance, rotation-only ≈ PT.**
  *Caveat you must state:* self-decoding comes back ≈ 0.99, above the real level, because
  the surrogate is full-rank — so this floor is measured at the wrong operating point.
- **PR-matched** (the rigorous successor): applies a random rotation *within the centered
  stimulus subspace*, which preserves the neuron covariance — and hence PR and
  self-decoding — to machine precision while destroying stimulus correspondence.
  Procrustes fits a neuron-space rotation and cannot undo a stimulus-space one.
  This is the fair floor. **Compare real PT against this, not against 1/50.**

### 3c. Correspondence-shuffled control on real data (all 4 populations)
| | |
|---|---|
| Compute | `NDC/code/decoding_scripts/procrustes_decoding_basic_trialShuffledWithMean.m` |
| Viz | `NDC/code/visualization_scripts/visualization_procrustes_decoding_basic_trialShuffledWithMean.m` |
| Outputs | `NDC/results/decoding_outputs/procrustes_decoding_basic_trialShuffledWithMean/{FR,KO}/{V1,V2}/` (6-column format, see script header) |

Self-consistent shuffle: the mean matrix is permuted by π, a rotation is fit onto the
target from that permuted matrix, and the **same** π is applied to the held-out test
trials. Scored twice — against target labels (stays HIGH → the rotation can force an
arbitrary correspondence) and against the shuffled identity (≈ 2% → correspondence really
was destroyed). This is the real-data demonstration of the same DOF argument.

### 3d. Constraining the DOF directly — PT vs PC dimension
| | |
|---|---|
| Compute | `ND/code/decoding_scripts/procrustes_decoding_basic_incremental_pc.m` (all 4 populations) |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_basic_incremental_pc.m` |
| Support | `ND/code/visualization_scripts/visualization_variance_vs_pc.m` |
| Outputs | `ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pc_results/<M>/<A>/` |

Fixes the full population and sweeps the number of retained **PC dimensions** k, so
dimension is *controlled* rather than measured. Each cue gets its own leakage-clean mean-PCA
basis (stim1's basis built from training trials only). Directly answers the reviewer's
"depends on the ratio of neurons to classes" — at small k the rotation has few parameters
and cannot overfit; the real/shuffle gap at small k is the honest signal.
Figure 1 = accuracy vs k (raw + normalized by self); Figure 2 = residual Procrustes shape
distance vs k. `visualization_variance_vs_pc.m` gives the cumulative signal-variance
spectrum so you can say how many PCs carry the geometry.

> **Status gap for point 3:** none of 3a–3d has saved figures.
> `Simulations/results/figures/`, `NDC/results/figures/`, and
> `ND/results/figures/incremental_pc/` are all empty. There is also **no visualization
> script for the simulation driver at all**, and `NDC/.../visualization_procrustes_decoding_basic.m`
> points at a stale leaf `procrustes_decoding_basic_null_floor_results` (the real leaves are
> `..._mmtmatched_results` and `..._prmatched_results`) — fix that path before running it.
> The null floors are also only run on **KO V1 (mmt)** and **KO V2 (prmatched)**; you will
> want the PR-matched floor on all four populations before making a quantitative claim.

---

## Point 4 — Compare to CCGP (Bernardi et al. 2020): generalization *without* rotating

Two levels, both already computed for all four populations.

### 4a. Trial-level CCGP (held-out stimuli, ECOC decoder)
| | |
|---|---|
| Compute | `ND/code/decoding_scripts/procrustes_decoding_cross_stimulus_generalization.m` |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_cross_stimulus_generalization.m` |
| Outputs | `ND/results/decoding_outputs/procrustes_decoding_cross_stimulus_generalization_results/HoldStim{2,10}/<M>/<A>/` |
| Figures | **none saved** (`ND/results/figures/cross_stimulus_generalization/` is empty) |

The rotation is fit on `50 − n_hold` stimuli and applied to the `n_hold` stimuli it never
saw; everything is scored on the held-out stimuli with a target-cue-trained ECOC. This is
the CCGP-style question: *is the cross-cue alignment a global property of the manifold, or
does it need to see every stimulus?* Both hold-out levels (2 and 10) are run.

Bracket to read `pt_gen` against (columns of `.acc`):
`1 self_decode · 2 no_transform · 3 pt_gen (MAIN) · 4 pt_ceiling (in-sample) ·
5 pt_floor (self-consistent shuffle) · 6 rand_rot · 7 chance`.
`pt_gen ≈ pt_ceiling` ⇒ the rotation generalizes. **Read the script header on `pt_floor`:**
it is an *overfitting-matched* null, not a chance floor — at large N it rides toward the
ceiling, and the ceiling↔floor gap vs neuron count is itself the overfitting diagnostic
(ties directly back to point 3).

### 4b. Mean-level CCGP (nearest centroid, no SVM, no trial noise)
| | |
|---|---|
| Compute | `ND/code/decoding_scripts/procrustes_meancentroid_cross_stimulus_generalization.m` |
| Outputs | `ND/results/decoding_outputs/procrustes_meancentroid_cross_stimulus_generalization_results/HoldStim2/<M>/<A>/` |
| Viz | **does not exist** |

Strips the decoder entirely: after applying the fit-on-F rotation, does each held-out
source *mean* land nearest its own target mean? Two reference sets — **GLOBAL** (nearest
among all 50 targets, chance 1/50) and **WITHIN** (nearest among the held-out targets only,
chance 1/n_hold) — each with three variants: `gen` (main), `no_transform` (how much
cue-invariance already sits in the raw responses — the same neurons are recorded across
cues, so this is the honest baseline), and `rand` (chance floor).
`.acc` columns: `1 gen_global · 2 notransform_global · 3 rand_global · 4 gen_within ·
5 notransform_within · 6 rand_within`.

**Framing for the reviewer:** this is the cleanest CCGP analogue — it answers the structural
question with no classifier and no trial noise, and `gen vs no_transform` is exactly the
"does the rotation add anything" comparison Bernardi et al. would ask.

> **Status gap:** no figures for either; no visualization script for 4b.

---

## Point 5 — Benchmark the pipeline on artificial data, varying #neurons, #classes, representation statistics

Entirely in `SIM/`. Same decoding math as `procrustes_decoding_basic.m`
(z-score per cue, MATLAB `procrustes` similarity transform, ECOC SVM); no neuron
subsampling — N is set directly and the error bar comes from `R_pop` independent
synthetic populations.

| Component | File |
|---|---|
| Shared driver | `SIM/code/decoding_scripts/procrustes_decoding_simulation.m` |
| Model 1 (null: i.i.d. random means per stimulus × rendering) | `SIM/code/generation_scripts/generate_model1_trial_data.m` |
| Model 3 (linear mixed selectivity: shared rank-2 shape latents, independent per-cue readout weights) | `SIM/code/generation_scripts/generate_model3_trial_data.m` |
| SNR calibration sweep | `SIM/code/generation_scripts/calibrate_weight_scale.m` |

**Runs on disk** (`SIM/results/decoding_outputs/`):
- `model1/rho0_N10_ws0.7_poisson`, `model1/rho0_N50_ws{1,10}_poisson`,
  `model1/rho0_N100_ws{0.7,1,2}_poisson` — **the neuron-count sweep** (N = 10 / 50 / 100)
  and two self-decoding operating points (matched-self ws=0.7 → self ≈ 0.48;
  near-ceiling ws=2 → self ≈ 0.97).
- `model3/rho0_N100_ws12_poisson` — locked operating point self ≈ 0.49, PR ≈ 1.98.
- `model1/calibration/`, `model3/calibration/` — `calibration_sweep.mat` **and
  `calibration_sweep.png`** (the only saved simulation figure).

**Key design points to state in the methods:**
- Poisson trials (count-faithful version of the reviewer's Gaussian mixture); positivity is
  handled by a **baseline offset, never rectification** — rectifying would curve the
  manifold and inflate PR, silently turning Model 3 into a nonlinear model.
- Models are matched on **self-decoding (~0.49), not on `weight_scale`**, so PT/PR
  differences are not SNR confounds. PR is the axis being measured, so it is left free.
- `rho` knob interpolates independent → shared cross-cue readout (rho=1 = trivial transfer).

**Reviewer-facing conclusion:** Model 1 is the reviewer's exact "hopefully this can be ruled
out" case — and it is **not** ruled out at N ≥ 50. That is the finding, and it is a
methods contribution rather than a retraction as long as you (a) report the null floor and
(b) shift the claims to relative comparisons (V1 vs V2, neural vs Gabor at matched PR,
real PT vs PR-matched floor).

> **Gaps:** Models 2/4 are algebraically redundant (skip, but say so). **Model 5** (EX gets
> an extra latent) and **Model 6** (nonlinear mixed selectivity — the model that decouples
> PR from transferability) are **not built**. There is **no stimulus-count (n_stim) sweep** —
> the reviewer explicitly asks for varying the number of stimulus classes, and everything is
> fixed at 50. There is **no simulation visualization script**.

---

## Point 6 — Functional/computational significance of the rotation; how could the brain use it?

Largely a **Discussion** point, but you have four concrete pieces of computational support:

1. **The rotation lives in very few dimensions**, so a downstream readout is not learning an
   N×N matrix. `ND/results/figures/incremental_pr/fig1_accuracy_pr_vs_N.png` and
   `fig5_saturation_four_views.png` — PR saturates around **8–11** for all populations
   (theoretical cap 49). Plus `visualization_variance_vs_pc.m` for the cumulative signal
   spectrum, and `procrustes_decoding_basic_incremental_pc.m` for how much PT survives at
   small k. The readout only has to rotate within a low-dimensional subspace.
2. **The alignment generalizes to stimuli it never saw** — the CCGP scripts under point 4.
   A rotation that must be re-fit per stimulus would be biologically useless; one fit on
   part of the manifold that places the rest correctly is a plausible fixed readout.
3. **The alignment is not tied to specific neurons.**
   `ND/code/decoding_scripts/procrustes_decoding_half_and_half_populations.m` (+ `_incremental`)
   and `procrustes_decoding_neuron_shuffled.m`, with figures
   `ND/results/figures/{FR,KO} half-half decoding with zscore only rotation.png`,
   `{FR,KO} incremental half-half ...png`, `{FR,KO} neuron shuffled ...png`.
4. **Hierarchy / feedforward-feedback framing** for "how would a higher area do this":
   `ND/code/decoding_scripts/procrustes_decoding_across_area.m` →
   `ND/results/figures/{FR,KO} across area decoding.png`; ANN analogue
   `ANN_decoding/.../procrustes_ann_decoding_along_hierarchy.m` →
   `ANN_decoding/results/figures/ALEXNET pool 2&5 feedforward and feedback decoding.png`.

**On the z-scoring sub-question** (the reviewer asks whether the readout would also have to
z-score): `geometric_metrics_basic.m` runs `raw`, `zscore`, and `pt` normalization modes
side by side — the `raw` figures are your direct answer to how much of the effect survives
without z-scoring. Be candid: for the distance metric, raw cross-rendering correlation is
much weaker (~0.06–0.11 vs ~0.19–0.38 z-scored).

---

## Point 7 — Pseudopopulation limits: missing correlation information; possible block-diagonal rotation structure

### 7a. How much noise correlation is actually there
| | |
|---|---|
| Compute | `ND/code/noise_correlation_scripts/noise_correlation_analysis.m` (random-shuffle control) |
| Compute | `ND/code/noise_correlation_scripts/noise_correlation_analysis_affine.m` (deterministic affine-shift control) |
| Outputs | `ND/results/noise_correlation_outputs/<M>/<A>/noise_correlation_results.mat` (all 4) |
| Figures | `ND/results/figures/noise_correlation_<M>_<A>/` and `ND/results/figures/nc_affine_<M>_<A>/` (all 4 populations) |

Per neuron pair × stimulus, Pearson r over the 10 trials (isolates noise from signal
correlation), plus the geometric-mean firing rate. Every pair is tagged **same-session vs
cross-session** by rebuilding the neuron→session map from `bad_channel` in the `ac` folder.
Figure inventory (prefix `nc_<M>_<A>_nan<mode>_`): example scatters (stage1), geo-mean–NC
correlation distributions, `mean_vs_variance_nc` scatterhists (v1 single-colour, v2 split by
session), `significance_real_vs_shuffle`, and `delta_mean_nc_by_session`.

**The result to report:** same-session pairs carry genuine positive noise correlations;
cross-session pairs sit at zero; shuffling collapses same-session down to the cross-session
baseline. So pseudopopulation construction removes most — but not all — correlation
structure, and the residual is characterizable. `nc_<M>_<A>_nanexclude_delta_mean_nc_by_session.png`
shows this most directly. Note the `nan_mode` flag (`exclude` vs `zero`) — use `exclude`
for the delta figures, because shuffling a constant vector stays constant and `zero` mode
would force delta to 0 artificially.

### 7b. Decoding on *true* simultaneously-recorded populations
| | |
|---|---|
| Compute | `ND/code/noise_correlation_scripts/procrustes_decoding_per_session.m` |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_per_session.m` (6 figures) |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_per_session_distribution_test.m` |
| Outputs | `ND/results/decoding_outputs/procrustes_decoding_per_session_results/{none,affine}/<M>/<A>/` (all 4, both modes) |
| Figures | **none saved** |

Slices the pseudopopulation back into real single-session populations (FR: V1 = 5, V2 = 7
sessions; KO same counts; ~5–15 neurons each) via `neuron_session_id`, and runs the same PT
decoding with and without the affine trial perturbation. Six figures: per-session decoding
for each mode, delta, self-vs-PT scatter, normalized delta with Wilcoxon+FDR, and a pooled
figure over both monkeys.

**The framing that was agreed ("Option C"):** the perturbation's effect is **not
unidirectional** — some sessions/pairs show affine < none (correlations were aiding
decoding), others affine > none (correlations were limiting it). This matches the
information-limiting vs signal-correlation dichotomy in the literature. **Do not claim a
single signed effect.** Pooled across stimulus pairs (N = 15 for V1, 21 for V2 per monkey),
nothing survives FDR over 12 tests — the closest is KO V2 PT(rot), raw p = 0.0064, q = 0.076.
Per-monkey per-pair tests are hopelessly underpowered (min two-sided signrank p = 2/2^N:
N=5 → 0.0625). Say this explicitly rather than hiding it.

### 7c. Helper files these depend on
`ND/code/noise_correlation_scripts/apply_trial_perturbation.m` (modes `none|shuffle|affine`)
and `build_affine_shift_config.m` (24 patterns `t → (a·p+b) mod 10`: forward a=1 b=0..9,
reverse a=9 b=0..9, stride-3 a=3 b=0..3, assigned sequentially within session so
simultaneously recorded neurons always get distinct shifts). Duplicated in
`ND/code/geometry_scripts/` so both stacks are self-contained.

> **Gap on the reviewer's specific example:** he suggests you might find **block-diagonal
> structure in the fitted rotation matrices** corresponding to simultaneously recorded sets.
> **Nothing in the repo inspects the rotation matrix `T` itself.** This would be a small,
> high-value new analysis: save `T` from `pro_decoding`, order neurons by session, and
> compare within-session vs cross-session block energy against a session-label shuffle.

---

## Point 8 — Noise correlations shape real geometry; simulate injecting structured NC (mean pairwise r = 0.05–0.15)

Same machinery as point 7, applied to the decoding and geometry pipelines.

| Analysis | Compute | Outputs / figures |
|---|---|---|
| PT decoding with NC broken (pseudopopulation) | `ND/code/noise_correlation_scripts/procrustes_decoding_trial_perturbed.m` | `ND/results/decoding_outputs/procrustes_decoding_trial_perturbed_results/affine/<M>/<A>/` (all 4). Viz: `ND/code/visualization_scripts/visualization_procrustes_decoding_trial_perturbed.m` (figures not saved) |
| Incremental PT + PR with NC broken | `ND/code/noise_correlation_scripts/procrustes_decoding_basic_incremental_pr_trial_perturbed.m` | `ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pr_trial_perturbed_results/affine/<M>/<A>/` (all 4) |
| Geometry with vs without NC | `ND/code/geometry_scripts/geometric_metrics_basic_nc_compare.m` | `ND/results/figures/simple_geometrics/nc_control/<M>_<A>/{cos,dist}_distributions_nc_{raw,zscore}.png` (**all 4 populations, both norms — saved**) |
| Within-cloud covariance structure | `geometric_metrics_pc1_diagnostic.m`, `geometric_metrics_pc1_angle.m` | see point 2 |

**The clean mathematical argument for the geometry side** (`geometric_metrics_basic_nc_compare.m`):
the affine shift only reorders trials *within* each (neuron, stimulus) block, so the block
mean — and therefore the entire trial-averaged geometry, its cross-rendering correlations,
and the black-triangle point estimates — is **invariant to machine precision**. Z-scoring is
also permutation-invariant. The figures confirm this: the trial-averaged triangles coincide
exactly. Where NC *does* show up is the **bootstrap spread**: removing NC makes the
per-neuron mean estimates fluctuate more independently, so the **within-condition ceiling
rises** (KO V1 z-scored distance: AC self 0.82 → 0.87, EC self 0.90 → 0.95) while the
observed cross-rendering correlation barely moves (0.26 → 0.27). Original and NC-removed
runs share the same RNG seed, so the per-category comparison is a paired Wilcoxon
signed-rank on matched bootstrap iterations (stars are on the figures).

**So the reviewer-3/8 headline is:** *cue-invariant geometry is a trial-averaged quantity and
is mathematically invariant to noise correlations; NC affects estimation reliability, not the
geometry. Single-trial decoding is where NC can matter, and there the effect is
session- and pair-dependent, not unidirectional (point 7b).*

> **Gap — the reviewer's literal ask is not built.** He asks for a **simulation injecting
> structured noise correlations of realistic magnitude (r = 0.05–0.15) into the
> pseudopopulation** and testing whether PT is robust or biased. What exists is the
> *opposite* operation (removing the residual NC that is there). The `rho` knob in the
> simulation generators is cross-cue *readout weight* correlation, not trial-to-trial NC.
> Adding a correlated-noise term to `generate_model1/3_trial_data.m` (e.g. a shared
> latent per trial with loading tuned to hit a target mean pairwise r) and re-running the
> driver at r = 0, 0.05, 0.10, 0.15 would answer it directly and is a modest addition.
> Meanwhile `noise_correlation_analysis.m` already gives you the **empirical** r
> distributions, which is what you need to justify (or contest) the 0.05–0.15 range.

---

## Point 9 — Neural vs Gabor at matched population size, with participation ratio

The most complete work stream — **compute, visualization, and saved figures all exist.**

| | |
|---|---|
| Neural compute | `ND/code/decoding_scripts/procrustes_decoding_basic_incremental_pr.m` → `ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pr_results/<M>/<A>/` (all 4) |
| Gabor compute | `GAB/code/decoding_scripts/procrustes_Gabor_resps_decoding_basic_incremental_pr.m` → `GAB/results/decoding_outputs/procrustes_Gabor_incremental_pr_results/even and odd combined/` |
| Viz (Figs 1–5) | `ND/code/visualization_scripts/visualization_procrustes_decoding_basic_incremental_pr.m` |
| **Figures** | `ND/results/figures/incremental_pr/fig1_accuracy_pr_vs_N.png` … `fig5_saturation_four_views.png` |

**PR definition** (state it in methods): `PR = (Σλ)²/Σλ²` over the covariance spectrum via SVD
of the **z-scored, 10-trial-averaged, stimulus-centered 50×N signal manifold**. No second
renormalization after averaging. PR is invariant to the Procrustes transform, so it is a
clean per-condition quantity. Uniform 10-step neuron sweeps so neural and Gabor match
exactly. **Gabor PR uses the mean estimated from generated Poisson trials, not the noiseless
filter rate** — otherwise Gabor gets an unfair clean-mean advantage.

What each figure says:
- **Fig 1 (accuracy + PR vs N):** at matched *neuron count*, Gabor wins on both axes
  (self ≈ 1.0, PR ≈ 11). This comparison is unfavourable and is **not** the right one —
  Gabor is near-noiseless. Show it, then explain why you move to matched PR.
- **Fig 2 (accuracy vs PR — the headline):** **KO V1 sits clearly above the Gabor curve**
  (higher rotation transfer at equal PR); FR V2 ≈ on it; FR V1 slightly below; **KO V2 below.**
  Population-specific, not a blanket neural win — say so.
- **Fig 3 (PT/self vs PR)** and **Fig 4 (PT/PR vs N):** normalized views; ratio-of-pooled-means
  with a percentile bootstrap over neuron-repeats.
- **Fig 5 (2×2 saturation):** PR(N) asymptotes are well constrained (PR_∞ ≈ FR V1 8.2,
  FR V2 ~8.5, KO V1 8.5, KO V2 10.3, Gabor 11.2 under `satexp`). **PT(N) does not converge** —
  the ~0.5 ceiling in panel (b) is a fit projection from extrapolation, not an observed
  plateau, and is fit-form dependent. Present it as a caveat, not a finding.

**The dissociation to lead with:** Gabor dominates *self*-decoding per dimension but loses
*transfer* per dimension to KO V1. Where neural geometry is efficient, the efficiency is
specific to transferable/cue-invariant structure, not to general decodability.

Also available as further baselines: `GAB/results/decoding_outputs/{complex, even and odd
combined, even and odd combined downsampled}/` with
`GAB/results/figures/*only-rotation-PT decoding.png`; and the full ANN suite in
`ANN_decoding/` (AlexNet/VGG, basic + incremental + hierarchy, figures saved).

> **Optional additions discussed but not built:** the `complex` Gabor model in the matched-PR
> figures, a per-pair (AC-EC/EC-EX/AC-EX) breakout, and reporting the `satexp` vs `mm`
> asymptote spread as an uncertainty band (`fit_method` in the viz config).

---

## Point 10 — Compare Procrustes residuals across the three rendering pairs (is AC-EX alignment tighter than genuine cue-invariance predicts?)

The reviewer's concern: AC is the first PC of the same patches that define EX, so AC-EX
alignment may reflect shared stimulus construction rather than cue-invariant coding.

### What exists that bears on it

**Decoder-free, already per-pair — this is your best evidence.**
`geometric_metrics_basic.m` (point 2) computes cross-rendering correlation separately for
`acec`, `acex`, `ecex`, each with its own null and its own per-condition ceilings.
In KO V1 z-scored the ordering is **AC-EC ≈ 0.38 > EC-EX ≈ 0.32 > AC-EX ≈ 0.19** —
i.e. AC-EX is the **weakest**, not the tightest, which is the direct empirical answer
to the reviewer's worry. Check whether that ordering holds in the other three populations
before you commit to the claim in writing.

**Procrustes residuals proper.** The only place a Procrustes residual is stored anywhere in
the repo is `.pdist` in
`ND/code/decoding_scripts/procrustes_decoding_basic_incremental_pc.m`
(→ `ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pc_results/<M>/<A>/`,
all 4 populations, all 6 pair files). Columns:
`1 d_before` (no rotation), `2 d_rot` (rotation only), `3 d_full` (MATLAB `procrustes` d
with optimal b, T, c) — all sharing MATLAB's stim2-target denominator, so they are
comparable across pairs.

**Per-pair PT decoding** for the same comparison at the accuracy level:
`ND/results/decoding_outputs/procrustes_decoding_basic_results/<M>/<A>/` with
`ND/code/visualization_scripts/visualization_procrustes_decoding_basic.m` →
`ND/results/figures/{FR,KO} with zscore only rotation.png` (already broken out by pair),
and the three-panel per-pair layout in the CCGP viz script.

> **Gap:** `visualization_procrustes_decoding_basic_incremental_pc.m` **pools `.pdist` across
> all six ordered pairs** (its header says so explicitly). The data to answer point 10 is
> already on disk — you only need to un-pool it into the three unordered pairs
> (`acec+ecac`, `ecex+exec`, `acex+exac`) and compare `d_full` / `d_rot` across them.
> That is a small edit to the existing `pack_data`/plotting loop, not a new analysis.
> There is **no dedicated residual-across-pairs script or statistic** yet.

---

## Cross-cutting status summary

**Computed, visualized, figures saved (cite directly):**
point 2 primary + pc1_angle + characteristic_trials · point 7a noise-correlation
characterization · point 8 geometry NC-compare · point 9 the whole neural-vs-Gabor PR suite ·
point 6 supporting half-half / neuron-shuffled / across-area / ANN figures ·
simulation calibration sweeps.

**Computed but no saved figures — run the viz before writing:**
point 3 (all four sub-analyses) · point 4 (both CCGP levels; 4b has no viz script at all) ·
point 5 (no simulation viz script at all) · point 7b per-session decoding ·
point 8 trial-perturbed decoding · point 10 per-pair residual breakout ·
`geometric_metrics_pc1_diagnostic.m`.

**Genuine analysis gaps, in rough priority order:**
1. **Point 8** — inject structured NC (r = 0.05–0.15) into the simulation. The reviewer's
   literal request; nothing built.
2. **Point 3** — PR-matched null floor is run on **KO V2 only**. Needs all four populations
   before any quantitative "real PT exceeds the floor by X" claim.
3. **Point 5** — no **stimulus-count** sweep; Models 5 and 6 not built.
4. **Point 7** — no inspection of the fitted rotation `T` for block-diagonal / session
   structure.
5. **Point 10** — no per-pair residual statistic (data exists, just needs un-pooling).

**Small fixes:** stale `file_location` in
`NDC/code/visualization_scripts/visualization_procrustes_decoding_basic.m`
(points at `procrustes_decoding_basic_null_floor_results`, which does not exist);
stale header comment in `visualization_procrustes_decoding_basic_incremental_pr.m`
("Figure 2 will be appended here" — all 5 figures are present).
