# Rebuttal analysis map — reviewer points 1–10 → scripts, outputs, figures

Regenerated 2026-08-08 from the `revising` branch by walking the repo (script inventory,
figure inventory, `decoding_outputs` inventory, and the `%{ … %}` header block of every
`.m`). All paths are relative to `decoding-analysis/`.

Shorthand: **ND** = `neuronal_decoding/`, **NDC** = `neuronal_decoding_control/`,
**SIM** = `Simulations/`, **GAB** = `Gabor_decoding/`, **ANN** = `ANN_decoding/`.
`<M>` ∈ {FR, KO}, `<A>` ∈ {V1, V2}.

Shared conventions:
- 8-column `.acc` order everywhere: `1 self(genAcc) · 2 no-transform · 3 full PT ·
  4 rand/shuffle null · 5 scale-only · 6 rotation-only · 7 translation-only · 8 chance ctrl`.
  **Column 6 (rotation-only) is the headline PT metric.** Chance = 0.02.
  (Exceptions with their own column maps: `_trialShuffledWithMean` = 6 cols,
  `_cross_stimulus_generalization` = 7 cols, `_meancentroid_…` = 6 cols. Headers document each.)
- Six ordered pair files per population: `acec ecac | ecex exec | acex exac`; complementary
  directions combine into three unordered pairs (EC-AC, EC-EX, AC-EX).
- Populations: FR V1 = 48, FR V2 = 112, KO V1 = 109, KO V2 = 146 neurons.
  Sessions: FR V1 = 5, FR V2 = 7, KO V1 = 5, KO V2 = 7. Window `[330 630]` ms, 50 stim × 10 trials.
- The basic pipeline runs `neuron_squence = 1`, so V2 is subsampled to the within-monkey
  V1 count (FR→48, KO→109). Three of four areas sit above n_stim = 50 — the overfitting regime.

---

## Point 1 — Terminology: "cue-invariant" vs "rendering-equivariant"

**No computation required.** Framing/wording change, plus the reviewer's aside that early
visual cortex *cannot* be truly rendering-invariant (or renderings would be indistinguishable)
belongs in the motivation.

If you want to quantify "how far from invariant": the gap between **no-transform transfer
(col 2)** and **rotation-only PT (col 6)** in
`ND/results/decoding_outputs/procrustes_decoding_basic_results/<M>/<A>/` *is* the amount of
transformation required. Already plotted in
`ND/results/figures/{FR,KO} with zscore only rotation.png`.

---

## Point 2 — Test the hypothesis directly with Procrustes-invariant quantities (distances, angles), no decoder

The **decoder-free geometry** work stream, `ND/code/geometry_scripts/`.

### Primary — `geometric_metrics_basic.m`

| | |
|---|---|
| Compute | `ND/code/geometry_scripts/geometric_metrics_basic.m` |
| Outputs | `ND/results/geometric_metrics_outputs/<M>/<A>/geom_{raw,zscore}_results.mat` (all 4 populations) + `geom_pt_results.mat` (**FR V1 and KO V1 only**) |
| Figures | `ND/results/figures/simple_geometrics/mean_of_trials/geometric_metrics_<M>_<A>/` — 20 PNGs per population: `{cos,dist}_{distributions,heatmap_avg,heatmap_rs,scatter_avg,scatter_rs}_{raw,zscore}.png` |

Builds, per rendering, the 50×50 **pairwise Euclidean distance** matrix `D` and the 50×50
**centroid-anchored cosine** matrix `Cos` (cosine of the angle at the population centroid
between stimuli *i* and *j*) from trial-averaged responses — exactly the Procrustes-invariant
quantities the reviewer names. It then correlates those matrices **across renderings**
(1225 unique upper-triangle entries, `triu(...,1)`, diagonal excluded) and brackets the
result with two references computed inside the *same* bootstrap so noise structure matches:
**Null** (stimulus-label shuffle nested in the bootstrap, ≈ 0) and **Ceiling** (two
independent within-rendering bootstrap resamples = split-half reliability of that rendering's
own geometry). Both trial-averaged and trial-resampled estimation modes run every time.

**Result:** cross-rendering correlation is clearly above null but well below ceiling — KO V1
z-scored: AC-EC ≈ 0.38, EC-EX ≈ 0.32, AC-EX ≈ 0.19, against ceilings 0.88–0.91 and null ≈ 0.
Partial, not perfect, preservation of shape — the honest decoder-free version of the claim.
`raw` distance is nearly flat (0.06–0.11); z-scoring is what reveals the structure, and the
decoding pipeline z-scores too, so it is consistent — but say so explicitly.

### Supporting / documented negatives

| Script | Figures | What it shows |
|---|---|---|
| `geometric_metrics_pc1_diagnostic.m` | **none saved** — both `saveas` lines commented (lines 268–269); `.mat` save also commented (line 148) | Per-stimulus PC1 variance fraction (λ₁/Σλ of the centered 10×N cloud) + full scree vs a within-neuron trial-shuffle null. **Observed scree ≈ null scree** → the pseudopopulation's within-cloud covariance is close to independent per-neuron variability. Doubles as a point 7/8 argument. |
| `geometric_metrics_pc1_angle.m` | `simple_geometrics/pc1_angle/<M>_<A>/angle_{within,across}_{raw,zscore}_{full,split}.png` — 8 per population, **all 4 saved**. `.mat` save commented (line 223) | Angle between each stimulus cloud's PC1 and the signal directions to other stimuli. Clean negative: z-score+split medians ~83–84° vs an 85.6° random-orientation null. The `full→split` (~2°) and `raw→zscore` (~6–7°) progressions are the two artifact controls — keep both, they *are* the argument. |
| `geometric_metrics_characteristic_trials.m` | `simple_geometrics/characteristic_trials/<M>_<A>/{cos,dist}_ksweep_{raw,zscore}_{global,own}.png` — 4 each for FR V1/FR V2/KO V2 (`global` only), 8 for KO V1 (`global`+`own`). `.mat` saved for all 4 | Represents each stimulus by the mean of k selected trials (far / near / random / full-mean / null), k = 1…10. Tests whether high-magnitude trials carry the cue-invariant geometry. **Refuted in every panel — far is worst**, ordering near ≥ random > far. Cue-invariant geometry lives in typical, not extreme, responses. |

> **Flags.** (a) `geometric_metrics_basic.m`'s header still says "*saveas calls are commented
> out during active development*" — **stale**: the `.png` saves are live (lines 402/428/475),
> only the `.fig` saves are commented. (b) `geom_pt_results.mat` exists only for FR V1 / KO V1,
> and **no `pt`-mode figures were saved for any population**. (c) `pc1_diagnostic` has no
> figure and no `.mat` — uncomment both if you want to cite the scree comparison.

---

## Point 3 — Does successful transfer decoding actually imply special geometric structure? (the rotation has ~N² parameters)

**The strongest new result in the revision: the PT overfitting artifact.** The rotation has
~N²/2 free parameters against 50 stimulus anchors; when N ≥ 50 it can align two
*structureless* clouds using whatever correspondence it is handed. Four independent lines,
all computed.

### 3a. Synthetic null (Model 1) — the reviewer's own "simplest case"

| | |
|---|---|
| Generator | `SIM/code/generation_scripts/generate_model1_trial_data.m` |
| Driver | `SIM/code/decoding_scripts/procrustes_decoding_simulation.m` |
| Outputs | `SIM/results/decoding_outputs/model1/rho0_N{10,50,100}_ws*_poisson/simulation_results.mat` (6 runs) |
| Figures | **none** — `SIM/results/figures/` is empty and SIM has **no `visualization_scripts/` folder at all** |

Means drawn i.i.d. per (neuron, stimulus, rendering) — no shared structure whatsoever.
**PT ≈ 0.385 at N=100, ≈ 0.029 (chance) at N=10.** Rotation-only ≈ full PT (0.398 vs 0.385),
confirming the ~N²/2-parameter rotation carries the artifact. The existing
shuffled-correspondence null misses it because shuffling maps to *wrong* labels → chance
regardless of overfitting. **The implicit 1/50 chance baseline is too lenient.**

### 3b. Real-data null floors (structureless surrogates at the real N)

| Script (`NDC/code/decoding_scripts/`) | Output leaf | Coverage |
|---|---|---|
| `procrustes_decoding_basic_null_floor_mmtmatched.m` | `NDC/results/decoding_outputs/procrustes_decoding_basic_null_floor_mmtmatched_results/KO/V1/` | **KO V1 only** |
| `procrustes_decoding_basic_null_floor_prmatched.m` | `…/procrustes_decoding_basic_null_floor_prmatched_results/KO/V2/` | **KO V2 only** |

- **Moment-matched:** replaces the 50 stimulus means with per-neuron Gaussian draws
  (destroys structure, preserves per-neuron signal scale), re-adds real residuals.
  KO V1: **PT ≈ 0.91–0.95, no-transform and rand at chance, rotation-only ≈ PT.**
  *Caveat that must be stated:* self-decoding returns ≈ 0.99, above the real level, because
  the surrogate is full-rank — the floor is measured at the wrong operating point.
- **PR-matched** (rigorous successor): applies a random rotation *within the centered stimulus
  subspace* (`Q = Bc·Wr·Bc'`), preserving neuron covariance — hence PR and self-decoding — to
  machine precision while destroying stimulus correspondence. Procrustes fits a *neuron*-space
  rotation and cannot undo a *stimulus*-space one. This is the fair floor; compare real PT
  against it, not against 1/50. Built-in check: its self-decoding must come back at the real
  level, not ~0.99.

### 3c. Correspondence-shuffled control on real data (all 4 populations)

| | |
|---|---|
| Compute | `NDC/code/decoding_scripts/procrustes_decoding_basic_trialShuffledWithMean.m` |
| Viz | `NDC/code/visualization_scripts/visualization_procrustes_decoding_basic_trialShuffledWithMean.m` (`saveas` commented, line 122) |
| Outputs | `NDC/results/decoding_outputs/procrustes_decoding_basic_trialShuffledWithMean/<M>/<A>/` — **all 4 populations** (6-column format) |
| Figures | **none saved** — NDC has no `results/figures/` directory at all |

Self-consistent shuffle: permutation π is applied to the mean matrix, a rotation-only
Procrustes transform is fit from that permuted matrix onto the target, and the **same** π is
applied to the held-out test trials. Scored twice — against target labels (stays HIGH → the
rotation can force an arbitrary correspondence) and against the shuffled identity labels
(≈ 2% → correspondence really was destroyed). The real-data demonstration of the DOF argument.

### 3d. Constraining the DOF directly — PT vs PC dimension

| | |
|---|---|
| Compute | `ND/code/decoding_scripts/procrustes_decoding_basic_incremental_pc.m` — **all 4 populations**, all 6 pair files |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_basic_incremental_pc.m` (3 figures; all `saveas` commented, lines 195/235/289) |
| Support | `ND/code/visualization_scripts/visualization_variance_vs_pc.m` (recomputes the signal spectrum from raw data; `saveas` commented, line 161) |
| Outputs | `ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pc_results/<M>/<A>/` |
| Figures | **none** — `ND/results/figures/incremental_pc/` is empty |

Fixes the full population and sweeps the number of retained **PC dimensions** k, so dimension
is *controlled* rather than measured. Each cue gets its own leakage-clean mean-PCA basis
(stim1's built from training trials only). This is the direct answer to "it will depend on the
ratio between the number of neurons and the number of stimulus classes": at small k the
rotation has few parameters and cannot overfit, so the real/shuffle gap at small k is the
honest signal. Fig 1 = accuracy vs k (raw + normalized by self); Fig 2 = residual Procrustes
shape distance vs k (`.pdist`); Fig 3 = ratio summary. `visualization_variance_vs_pc.m` gives
cumulative signal-variance vs k so you can state how many PCs carry the geometry
(also answers Reviewer 2's minor concern 6).

> **Status gap for point 3: none of 3a–3d has a saved figure.** `SIM/results/figures/`,
> `NDC/results/figures/` (nonexistent), and `ND/results/figures/incremental_pc/` are all empty.
> Null floors cover only **KO V1 (mmt)** and **KO V2 (PR-matched)** — you need the PR-matched
> floor on all four before any quantitative "real PT exceeds the floor by X" claim.

---

## Point 4 — Compare to CCGP (Bernardi et al. 2020): generalization *without* rotating

Two levels, both computed for all four populations.

### 4a. Trial-level CCGP (held-out stimuli, ECOC decoder)

| | |
|---|---|
| Compute | `ND/code/decoding_scripts/procrustes_decoding_cross_stimulus_generalization.m` |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_cross_stimulus_generalization.m` (`saveas` commented, line 145) |
| Outputs | `ND/results/decoding_outputs/procrustes_decoding_cross_stimulus_generalization_results/HoldStim{2,10}/<M>/<A>/` — **both hold-out levels, all 4 populations** |
| Figures | **none saved** — `ND/results/figures/cross_stimulus_generalization/` empty |

The rotation (rotation-only, no scale/translation) is fit on `50 − n_hold` stimulus means and
applied to the `n_hold` stimuli it never saw; everything is scored on the held-out stimuli
with a target-cue-trained ECOC. The CCGP question: *is the cross-cue alignment a global
property of the manifold, or does it need to see every stimulus?*

`.acc` columns: `1 self_decode · 2 no_transform · 3 pt_gen (MAIN) · 4 pt_ceiling (in-sample) ·
5 pt_floor (self-consistent shuffle) · 6 rand_rot · 7 chance`. `pt_gen ≈ pt_ceiling` ⇒ the
rotation generalizes. **Read the header on `pt_floor`:** it is an *overfitting-matched* null,
not a chance floor — at large N it rides toward the ceiling, and the ceiling↔floor gap vs
neuron count is itself the overfitting diagnostic, tying straight back to point 3.
Note `self_decode` is constant across the trial loop within a (partition, neuron) cell —
the viz correctly aggregates per-partition first.

### 4b. Mean-level CCGP (nearest centroid, no SVM, no trial noise)

| | |
|---|---|
| Compute | `ND/code/decoding_scripts/procrustes_meancentroid_cross_stimulus_generalization.m` |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_meancentroid_cross_stimulus_generalization.m` (dated 2026-08-06; `saveas` commented, line 132) |
| Outputs | `ND/results/decoding_outputs/procrustes_meancentroid_cross_stimulus_generalization_results/HoldStim2/<M>/<A>/` — **all 4 populations, HoldStim2 only** |
| Figures | **none saved** — `ND/results/figures/meancentroid_cross_stimulus_generalization/` empty |

Strips the decoder entirely: after applying the fit-on-F rotation, does each held-out source
*mean* land nearest its own target mean? Two reference sets — **GLOBAL** (nearest among all 50
targets, chance 1/50) and **WITHIN** (nearest among the held-out targets only, chance
1/n_hold) — each with three variants: `gen` (main), `no_transform` (how much cue-invariance
already sits in the raw responses — the same neurons are recorded across cues, so this is the
honest baseline), `rand` (floor). `.acc` columns: `1 gen_global · 2 notransform_global ·
3 rand_global · 4 gen_within · 5 notransform_within · 6 rand_within`.

This is the cleanest CCGP analogue: no classifier, no trial noise, and `gen vs no_transform`
is exactly the "does the rotation add anything" comparison Bernardi et al. would ask.
`pt_ceiling`/`pt_floor` are deliberately omitted here — at the mean level with rotation-only
they saturate to ~1 whenever N ≥ n_hold (documented in the header).

> **Status gap:** neither level has saved figures, but **both now have viz scripts** — the
> mean-level viz was added 2026-08-06 (an earlier version of this map said it did not exist).

---

## Point 5 — Benchmark the pipeline on artificial data, varying #neurons, #classes, representation statistics

Entirely in `SIM/`. Decoding math identical to `procrustes_decoding_basic.m` (z-score per cue,
MATLAB `procrustes` similarity transform, ECOC SVM); no neuron subsampling — N is set directly
and the error bar comes from `R_pop` independent synthetic populations (one population = one
complete synthetic experiment).

| Component | File |
|---|---|
| Shared driver (model chosen by `cfg.model_fn`; infers `n_cue`) | `SIM/code/decoding_scripts/procrustes_decoding_simulation.m` |
| Model 1 — null: i.i.d. random means per (neuron, stimulus, cue) | `SIM/code/generation_scripts/generate_model1_trial_data.m` |
| Model 3 — linear mixed selectivity: shared rank-2 shape latents, independent per-cue readout | `SIM/code/generation_scripts/generate_model3_trial_data.m` |
| SNR calibration sweep (model-agnostic) | `SIM/code/generation_scripts/calibrate_weight_scale.m` |

**Runs on disk** (`SIM/results/decoding_outputs/`):
- `model1/rho0_N10_ws0.7_poisson`, `model1/rho0_N50_ws{1,10}_poisson`,
  `model1/rho0_N100_ws{0.7,1,2}_poisson` — the **neuron-count sweep** (N = 10 / 50 / 100) at
  two self-decoding operating points (matched-self ws = 0.7 → self ≈ 0.48; near-ceiling ws = 2
  → self ≈ 0.97).
- `model3/rho0_N100_ws12_poisson` — locked operating point self ≈ 0.49, PR ≈ 1.98.
- `model1/calibration/`, `model3/calibration/` — `calibration_sweep.mat` **and
  `calibration_sweep.png`** (the only saved simulation figures; `saveas` in
  `calibrate_weight_scale.m` line 227 is live).

**Design points for the methods:** Poisson trials (count-faithful version of the reviewer's
Gaussian-mixture request); positivity handled by a **baseline offset, never rectification** —
rectifying curves the manifold and inflates PR, silently turning Model 3 into a nonlinear
model, so the floor must essentially never engage. Models are matched on **self-decoding
(~0.49), not on `weight_scale`**, so PT/PR differences are not SNR confounds; PR is the axis
being measured and is left free. The `rho` knob interpolates independent → shared cross-cue
*readout* (rho = 1 = trivial transfer).

**Reviewer-facing conclusion:** Model 1 is the reviewer's exact "hopefully this can be ruled
out" case — and it is **not** ruled out at N ≥ 50. That is a methods contribution rather than
a retraction, provided you (a) report the null floor and (b) shift claims to relative
comparisons (V1 vs V2, neural vs Gabor at matched PR, real PT vs the PR-matched floor).

> **Gaps.** Models 2/4 are algebraically redundant (skip, but say so). **Model 5** (EX gets an
> extra latent) and **Model 6** (nonlinear mixed selectivity — the one that decouples PR from
> transferability) are **not built**. There is **no stimulus-count (n_stim) sweep** — the
> reviewer explicitly asks to vary the number of stimulus classes and everything is fixed at
> 50. There is **no simulation visualization script** (SIM has no `visualization_scripts/`).

---

## Point 6 — Functional/computational significance of the rotation; how could the brain use it?

Largely a **Discussion** point, but four concrete pieces of computational support already exist:

1. **The rotation lives in very few dimensions**, so a downstream readout is not learning an
   N×N matrix. `ND/results/figures/incremental_pr/fig1_accuracy_pr_vs_N.png` and
   `fig5_saturation_four_views.png` — PR saturates around **8–11** for all populations
   (theoretical cap 49). Plus `visualization_variance_vs_pc.m` for the cumulative signal
   spectrum and `procrustes_decoding_basic_incremental_pc.m` for how much PT survives at small
   k. The readout only has to rotate within a low-dimensional subspace.
2. **The alignment generalizes to stimuli it never saw** — the CCGP scripts under point 4. A
   rotation that must be re-fit per stimulus would be biologically useless; one fit on part of
   the manifold that places the rest correctly is a plausible fixed readout.
3. **The alignment is not tied to specific neurons.**
   `ND/code/decoding_scripts/procrustes_decoding_half_and_half_populations.m` (+ `_incremental`)
   and `procrustes_decoding_neuron_shuffled.m` →
   `ND/results/figures/{FR,KO} half-half decoding with zscore only rotation.png`,
   `{FR,KO} incremental half-half …png`, `{FR,KO} neuron shuffled …png`, `{FR,KO} neuron shuffled.png`.
4. **Hierarchy / feedforward–feedback framing** for "how would a higher area do this":
   `ND/code/decoding_scripts/procrustes_decoding_across_area.m` →
   `ND/results/figures/{FR,KO} across area decoding.png`; ANN analogue
   `ANN/code/decoding_scripts/procrustes_ann_decoding_along_hierarchy.m` →
   `ANN/results/figures/ALEXNET pool 2&5 feedforward and feedback decoding.png`.

**On the z-scoring sub-question:** `geometric_metrics_basic.m` runs `raw`, `zscore` and `pt`
normalization side by side — the `raw` figures are the direct answer to how much survives
without z-scoring. Be candid: for distance, raw cross-rendering correlation is much weaker
(~0.06–0.11 vs ~0.19–0.38 z-scored).

Also relevant: `PT PCA sactter visualization/visualize_procrustesRotation_alignment_pca.m`
(+ `spike_loader.m`) — rotation-only alignment of AC and EC onto EX, then a shared PCA, drawn
as a 2×3 before/after scatter with stimulus thumbnails. No `saveas`; note it resolves data via
`fullfile(pwd, monkey, vp, …)`, so it only runs from inside that folder with a data copy there.

---

## Point 7 — Pseudopopulation limits: missing correlation information; possible block-diagonal rotation structure

### 7a. How much noise correlation is actually there

| | |
|---|---|
| Compute | `ND/code/noise_correlation_scripts/noise_correlation_analysis.m` (random-shuffle control, 100 iterations) |
| Compute | `ND/code/noise_correlation_scripts/noise_correlation_analysis_affine.m` (deterministic affine-shift control, single pass) |
| Outputs | `ND/results/noise_correlation_outputs/<M>/<A>/noise_correlation_results.mat` — all 4. **No affine `.mat`**: the save in the affine script is commented (line 821); the shuffle script's save is also commented now (line 770), so the four `.mat` files are from an earlier run |
| Figures | `ND/results/figures/noise_correlation_<M>_<A>/` (shuffle) and `ND/results/figures/nc_affine_<M>_<A>/` (affine) — **all 4 populations, saved** |

Per neuron pair × stimulus, Pearson r over the 10 trials (isolates noise from signal
correlation), plus the geometric-mean firing rate. Every pair is tagged **same-session vs
cross-session** by rebuilding the neuron→session map from `bad_channel` in the `ac` folder.
Saved figure inventory (prefix `nc_<M>_<A>_` / `nc_affine_<M>_<A>_`): stage-1 example scatters,
geo-mean–NC correlation distributions (original vs shuffle), `mean_vs_variance_nc`
scatterhists, `significance_real_vs_shuffle`, and `delta_mean_nc_by_session`.
Coverage detail: FR V1, FR V2 and KO V2 have both `nanexclude` and `nanzero` variants; **KO V1
has `nanexclude` only**. The affine folders are `nanexclude` only (4 PNGs each).

The affine construction (also used by every trial-perturbation script): 24 permutations
`t → (a·p + b) mod 10` — forward `a=1, b=0..9`; reverse `a=9, b=0..9`; stride-3 `a=3, b=0..3` —
assigned sequentially within a recording session so simultaneously recorded neurons always get
distinct shifts. 236/276 pattern pairs have zero positional overlap, 40 have exactly 2/10.

**The result to report:** same-session pairs carry genuine positive noise correlations,
cross-session pairs sit at zero, and decorrelating collapses same-session down to the
cross-session baseline. Pseudopopulation construction removes most — not all — correlation
structure, and the residual is characterizable.
`nc_<M>_<A>_nanexclude_delta_mean_nc_by_session.png` shows this most directly. Use `nan_mode =
exclude` for the delta figures: in `zero` mode a shuffled constant vector stays constant and
forces delta to 0 artificially.

### 7b. Decoding on *true* simultaneously-recorded populations

| | |
|---|---|
| Compute | `ND/code/noise_correlation_scripts/procrustes_decoding_per_session.m` |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_per_session.m` (6 figures; all `saveas` commented, lines 93–124) |
| Viz | `ND/code/visualization_scripts/visualization_procrustes_decoding_per_session_distribution_test.m` (per-cell Mann-Whitney AUC; prints tables) |
| Outputs | `ND/results/decoding_outputs/procrustes_decoding_per_session_results/{none,affine}/<M>/<A>/` — **all 4 populations, both modes** |
| Figures | **none saved** |

Slices the pseudopopulation back into real single-session populations (5 / 7 sessions,
~5–15 neurons each) via `neuron_session_id`, no neuron downsampling, and runs the same PT
decoding with and without the affine perturbation. Figures: per-session decoding per mode,
delta, self-vs-PT scatter, normalized delta with Wilcoxon + FDR, and a pooled figure across
both monkeys.

**The agreed framing ("Option C"):** the perturbation's effect is **not unidirectional** —
some sessions/pairs show affine < none (correlations were aiding decoding), others affine >
none (correlations were limiting it), matching the information-limiting vs signal-correlation
dichotomy. **Do not claim a single signed effect.** Pooled across stimulus pairs (N = 15 for
V1, 21 for V2 per monkey), nothing survives FDR over 12 tests — closest is KO V2 PT(rot),
raw p = 0.0064, q = 0.076. Per-monkey per-pair tests are hopelessly underpowered (min
two-sided signrank p = 2/2^N: N=5 → 0.0625). State this explicitly rather than hiding it.
Also state the pseudoreplication caveat: the 50–100 trial-resamples within a session are not
independent biological observations, which is why effect sizes are aggregated at session level.

### 7c. Shared helpers

`apply_trial_perturbation.m` (modes `none | shuffle | affine`; reorders trials within each
(neuron, stimulus) block, preserving marginals exactly) and `build_affine_shift_config.m`
(builds the 24 patterns + per-neuron session/pattern assignment). Present in **both**
`ND/code/noise_correlation_scripts/` and `ND/code/geometry_scripts/` so each stack is
self-contained — keep the two copies in sync if either is edited.

> **Gap on the reviewer's specific example:** he suggests you might find **block-diagonal
> structure in the fitted rotation matrices** corresponding to simultaneously recorded sets.
> **Nothing in the repo inspects the rotation matrix `T` itself** — every script discards it
> after use. Small, high-value new analysis: save `T` from `pro_decoding`, order neurons by
> session, compare within-session vs cross-session block energy against a session-label shuffle.

---

## Point 8 — Noise correlations shape real geometry; simulate injecting structured NC (mean pairwise r = 0.05–0.15)

Same machinery as point 7, applied to the decoding and geometry pipelines.

| Analysis | Compute | Outputs / figures |
|---|---|---|
| PT decoding with NC broken (pseudopopulation) | `ND/code/noise_correlation_scripts/procrustes_decoding_trial_perturbed.m` | `ND/results/decoding_outputs/procrustes_decoding_trial_perturbed_results/affine/<M>/<A>/` — all 4, **`affine` mode only** (no `none` leaf; the `none` baseline is `procrustes_decoding_basic_results`). Viz: `ND/code/visualization_scripts/visualization_procrustes_decoding_trial_perturbed.m` — **no figures saved** |
| Incremental PT + PR with NC broken | `ND/code/noise_correlation_scripts/procrustes_decoding_basic_incremental_pr_trial_perturbed.m` | `ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pr_trial_perturbed_results/affine/<M>/<A>/` — all 4. **No dedicated viz script** |
| Geometry with vs without NC | `ND/code/geometry_scripts/geometric_metrics_basic_nc_compare.m` | `ND/results/figures/simple_geometrics/nc_control/<M>_<A>/{cos,dist}_distributions_nc_{raw,zscore}.png` — **all 4 populations, both norms, saved.** `.mat` save commented (line 112) |
| Within-cloud covariance structure | `geometric_metrics_pc1_diagnostic.m`, `geometric_metrics_pc1_angle.m` | see point 2 |

**The clean mathematical argument on the geometry side** (`geometric_metrics_basic_nc_compare.m`):
the affine shift only reorders trials *within* each (neuron, stimulus) block, so the block
mean — and therefore the whole trial-averaged geometry, its cross-rendering correlations, and
the black-triangle point estimates — is **invariant to machine precision**. Z-scoring is
permutation-invariant too. The figures confirm it: the trial-averaged triangles coincide
exactly. Where NC *does* appear is the **bootstrap spread**: removing NC makes per-neuron mean
estimates fluctuate more independently, so the **within-condition ceiling rises** (KO V1
z-scored distance: AC self 0.82 → 0.87, EC self 0.90 → 0.95) while the observed cross-rendering
correlation barely moves (0.26 → 0.27). Both versions share one RNG seed, so each category is
compared by a paired Wilcoxon signed-rank on matched bootstrap iterations (stars on the figures).

Note the same invariance is a built-in **sanity check** in the incremental-PR-perturbed script:
PR is computed on the trial-averaged manifold, so it must be unchanged by the perturbation —
any PR difference at matched neuron count indicates a bug, not an NC effect. Its header says so.

**Headline for points 7/8:** *cue-invariant geometry is a trial-averaged quantity and is
mathematically invariant to noise correlations; NC affects estimation reliability, not the
geometry. Single-trial decoding is where NC can matter, and there the effect is session- and
pair-dependent, not unidirectional.*

> **Gap — the reviewer's literal ask is not built.** He asks for a **simulation injecting
> structured noise correlations of realistic magnitude (r = 0.05–0.15) into the
> pseudopopulation** and testing whether PT is robust or systematically biased. What exists is
> the *opposite* operation (removing residual NC). The `rho` knob in the simulation generators
> is cross-cue *readout weight* correlation, not trial-to-trial NC. Adding a correlated-noise
> term to `generate_model1/3_trial_data.m` (a shared per-trial latent with loading tuned to a
> target mean pairwise r) and re-running the driver at r = 0, 0.05, 0.10, 0.15 would answer it
> directly and is a modest addition. `noise_correlation_analysis.m` already supplies the
> **empirical** r distributions needed to justify (or contest) the 0.05–0.15 range.

---

## Point 9 — Neural vs Gabor at matched population size, with participation ratio

The most complete work stream — **compute, visualization, and saved figures all exist.**

| | |
|---|---|
| Neural compute | `ND/code/decoding_scripts/procrustes_decoding_basic_incremental_pr.m` → `ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pr_results/<M>/<A>/` (all 4) |
| Gabor compute | `GAB/code/decoding_scripts/procrustes_Gabor_resps_decoding_basic_incremental_pr.m` → `GAB/results/decoding_outputs/procrustes_Gabor_incremental_pr_results/even and odd combined/` (**this filter model only**) |
| Gabor trial data | `GAB/code/decoding_scripts/generate_gabor_trial_data.m` (complex = Poisson on filter energy; odd/even = half-wave rectified then Poisson; combined = concatenation) |
| Viz (Figs 1–5) | `ND/code/visualization_scripts/visualization_procrustes_decoding_basic_incremental_pr.m` |
| **Figures (saved)** | `ND/results/figures/incremental_pr/fig1_accuracy_pr_vs_N.png`, `fig2_accuracy_vs_pr.png`, `fig3_normalized_transfer_vs_pr.png`, `fig4_transfer_per_dim_vs_N.png`, `fig5_saturation_four_views.png` |

**PR definition** (state it in methods): `PR = (Σλ)²/Σλ²` over the covariance spectrum via SVD
of the **z-scored, 10-trial-averaged, stimulus-centered 50×N signal manifold**, computed once
per neuron-sampling repeat, no second renormalization after averaging. PR is invariant to the
Procrustes transform, so it is a clean per-condition quantity. Uniform 10-unit neuron sweeps so
neural and Gabor population sizes match exactly. **Gabor PR uses the mean estimated from the
generated Poisson trials, not the noiseless filter rate** — otherwise Gabor gets an unfair
clean-mean advantage. To pool a condition's PR, gather `pr_stim1` from pairs where stim1 is
that condition and `pr_stim2` from pairs where stim2 is (each condition appears in 4 of 6 pairs).

What each figure says:
- **Fig 1 (accuracy + PR vs N):** at matched *neuron count*, Gabor wins on both axes
  (self ≈ 1.0, PR ≈ 11). Unfavourable and **not** the right comparison — Gabor is near-noiseless.
  Show it, then explain the move to matched PR.
- **Fig 2 (accuracy vs PR — the headline):** **KO V1 sits clearly above the Gabor curve**
  (higher rotation transfer at equal PR); FR V2 ≈ on it; FR V1 slightly below; **KO V2 below.**
  Population-specific, not a blanket neural win — say so.
- **Figs 3–4 (normalized):** PT/self vs PR and transfer-per-dimension vs N; ratio-of-pooled-means
  with a percentile bootstrap over neuron-repeats (the resampling unit). Fig 3 also overlays the
  correspondence-shuffled control (col 4) against true PT (col 6) — the gap isolates the part of
  transfer that requires true stimulus correspondence, which links this figure to point 3.
- **Fig 5 (2×2 saturation):** PR(N) asymptotes are well constrained under `satexp`
  (PR_∞ ≈ FR V1 8.2, FR V2 ~8.5, KO V1 8.5, KO V2 10.3, Gabor 11.2). **PT(N) does not converge** —
  the ~0.5 ceiling is a fit projection from extrapolation, not an observed plateau, and is
  fit-form dependent (`fit_method` switches `satexp` ↔ `mm`). Present it as a caveat.

**Lead with the dissociation:** Gabor dominates *self*-decoding per dimension but loses
*transfer* per dimension to KO V1. Where neural geometry is efficient, the efficiency is
specific to transferable/cue-invariant structure, not to general decodability.

Further baselines available: `GAB/results/decoding_outputs/{complex, even and odd combined,
even and odd combined downsampled}/` with `GAB/results/figures/*only-rotation-PT decoding.png`
(2 saved); and the full ANN suite — `ANN/results/decoding_outputs/procrustes_ann_decoding_
{basic,basic_incremental}_results/{alexnet,vgg}/pool*/` and `…_along_hierarchy_results/`, with
5 saved figures in `ANN/results/figures/`.

> **Optional additions discussed but not built:** the `complex` Gabor model in the matched-PR
> figures; a per-pair (AC-EC / EC-EX / AC-EX) breakout; reporting the `satexp` vs `mm`
> asymptote spread as an explicit uncertainty band.

---

## Point 10 — Compare Procrustes residuals across the three rendering pairs (is AC-EX alignment tighter than genuine cue-invariance predicts?)

The concern: AC is the first PC of the same patches that define EX, so AC-EX alignment may
reflect shared stimulus construction rather than cue-invariant coding.

**Decoder-free, already per-pair — the best existing evidence.** `geometric_metrics_basic.m`
(point 2) computes the cross-rendering correlation separately for `acec`, `acex`, `ecex`, each
with its own null and its own per-condition ceilings. In KO V1 z-scored the ordering is
**AC-EC ≈ 0.38 > EC-EX ≈ 0.32 > AC-EX ≈ 0.19** — AC-EX is the **weakest**, not the tightest,
which is the direct empirical answer to the worry. **Check that the ordering holds in the other
three populations before committing to it in writing** — the `geom_zscore_results.mat` files
for all four are on disk, so this is a read, not a re-run.

**Procrustes residuals proper.** The MATLAB `procrustes` residual `d` is computed inside almost
every decoding script but **discarded** — the only place it is stored is `.pdist` in
`ND/code/decoding_scripts/procrustes_decoding_basic_incremental_pc.m` →
`ND/results/decoding_outputs/procrustes_decoding_basic_incremental_pc_results/<M>/<A>/`
(all 4 populations, all 6 pair files). Columns: `1 d_before` (no rotation), `2 d_rot`
(rotation only), `3 d_full` (optimal b, T, c) — all sharing MATLAB's stim2-target denominator,
so they are comparable across pairs. Caveat from the header: these are *residual* distances
after each cue's own top-k PCA basis, not the full cross-cue rotation magnitude.

**Per-pair PT decoding** for the same comparison at the accuracy level:
`ND/results/decoding_outputs/procrustes_decoding_basic_results/<M>/<A>/` with
`ND/code/visualization_scripts/visualization_procrustes_decoding_basic.m` →
`ND/results/figures/{FR,KO} with zscore only rotation.png` (already broken out by pair); and
the three-panel per-pair layout in both CCGP viz scripts.

> **Gap:** `visualization_procrustes_decoding_basic_incremental_pc.m` **pools `.pdist` across
> all six ordered pairs**. The data to answer point 10 is already on disk — un-pool it into the
> three unordered pairs (`acec+ecac`, `ecex+exec`, `acex+exac`) and compare `d_full` / `d_rot`
> across them. Small edit to the existing plotting loop, not a new analysis. There is still
> **no dedicated residual-across-pairs script or statistic**, and no full-population (non-PCA)
> residual is saved anywhere.

---

## Cross-cutting status summary

### (a) Computed **and** visualized, figures on disk — cite directly

| Point | Figures |
|---|---|
| 2 | `simple_geometrics/mean_of_trials/geometric_metrics_<M>_<A>/` (20 each, all 4) · `simple_geometrics/pc1_angle/<M>_<A>/` (8 each, all 4) · `simple_geometrics/characteristic_trials/<M>_<A>/` (4 each; 8 for KO V1) |
| 6 | `{FR,KO} half-half …`, `{FR,KO} incremental half-half …`, `{FR,KO} neuron shuffled …`, `{FR,KO} across area decoding.png` · `ANN/results/figures/` (5) |
| 7a | `noise_correlation_<M>_<A>/` (19–26 each, all 4) · `nc_affine_<M>_<A>/` (4 each, all 4) |
| 8 | `simple_geometrics/nc_control/<M>_<A>/` (4 each, all 4) |
| 9 | `incremental_pr/fig1…fig5` · `GAB/results/figures/` (2) |
| 5 | `SIM/results/decoding_outputs/model{1,3}/calibration/calibration_sweep.png` (only saved sim figures) |
| 1/10 | `{FR,KO} with zscore only rotation.png` (per-pair PT bars) |

### (b) Computed but figures not saved — run the viz before writing

All of these have `.mat` outputs on disk and a working viz script; only `saveas` is commented
(or, for NDC, the figures directory does not exist yet).

| Point | Analysis | Empty/absent figure target |
|---|---|---|
| 3a / 5 | simulation driver runs (Models 1, 3) | `SIM/results/figures/` — **and no viz script exists** |
| 3b | both null floors | NDC has **no `results/figures/`** dir; no viz wired to these leaves |
| 3c | `trialShuffledWithMean` | NDC, `saveas` commented |
| 3d | PT vs PC dimension (3 figs) + `variance_vs_pc` | `ND/results/figures/incremental_pc/` (empty) |
| 4a | trial-level CCGP, HoldStim2 + HoldStim10 | `ND/results/figures/cross_stimulus_generalization/` (empty) |
| 4b | mean-level CCGP, HoldStim2 | `ND/results/figures/meancentroid_cross_stimulus_generalization/` (empty) |
| 7b | per-session decoding (6 figs) + distribution test | `saveas` commented (lines 93–124) |
| 8 | trial-perturbed PT decoding | `saveas` absent in that viz script |
| 2 | `pc1_diagnostic` scree | `saveas` commented (lines 268–269) |
| 8 | incremental-PR trial-perturbed | **no viz script** (the unperturbed PR viz would need a path switch) |

Also unsaved as `.mat`: `nc_compare_*`, `pc1_angle_*`, `pc1_diagnostic_*`
(`geometric_metrics_outputs/`), and `noise_correlation_results_affine.mat`
(`noise_correlation_outputs/`) — the figures exist but the underlying numbers are not persisted,
so any number quoted from them currently requires a re-run.

### (c) Genuine analysis gaps — the reviewer's ask is not built

1. **Point 8** — inject structured NC (mean pairwise r = 0.05–0.15) into the simulation and
   test PT robustness. The reviewer's literal request; nothing built. Only NC *removal* exists.
2. **Point 3/5** — PR-matched null floor run on **KO V2 only**; moment-matched on **KO V1 only**.
   Needs the PR-matched floor on all four populations before any quantitative floor claim.
3. **Point 5** — no **stimulus-count** sweep (n_stim fixed at 50 everywhere, but the reviewer
   asks explicitly for varying the number of classes). Models 5 and 6 not built; Model 6 is the
   one that decouples PR from transferability.
4. **Point 7** — no inspection of the fitted rotation matrix `T` for block-diagonal /
   session structure. `T` is computed everywhere and saved nowhere.
5. **Point 10** — no per-pair Procrustes-residual statistic. Data exists (`.pdist`), needs
   un-pooling; no full-population residual is stored at all.
6. **Point 6** — the "how is it implemented in a network" question has no computational
   component and is not planned as one; it is Discussion-only.

### (d) Stale paths and broken copies — fix before running

- **`NDC/code/visualization_scripts/visualization_procrustes_decoding_basic.m`** points at
  `procrustes_decoding_basic_null_floor_results`, which **does not exist**. The real leaves are
  `..._null_floor_mmtmatched_results` and `..._null_floor_prmatched_results`. It also has V2
  loads commented out (set up for the KO V1 mmt run only).
- **`NDC/code/visualization_scripts/visualization_procrustes_decoding_basic_incremental_pr.m`**
  resolves `neural_base` to `NDC/results/decoding_outputs/procrustes_decoding_basic_incremental_pr_results`,
  which does not exist under NDC — it would warn-and-skip all four neural populations and plot
  Gabor alone. This copy has diverged from the ND original: fig3/fig4/fig5 `saveas` are **live**
  here (writing into a not-yet-existing `NDC/results/figures/incremental_pr/`) while the ND copy
  has all five commented. Prefer the ND copy; the NDC copy needs its base path repointed or deleting.
- **Most other `NDC/code/visualization_scripts/` files** (`…_across_area`, `…_basic_incremental`,
  `…_half_and_half_populations`, `…_half_populations_incremental`, `…_neuron_shuffled`,
  `…_per_session`, `…_trial_perturbed`) are copies of the ND scripts pointing at
  `NDC/results/decoding_outputs/<leaf>` — none of those leaves exist under NDC. Only
  `visualization_procrustes_decoding_basic_trialShuffledWithMean.m` targets a real NDC leaf.
- **`geometric_metrics_basic.m` header** claims `saveas` is commented out; the `.png` saves are
  in fact live (only `.fig` is commented).
- **`visualization_procrustes_decoding_basic_incremental_pr.m` header** still says "Figure 2
  (accuracy vs PR) will be appended here" — all five figures are present.
- **`PT PCA sactter visualization/`** resolves data with `fullfile(pwd, …)`, so it only runs
  from inside that folder with a copy of the data tree there. (Folder name has a typo: "sactter".)
- Duplicated helpers to keep in sync: `multiclass_svmloader_PT.m` (4 copies),
  `apply_trial_perturbation.m` and `build_affine_shift_config.m` (2 copies each).
