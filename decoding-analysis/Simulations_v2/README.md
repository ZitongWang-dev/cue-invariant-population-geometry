# Simulations_v2 — PT transfer-decoding benchmark with the corrected control

Second version of the synthetic benchmark for reviewer points 3 and 5 (Reviewer 1,
J Neurosci). It replaces `Simulations/` (v1), which is kept untouched for reference.
Everything here is MATLAB (R2025b, Statistics and Parallel Computing toolboxes).

## What changed relative to `Simulations/` (v1)

| | v1 (`Simulations/`) | v2 (this folder) |
|---|---|---|
| Control | "rand null": transform fit from *shuffled* means, applied to *unshuffled* test trials → always chance, blind to overfitting | self-consistent correspondence shuffle of the real-data pipeline: one permutation π applied to **both** the stim1 training means and the held-out stim1 trials, rotation-only fit, scored against target labels (= the overfitting floor). Same as `neuronal_decoding_control/..._trialShuffledWithMean.m`, `neuronal_decoding/..._incremental_pr.m`, `..._incremental_pc.m`. Its identity-label score (column 9) must sit at chance. |
| Sweep axes | N ∈ {10, 50, 100}, n_stim fixed at 50, one N per run | N ∈ {5, 10, 20, …, 100} × n_stim ∈ {5, 10, 20, 30, 40, 50}, both models, one run, resumable |
| Means (Model 1) | Gaussian deviations + large constant positivity offset (shallow tuning, count set indirectly) | Gamma-distributed random means with each neuron's mean count and tuning depth drawn from the **measured** distributions (`calibrate_from_real_data.m`); positive by construction, no offset, no floor |
| Means (Model 3) | 10×5 latent grid (50 stimuli only), Gaussian weights, offset ≈ 5.5·ws·M_grid → count ≈ 140 needed for data-like self-decoding | random whitened 2-D latents (any n_stim), readout weights scaled so each neuron's tuning depth is exactly the same as in Model 1 with the same seed; positivity guaranteed by bounding the depth, never by rectification |
| Operating point | `weight_scale` searched for a target self-decoding | single-neuron statistics fixed at the empirical values (median count 9.4, log-SD 0.66; tuning CV 0.14, log-SD 0.65; Fano ≈ 1.08 → Poisson trials); self-decoding is then a *result* of geometry, N and n_stim |
| Visualization | none | `visualization_scripts/visualization_simulation_sweep.m` |
| Output | one `simulation_results.mat` per (model, N, weight_scale) folder | **one structured file per sweep** (`results/decoding_outputs/<tag>_results.mat`): both models, every grid point, every population, plus a summary table |

Same as v1: z-score per cue, MATLAB `procrustes` similarity fit on trial-averaged
means (stim2: 10 trials, stim1: 9 training trials), one held-out trial per stimulus,
ECOC linear SVM trained on all stim2 trials, self-decoding by 10-fold CV, all N
synthetic neurons used (no subsampling), 2 cues → 2 ordered pairs, error bars from
independent synthetic populations.

## Layout

```
Simulations_v2/
  code/
    generation_scripts/
      calibrate_from_real_data.m      measures count / tuning-depth / Fano of the 4 recorded populations
      generate_model1_trial_data.m    Model 1: independent random means (null)
      generate_model3_trial_data.m    Model 3: shared D-dim latents (n_latent, default 2), independent per-cue linear readout
    decoding_scripts/
      pt_transfer_decode_population.m          one population, all pairs, corrected control (+ optional PC constraint)
      procrustes_decoding_simulation_sweep.m   driver: models × N × n_stim → one structured results file, resumable
    visualization_scripts/
      visualization_simulation_sweep.m         heatmaps, curves vs N, N/n_stim comparison, real-vs-model overlay, text summary
      figure_results_section.m                 paper-style Figures 1-2 + every number quoted in the results write-up
  results/
    calibration/empirical_operating_point.mat
    decoding_outputs/<tag>_results.mat         ALL results of a sweep (see "Results file" below)
    decoding_outputs/<tag>_run.log             console log of the batch run
      sweep_N_S                                first run: Model 1 + Model 3 with a 2-D latent, 5-50 classes
      sweep_N_S_model3_D5                      Model 3 with a 5-D latent, 10-50 classes, same neurons as Model 1
    figures/calibration_single_neuron_stats.png, <tag>/*.png, <tag>/summary_<tag>.txt
    figures/results_section/                   fig1_synthetic_benchmark.png, fig2_recorded_vs_models.png, results_section_numbers.txt
```

## Running

1. (Optional, already done) `calibrate_from_real_data.m` — prints and saves the empirical operating point.
2. `procrustes_decoding_simulation_sweep.m` — edit the CONFIGURATION block, or set
   `cfg_override = struct(...)` in the workspace first for a quick test, e.g.
   `cfg_override = struct('R_pop',2,'N_list',[10 100],'S_list',[5 50],'trial_repeat',2,'sweep_tag','smoke')`.
   The results file is rewritten after every grid point (temporary file, then replace), so an
   interrupted run loses at most one point: rerunning the same configuration resumes and skips the
   points already marked done. A *different* configuration under the same `sweep_tag` stops with an
   error instead of mixing results; set `cfg.overwrite = true` or pick a new tag.
   Cost is dominated by the number of classes (n_stim = 50 ≈ 4 min per point with 16 populations on 8 workers).
3. `visualization_simulation_sweep.m` — reads the results file (works on a partial run).
   Plot another sweep with `viz_override = struct('tag', '<tag>')` in the workspace.

The Model 3 run with a 5-D latent was launched from `code/decoding_scripts` with:

```matlab
cfg_override = struct('models', {{'model3'}}, 'model_fns', {{@generate_model3_trial_data}}, ...
    'model_descriptions', {{'linear mixed selectivity (shared 5-D latent geometry, independent per-cue readout)'}}, ...
    'S_list', [10 20 30 40 50], 'seed_S_list', [5 10 20 30 40 50], ...
    'gen', struct('n_latent', 5, 'cap_rule', 'exact'), 'sweep_tag', 'sweep_N_S_model3_D5');
procrustes_decoding_simulation_sweep
```

The 5-class row is dropped because a 5-D latent needs at least 6 stimuli. `seed_S_list` numbers the
seeds on the first run's grid, and running Model 3 alone uses the first model's seed block, so every
Model 3 population has exactly the same neurons (rates and tuning depths) as the Model 1 population at
the same grid point. `cap_rule = 'exact'` caps each neuron's tuning depth at the largest value that keeps
its own rates non-negative (and at 0.35, as in Model 1). The first run used the older `'radius'` rule,
which reproduces its 2-D Model 3 results exactly but caps too many neurons once the latent has more
dimensions.

Useful variants (set in the driver): `cfg.pc_dims = k` projects each cue into the top-k of its own
mean-PCA basis before decoding (the `incremental_pc` construction; rotation DOF = k(k−1)/2);
`cfg.gen.cv_median = …` raises or lowers tuning depth for a matched-self-decoding variant;
`cfg.gen.rho` interpolates independent → shared cross-cue means/readouts.

## Output columns (per population, per ordered pair)

`1 self · 2 no-transform · 3 PT full (b·X·T + c) · 4 control · 5 scale-only · 6 PT rotation-only · 7 translation-only · 8 chance · 9 control vs identity labels`

Column 6 is the headline PT and column 4 its matched control (both rotation-only). The reported
ratio is PT/control computed per population, then averaged. The results file also stores the
noise-free signal PR (generator), the empirical PR of the z-scored condition means (the real-data
definition), mean count, and the in-sample Procrustes residuals of the true and shuffled fits.

## Results file

`results/decoding_outputs/<tag>_results.mat` holds one variable, `sweep`:

| Field | Content |
|---|---|
| `info`, `created`, `updated` | how to read the file; timestamps |
| `cfg` | run configuration (generator names instead of function handles) |
| `models`, `N_list`, `S_list` | model names; neuron counts = dim 1 and stimulus counts = dim 2 of every grid array |
| `metric_names`, `pair_labels` | the 9 accuracy columns above; ordered cue pairs (`'ac->ec'`, `'ec->ac'`) |
| `summary` | table, one row per model × N × n_stim: mean and SEM of every metric, PT/control ratio, PR, mean count, `done` |
| `model1`, `model3` | per model: `done` [nN×nS]; `mean`, `sem` [nN×nS×9]; `ratio_mean`, `ratio_sem`, `ratio_of_means`, `PR_true`, `PR_emp`, `count`, `k_eff`, `runtime_s` [nN×nS]; `pdist_mean` [nN×nS×2] |
| `model*.per_pop` | every population: `seeds`, `ratio`, `mean_count` [nN×nS×R]; `pop_metrics` [nN×nS×R×9]; `pair_metrics` [nN×nS×R×pair×9]; `pair_pdist` [nN×nS×R×pair×2]; `signal_PR`, `PR_emp` [nN×nS×R×cue] |

```matlab
S = load('results/decoding_outputs/sweep_N_S_results.mat');  sweep = S.sweep;
sweep.summary                                          % browse everything as one table
iN = find(sweep.N_list == 100);  iS = find(sweep.S_list == 50);
sweep.model1.mean(iN, iS, 6)                           % PT (rotation-only), N = 100, 50 stimuli
squeeze(sweep.model1.per_pop.pop_metrics(iN, iS, :, 6))    % the 16 populations behind it
```

## Reading the results

* Model 1 (no shared structure): PT should track the control everywhere (ratio ≈ 1). Wherever PT is
  high, the control is equally high — that is the rotation's overfitting, now measured rather than
  hidden by a chance baseline.
* Model 3 (shared rank-2 geometry): PT should exceed the control (ratio > 1); the margin is largest
  where the rotation has few degrees of freedom relative to the number of classes.
* `real_vs_models_nstim50.png` overlays the recorded populations' incremental_pr curves
  (self, PT, control vs N, 50 stimuli, same decoder and control) on the two models' n_stim = 50 rows.

## Results of the first full run (`sweep_N_S`; 16 populations per point, 10 holdouts, 10-fold CV)

Paper-style figures and the write-up numbers come from `code/visualization_scripts/figure_results_section.m`,
which writes `results/figures/results_section/` (Figure 1, Figure 2, and `results_section_numbers.txt` with
every number quoted below). All data: `results/decoding_outputs/sweep_N_S_results.mat`. Ratios below are
ratios of mean accuracies (mean PT / mean control).

**The corrected control works.** For Model 1 (no shared structure) PT/control has a median of 1.00 over the
66 grid points (0.87–1.15; 0.95–1.05 for N ≥ 20). At 50 classes, PT rises from 0.03 (N = 5) to 0.30 (N = 50)
and 0.75 (N = 100), and the control from 0.03 to 0.31 and 0.74. Scored against the permuted identities, the
control stays within 0.02 of chance, so the correspondence really is destroyed.

**Under the null, the floor is set by neurons per class.** Model 1's control/self is 0.28–0.50 at
N/n_stim ≤ 0.25, 0.49–0.59 at N = n_stim for 20–50 classes (0.72–0.84 for 5–10 classes), 0.75–0.91 at
N = 2·n_stim, and 0.97–1.00 from N/n_stim = 6 upward. The floor exceeds chance even when N < n_stim
(N = 10 with 50 classes: 0.04 against 0.02).

**Shared geometry is detected above the floor.** Model 3's PT matches its self-decoding (PT/self 0.85–1.07)
while its control stays low (N = 100, 50 classes: self 0.17, PT 0.17, control 0.05). PT/control grows with the
number of classes (means over N: 1.23, 1.73, 2.36, 2.83, 2.97, 3.39 for 5–50 classes) and changes little with
N beyond 20 (3.30–3.93 at 50 classes). Its floor does not rise with N (control/self 0.26–0.36 at 50 classes),
and its participation ratio stays at 2.3–3.5, whereas Model 1's grows from 3.9 to 30.

**Recorded populations at 50 classes lie between the two models** (incremental_pr results, same decoder and
control). At full size: self 0.26–0.49, PT 0.11–0.26, control 0.09–0.17, PT/control 1.26–1.50 (1.15–1.57
across subsampled N, flat in N), participation ratio 7.0–9.5. Model 1 reaches the recorded self-decoding
range at N = 20–40 with PT/control 0.98–1.05, so the recorded margin is unlikely to be a signal-to-noise effect.

Operating-point caveat: single-neuron statistics were matched to the recordings but self-decoding was left
free, so Model 1 self-decodes above and Model 3 below the recorded populations at the same N. A matched-self
Model 3 (deeper tuning via `cfg.gen.cv_median`, or more latents via `cfg.gen.n_latent`) has not been run.

## Model 3 with a 5-D latent (`sweep_N_S_model3_D5`, 2026-10-07; 108 min)

This is now the headline Model 3: the paper figures and the write-up take Model 3 from this run and
Model 1 from `sweep_N_S`. Numbers below come from `results/figures/results_section/results_section_numbers.txt`.

- **Transfer.** PT matches self-decoding from 30 neurons upward (PT/self 0.92–1.00) and reaches 0.62–0.80
  with 5–10 neurons, where each cue's random readout distorts the 5-D geometry differently.
- **Margin over the floor.** PT/control grows with the number of classes: means over N of 1.33, 2.10, 2.83,
  3.40 and 4.01 for 10–50 classes, and 3.59–4.91 at 50 classes for N ≥ 20. With 10 classes it is only
  1.31–1.36, because ten points in five dimensions leave the rotation room to fit an arbitrary pairing.
  The floor does not rise with N (control/self 0.20–0.30 at 50 classes).
- **Matched operating point.** At 50 classes and N = 40–100, self-decoding (0.28–0.41) and participation
  ratio (7.1–7.8) fall within the recorded ranges (0.26–0.49 and 7.0–9.5), while PT/control is 4.3–4.9
  against 1.26–1.50 in the recordings. A fully shared geometry at the recordings' operating point would
  therefore give a margin about three times larger than observed.
