# Simulations_v3_partial_sharing — PT transfer decoding with partially shared structure

Third simulation folder for reviewer points 3 and 5 (Reviewer 1, J Neurosci). It adds the
reviewer's "more interesting" family: mean responses built from contour-shape parameters that are
**shared** across renderings plus appearance parameters that are **specific** to each rendering.
`Simulations_v2/` (Model 1 null, fully shared Model 3) is unchanged and supplies the decoding
core and the α = 1 endpoint. MATLAB R2025b, Statistics and Parallel Computing toolboxes.

## The model family (fixed-K block design)

Every cue has the same number K of latent dimensions; K_s of them are shared across cues and the
other K − K_s are drawn independently per cue:

    alpha = K_s / K          (K = 5:  alpha = 0, 0.2, 0.4, 0.6, 0.8, 1)

Each cue's latent is whitened as a whole, so **every within-cue quantity has the same
distribution at every α** — each neuron's mean count and tuning depth, the signal spectrum, the
participation ratio, and therefore self-decoding. Only the relation between cues changes. This
is deliberate: a √α / √(1−α) mixture of a shared and a cue-specific subspace would change the
per-cue PR with α (PR = 1/(α²/D_s + (1−α)²/d_c), e.g. 5 → 10 → 5 for 5 + 5 dimensions), and the
correspondence-shuffled floor depends on dimensionality, so "more sharing" and "a different
floor" would be confounded.

- **α = 1** is `Simulations_v2` Model 3 with a 5-D latent (`sweep_N_S_model3_D5`), reproduced bit
  for bit (same seeds, same random-number stream).
- **α = 0** is a null with no shared structure *at the recorded dimensionality* (PR ≈ 7–8 at
  50 classes), unlike Model 1 (PR ≈ 30 at N = 100).
- K = 5 is the latent size at which the fully shared model matches the recordings' self-decoding
  and PR at 50 classes (see `Simulations_v2/README.md`).

Generative model (neuron i, stimulus s, cue c): `lambda = r_i + r_i·cv_i·u_{i,c}·z_{s,c}`,
Poisson trials. r_i and cv_i are drawn from the measured distributions (median count 9.4,
tuning CV 0.14), u_{i,c} is a random unit readout direction per neuron and cue, and
z_{s,c} = whiten([x_s(1:K_s), y_{s,c}]) with y independent across cues. Cholesky whitening keeps the
first K_s columns identical in every cue and orthogonalizes the cue-specific columns against them.
Positivity comes from capping each neuron's tuning depth (no offset, no rectification).

## Layout

```
Simulations_v3_partial_sharing/
  code/
    generation_scripts/
      generate_partial_sharing_trial_data.m       K latent dims per cue, K_s shared (alpha = K_s/K)
    decoding_scripts/
      procrustes_decoding_partial_sharing_sweep.m driver: alpha x N x n_stim -> one structured results file
      real_data_stimulus_subsampling.m            recorded populations on random subsets of C stimuli
    visualization_scripts/
      visualization_partial_sharing_sweep.m       exploratory figures + text summary (as v2's visualization_simulation_sweep.m)
      figure_partial_sharing_results.m            paper-style Figures 1-2 + every quoted number (as v2's figure_results_section.m)
  results/
    decoding_outputs/sweep_alpha_N_S_results.mat  ALL results of the alpha sweep (variable `sweep`)
    decoding_outputs/real_stimsub_results.mat     ALL results of the recorded stimulus subsampling (variable `real`)
    decoding_outputs/*_run.log                    console logs of the batch runs
    figures/sweep_alpha_N_S/                      exploratory figures, summary_sweep_alpha_N_S.txt
    figures/results_section/                      fig1_partial_sharing_calibration.png, fig2_recorded_vs_partial_sharing.png,
                                                  partial_sharing_numbers.txt
```

Decoding uses `Simulations_v2/code/decoding_scripts/pt_transfer_decode_population.m` unchanged,
for the synthetic populations and the recorded stimulus subsets alike (both scripts add it to the
path): z-score per cue, MATLAB `procrustes` on trial-averaged means, one held-out trial per
stimulus, ECOC linear SVM, self-decoding by 10-fold CV, and the self-consistent
correspondence-shuffled control. Do not edit that file without rerunning both folders.
Quote numbers from `results/figures/results_section/partial_sharing_numbers.txt`, not from this README.

## Run (`sweep_alpha_N_S`, 2026-10-07 23:45 → 2026-10-08 07:58, 494 min, complete: 275 computed + 55 imported points)

| | |
|---|---|
| α | 0, 0.2, 0.4, 0.6, 0.8 computed; 1 imported from `Simulations_v2/.../sweep_N_S_model3_D5_results.mat` |
| N | 5, 10, 20, …, 100 |
| n_stim | 10, 20, 30, 40, 50 (a 5-D latent needs more than 5 stimuli) |
| populations | 16 per point, trial_repeat 10, n_folds 10 (as v2) |
| seeds | numbered as in v2 and independent of α: every population has the same neurons as the Model 1 (`sweep_N_S`) and fully shared Model 3 populations at that grid point, and the α values are paired |
| order | 50 classes first, then 40, …, 10; α innermost, so any finished grid point has every α |
| cost | ≈ 12 / 45 / 95 / 170 / 268 s per point for 10 / 20 / 30 / 40 / 50 classes, nearly independent of N → ≈ 108 min per α, ≈ 9 h in total |

Launched from `code/decoding_scripts` with `matlab -batch procrustes_decoding_partial_sharing_sweep`
(default configuration). The file is rewritten after every point, so an interrupted run resumes
when the same command is rerun.

Checks done before launch (2026-10-07):
- the generator with K_s = K returns exactly the data of `generate_model3_trial_data.m`
  (n_latent 5, cap 'exact') for every seed and size tested;
- across α (N = 100, 50 classes, 40 seeds) the signal PR (≈ 4.0), median tuning depth (0.138)
  and fraction of depth-capped neurons (0.08) do not change, and the singular values of the
  cross-cue latent covariance are K_s ones plus chance-level values (0.52, 0.38, … at α = 0);
- α = 1 decoded inside the driver equals the imported v2 populations exactly (max difference 0).
  The decoder is only reproducible single-threaded (pool workers are): on the multithreaded client
  the Procrustes fit can land on a different but equally good rotation when N ≥ n_stim, which
  changes single-population accuracies by a few percent.

## Results file

`results/decoding_outputs/sweep_alpha_N_S_results.mat`, one variable `sweep`:

| Field | Content |
|---|---|
| `info`, `created`, `updated`, `cfg` | how to read the file; timestamps; run configuration |
| `N_list`, `S_list`, `alpha_list`, `n_shared_list`, `n_latent` | grid axes: dims 1, 2, 3 of every grid array |
| `metric_names`, `pair_labels` | the 9 accuracy columns (as v2); ordered cue pairs |
| `done`, `source` | [nN×nS×nA]; source 1 = computed here, 2 = imported (α = 1) |
| `mean`, `sem` | [nN×nS×nA×9] over populations, of pair-averaged accuracies |
| `ratio_mean`, `ratio_sem`, `ratio_of_means` | PT_rotation / control |
| `headroom_of_means`, `headroom_mean`, `headroom_sem` | (PT − control) / (self − control): the share of the gap between the shuffled floor and self-decoding that PT recovers; 0 = no shared structure, 1 = full transfer |
| `PR_true`, `PR_emp`, `count`, `frac_cv_capped`, `pdist_mean`, `runtime_s` | as v2, per grid point |
| `per_pop` | every population: `seeds` [nN×nS×R]; `pair_metrics`, `pop_metrics`, `ratio`, `headroom`, `mean_count`, `frac_cv_capped`, `pair_pdist`, `signal_PR`, `PR_emp`, `latent_cross_sv` [nN×nS×nA×R×…] |
| `summary` | table, one row per α × n_stim × N |

```matlab
S = load('results/decoding_outputs/sweep_alpha_N_S_results.mat');  sweep = S.sweep;
iN = find(sweep.N_list == 100);  iS = find(sweep.S_list == 50);
squeeze(sweep.mean(iN, iS, :, 6)).'            % PT_rotation at alpha = 0 ... 1
squeeze(sweep.headroom_of_means(iN, iS, :)).'  % share of the floor-to-self gap recovered
```

## Real-data stimulus subsampling (`real_stimsub`, 2026-10-08 11:25 → 14:06, 160 min, complete: 40 cells)

`code/decoding_scripts/real_data_stimulus_subsampling.m` decodes the four recorded populations
(FR V1 48, FR V2 112, KO V1 109, KO V2 146 neurons; [330 630] ms counts; ac, ec, ex; all 6 ordered
pairs) on random subsets of C ∈ {10, 20, 30, 40, 50} stimuli, with the same decoder and control as
the synthetic populations. Each subset is z-scored on its own trials, like a synthetic population of C
stimuli (the recordings' own incremental_pr script z-scores over all 50 stimuli and refits the SVM per
holdout; its 50-class accuracies agree with this script's to within 0.005, headroom to within 0.01).

| | |
|---|---|
| modes | `full_N`: all neurons; `matched_ratio`: N = round(N_full·C/50) random neurons per subset (N/C held at the recorded ratio) |
| repeats | 16 random stimulus subsets per C < 50; at C = 50 (all stimuli, all neurons, identical in both modes) 8 repeats of the holdouts and CV folds, computed once |
| silent neurons | a neuron with no spikes in a cue's selected trials stays in (MATLAB `zscore` maps it to zeros, as in the recordings' pipeline); FR V2 has 7 over all 50 stimuli; counted as `n_silent` |
| cost | ≈ 7–8 min per population at C = 50, ≈ 9 min per cell at C = 40; 160 min in total |
| output | `results/decoding_outputs/real_stimsub_results.mat` (variable `real`: [population × mode × C × …] arrays, per-subset metrics, stimulus and neuron indices, seeds, summary table), resumable |

Prediction from the family: headroom stays roughly flat in C (for C ≥ 20), while PT/control rises
with C at fixed α.

## How the results are read

The simulations calibrate what PT can show; they do not estimate "how much is shared" in the
brain. α is defined within this family (linear readout, K = 5, independent per-cue readouts), and
changing any of those would move the curves. The defensible statement is of the form: within this
family, at matched self-decoding and PR, the recorded populations are consistent with α in [a, b]
and inconsistent with α = 0.

First results (numbers in `results/figures/results_section/partial_sharing_numbers.txt`):
- Within-cue statistics do not depend on α (self-decoding, PR and depth capping flat across α at
  every grid point); α = 0 gives PT/control ≈ 1 and headroom ≈ 0.
- Headroom depends on α and hardly on the number of classes (C ≥ 20) or neurons (N ≥ 20); it is
  convex in α (≈ 0.09, 0.22, 0.38, 0.61, 0.93 for α = 0.2 … 1). PT/control at fixed α rises
  steeply with the number of classes, so headroom is the axis on which recordings and models are compared.
- At 50 classes the recorded headroom (0.15–0.25, same decoder) corresponds to α ≈ 0.30–0.45 within this family.
- On random subsets of 20–50 stimuli the recorded headroom stays flat (all neurons: FR V1 0.14–0.15,
  FR V2 0.15–0.18, KO V1 0.25–0.28, KO V2 0.21–0.23; same with N/C held), as the family predicts, while
  PT/control falls with fewer classes (KO V1 1.47 → 1.18). 10 classes is too noisy to read.
- Open mismatch: the recordings' floor (control / self, 0.33–0.45 at full size) is higher than the
  family's (0.20–0.24) at matched self-decoding and PR, and rises with N like Model 1's. Read through
  PT/control, which does not correct for the floor, the recordings would imply only α ≈ 0.18–0.29.

## Not run yet

- Shared-readout variant (`rho = 1`: each neuron keeps the same filter across cues, from the
  2026-09-25 meeting), three cues with a different α per pair (Model 5), PC-constrained runs
  (`cfg.pc_dims`), and latents with a decaying spectrum instead of K equal dimensions (a candidate
  explanation of the floor mismatch).
