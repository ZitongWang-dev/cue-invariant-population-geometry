# Noise-correlation injection (rebuttal Point 8): rationale and methods

Zitong Wang · 2026-10-24 · code: `neuron_decoding_rev/code/decoding_scripts/run_nc_injection.m`,
figures: `neuron_decoding_rev/code/visualization_scripts/plot_nc_injection.m`

## 1. What the reviewer asks

Reviewer 2, major point 1: the pseudo-population removes the noise correlations (NC) that shape
"local manifold curvature, effective dimensionality, and the orientation of trial-to-trial
variability". The reviewer proposes "a simulation injecting structured noise correlations of
realistic magnitude (mean pairwise r = 0.05-0.15) into the pseudo-population", to test whether PT
decoding is robust or biased by the independence assumption.

**We inject into the recorded pseudo-population, not into the simulated models.** "Simulation" in
the review means the simulated perturbation of the real data. Three places in our materials still
say otherwise and must change after the run:

- Prof. Lee's third yellow paragraph of Point 8 ("injected ... into the synthetic populations ...
  organized according to the sizes of the simultaneously recorded groups");
- his cyan note "Add noise to Model 1 / Model 3";
- the meeting notes ("Injected noise of different levels into simulations").

Why the reviewer's range fits our data: same-session pairs already carry NC of 0.11–0.20 (mean
0.135), but most of the pseudo-population's pairs (81% in KO V1) come from different sessions and
carry 0.01–0.04, so the population mean is only 0.03–0.05. A mean of 0.05–0.15 over all pairs is what the
pseudo-population would look like if every pair behaved like a same-session pair.

## 2. Design

| Choice | Decision | Why |
| --- | --- | --- |
| Populations | FR V1 (48), FR V2 (112), KO V1 (109), KO V2 (146), whole | Each population is compared with itself before and after injection, so no V1/V2 size matching is needed. Full size has the most pairs, so it is the strictest test. No neuron-draw variability. |
| Conditions | original; affine (NC removed); injected mean NC 0, 0.05, 0.10, 0.15 | Original = the real data; affine = the clean start (Point 7's decorrelation); four injected levels span the reviewer's range. Injected 0 differs from affine: its pairwise NC still spreads around 0. |
| Space | Full neuron space | The reviewer's concern is single-trial decoding in the space actually used. |
| Pairs | All 6 ordered cue pairs, pooled | As the main analyses. |
| Seeds | 10 | Each seed is one injection draw with its own fold partition. |

## 3. Method

**Step 1, clean start.** Remove the existing NC with the deterministic affine trial shift of Point 7
(`noise_correlation_analysis_affine.m`): neuron j of a session gets permutation j of
t → (a·t + b) mod 10, so neurons recorded together lose their trial alignment. Each neuron keeps its
own counts per stimulus.

**Step 2, build a correlation matrix.** The goal is a matrix whose off-diagonal entries (pairwise
NC) are roughly Gaussian with a chosen mean and SD.

- Drawing every entry independently does not work. (a) The entries of a correlation matrix must be
  mutually consistent (positive semidefinite: every weighted sum of neurons needs a variance ≥ 0).
  Independent draws stay valid only if SD < (1 − mean)/(2√N), i.e. below 0.035 for KO V2 and
  0.06 for FR V1, narrower than the measured spread (0.093). (b) Independent entries have no
  shared modes, so the noise stays nearly round and barely changes the "orientation" or
  "dimensionality" the reviewer cites.
- Instead, build R from a few shared fluctuations (a factor model):
  R = L·Lᵀ + D, with L = [a·1, b·B], B (N × k) drawn from N(0, 1), D = diag(1 − a² − b²·|Bᵢ|²).
  In words: on every trial, each neuron's deviation = a × one common fluctuation + its random
  sensitivities (b·B) to k other shared fluctuations + private noise (D). The NC of a pair is then
  rᵢⱼ = a² + b²·Σₘ BᵢₘBⱼₘ: roughly Gaussian across pairs with **mean a²** and **SD ≈ b²√k**.
  R is valid by construction, and the k + 1 shared modes give the noise real structure.
- k = 5. A neuron whose shared variance a² + b²|Bᵢ|² would exceed 0.95 at the top level has its
  loadings shrunk (rare: 0–1 neurons per seed in the test).
- **One R per population and seed, for all 50 stimuli and all 3 renderings.** A pair's NC is
  treated as a property of the two neurons. Supported by the data: per pair, the variance of NC
  across stimuli is 0.108–0.110, about 1/9, which is what 10-trial sampling noise alone gives.
- **Paired levels:** B and b are the same at every level; only a changes.

**Step 3, impose R on the real counts (Iman–Conover).** For every rendering and stimulus, draw 10
"trials" from N(0, R) and reorder each neuron's 10 real counts so their ranks follow its column.
Every neuron keeps exactly its own counts per stimulus, so means, variances, z-scoring and the
trial-averaged geometry are unchanged; only how trials line up across neurons changes. The result
does not depend on the starting trial order, so no original NC carries over.

**Step 4, set the spread from the data.** Target SD = SD across same-session pairs of each pair's
mean NC (over 50 stimuli × 3 renderings), minus the estimation noise of that mean, pooled over the
four populations. Measured: **0.093** (3,768 pairs; raw SD 0.097, estimation noise 0.028; mean
0.135). Fixed across levels.

**Step 5, calibrate.** Ten trials and integer counts weaken the achieved correlation (about 0.8×), so
a² and b² are rescaled over 3 rounds on 2 separate draws until the achieved mean (top level) and SD
match the targets.

**Step 6, decode.** `pt_transfer_core.m` unchanged: z-score per side, 10-fold row partition shared by
source and target, ECOC linear SVM, rotation-only PT, corrected control (10 permutations per fold).
The partition and permutations depend on (population, seed, pair) only, so all six conditions of a
seed are decoded on the same folds.

## 4. Measurements

- **NC metric (as Point 7):** Pearson correlation of two neurons' counts over the 10 repeats of a
  stimulus; per pair, the mean over 50 stimuli × 3 renderings (silent neurons excluded). Reported for
  all pairs, same-session and cross-session pairs, plus the spread (estimation noise removed).
- **Noise participation ratio:** PR of the correlation matrix of the residuals (trial minus its
  stimulus mean), mean over renderings: the dimensionality of the trial-to-trial noise.
- **Decoding:** self-decoding, PT, control, PT/control, headroom = (PT − control)/(self − control),
  pooled means over pairs and folds, 95% CIs by bootstrap over seeds.
- **Readout fixed in advance:** Δheadroom of each injected level versus the affine condition,
  paired over seeds, with its 95% CI. Scale for "small": in the partial-sharing calibration
  (Point 5), headroom rises by ~0.15 per 0.2 of shared fraction α near the recorded range.
- The original-vs-affine comparison also redoes the Point 7-II pseudo-population result with the
  corrected control (the version in the summary still uses the old control).

## 5. Running

```matlab
% quick check of a machine (~2 min, writes results/decoding_outputs/nc_injection_test/)
cfg_override = struct('populations', {{'FR','V1'}}, 'n_seeds', 1, 'n_perm', 2, ...
    'pairs', {{'ac','ec'}}, 'tag', 'nc_injection_test');
run_nc_injection
clear cfg_override

% full run: 4 populations x 6 conditions x 6 pairs x 10 seeds = 1,440 core calls, ~2.3 h on 8 workers
run_nc_injection

% figures
plot_nc_injection            % or: nc_tag = 'nc_injection_test'; plot_nc_injection
```

- Output: `neuron_decoding_rev/results/decoding_outputs/nc_injection/<monkey>_<area>.mat` (one
  `res` per population; field list in the driver header). Saved after every seed; re-running
  resumes. To split across machines: `cfg_override = struct('populations', {{'KO','V2'; 'FR','V1'}})`.
  The target SD is measured on all four populations on every machine, so it is identical everywhere.
- Figures (`results/figures/09_nc_injection/`): `nc_injection_metrics` (per condition),
  `nc_injection_change` (paired change from affine vs achieved NC), `nc_injection_check`
  (achieved NC, spread, noise PR, per-pair NC distributions), `nc_injection.csv`.

## 6. Quick test (2026-10-24, FR V1, one pair ac→ec, 2 seeds, 2 permutations)

- Calibration converged in 3 rounds (round 1 achieved 0.119 for a target of 0.150).
- Achieved mean NC 0.000 / 0.051 / 0.101 / 0.150, same- and cross-session alike (0.142 / 0.152
  at the top level); spread 0.090–0.104.
- Noise PR: affine 39.7, original 26.2, injected 34 → 20 from level 0 to 0.15. The injection really
  changes the noise dimensionality; level 0.10 matches the original data's.
- The affine data keep a residual spread of 0.036 (mean −0.003): the deterministic shift preserves
  some trial structure (e.g. slow drifts). The injected conditions do not inherit it.
- Decoding numbers from this test are too noisy to read (one pair, two seeds).

## 7. Open choices (not in the first run)

- A tuning-dependent variant: higher NC for more similarly tuned pairs (fit r_noise ≈ c + d·r_signal
  on same-session pairs).
- k = 5 PC space as well (doubles the cost).
- V2 subsampled to the V1 count, one new neuron draw per seed (+~50% cost), for numbers matching Fig 3b.
- Model 1 with injected NC, to show the control still tracks PT under correlated noise (probably
  unnecessary: in the recorded data the control sees the same noise as PT).

## 8. Text to update after the run

- Point 8, Prof. Lee's third paragraph: "synthetic populations" → the recorded pseudo-population;
  describe the factor-model structure and the data-matched spread; fill in [RESULT].
- Point 8, the cyan to-do and the "overall story" line (currently "add NC to synthetic data").
- Point 7-II: replace the old-control numbers with the original-vs-affine result from this run.
