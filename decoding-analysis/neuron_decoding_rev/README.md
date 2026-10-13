# neuron_decoding_rev — main-figure decoding with the corrected control

Revision reruns of the neural decoding behind the main figures of the J Neurosci
resubmission (Figs 3b, 4a, 4b, 4c, 12, 13 and their PC-dimension versions). Every
analysis scores Procrustes transfer (PT) against the **corrected control**: one
permutation of the stimuli is applied to both the source class means and the held-out
source trials, so the control has exactly the rotation's degrees of freedom. The old
scripts in `neuronal_decoding/` shuffled only the means and stayed at chance by
construction. The plan behind this folder is the Claude Doc "Main-Figure Rerun Plan".

Data are read in place from `../neuronal_decoding/neuronal_data` (not copied).

## Layout

```
neuron_decoding_rev/
  code/
    decoding_scripts/
      load_population.m              one monkey/area: [500 x N] spike counts per rendering + labels
      pt_transfer_core.m             THE decoder: one source -> target comparison (full space or PC-k)
      population_pr.m                participation ratio: raw (rebuttal definition) and noise-corrected
      transfer_metrics.m             pooled means, PT/control, PT - control, headroom
      multiclass_svmloader_PT.m      local copy of the trial loader
      run_cue_transfer_incremental.m Fig 3b (one point) + Fig 4c: cue transfer vs N, full space
      run_halves_incremental.m       Fig 4a (one point) + Fig 13: within-area halves vs half size
      run_across_area.m              Fig 4b + Fig 12: V1 <-> V2, both directions, 3 x 3 renderings
      run_pc_sweep.m                 k = 1..49 sweep for cue transfer, halves and across-area
      run_all.m                      MASTER: runs the drivers back to back, unattended, with a log
      run_driver.m                   runs one driver in its own workspace (used by run_all)
      progress_queue.m               progress lines inside the parfor loops
    visualization_scripts/           one script per figure series; each writes results/figures/<NN_series>/
      plot_sanity_old_vs_new.m           01  rerun vs published pipeline, bar by bar
      plot_paper_style_bars.m            02  Figs 3b, 4a, 4b, 12 in the published style + corrected control (full, k = 5)
      plot_paper_style_curves.m          03  Figs 4c, 13 in the published style + corrected control
      plot_normalized_metric_bars.m      04  PT/control, headroom, PT - control: full space vs k = 5, per pair, matched size
      plot_normalized_metrics_vs_size.m  05  the same metrics vs number of neurons / half size
      plot_pc_sweep_accuracy.m           06  self, PT, control, non-PT vs number of PCs
      plot_pc_sweep_metrics.m            07  PT/control, headroom, PT - control vs number of PCs
      plot_dimensionality.m              08  participation ratio and variance spectra
      helpers/                           figure_settings, load_all_results, plot_grouped_bars, bootstrap_metrics,
                                         to_resampling_units, mean_sem, save_figure, save_table
  results/
    decoding_outputs/<analysis>/<monkey>_<area>.mat   (across-area: <monkey>.mat)
    logs/run_all_<set>_<date>.log
    figures/<NN_series>/             figures (.png, .fig) + a .csv of every plotted number;
                                     index in results/figures/README.md
```

## How every comparison is decoded (`pt_transfer_core.m`)

1. Z-score each side per neuron over its 500 trials (as published).
2. **10-fold trial partition**: within each stimulus the 10 repeats are shuffled once and
   repeat k goes to fold k. Folds move whole rows, and source and target share the
   partition, so within-session noise correlations are kept.
3. For fold j, using trials outside fold j only: target class means, (PC mode) target
   PCA basis and the ECOC linear SVM (fitcecoc defaults); source class means and (PC
   mode) source basis; Procrustes fit from source means onto target means.
4. The fold-j model is tested on the **target's** fold-j trials (self-decoding, so the
   10 fold models are a 10-fold CV) and on the **source's** fold-j trials:
   no transform (non-PT), rotation only (`pt_rot`, the headline), full `b*X*T + c`
   (`pt_full`), and the corrected control (`ctrl`, n_perm permutations per fold;
   `ctrl_id` scores the same predictions against the permuted identities, ~1/50).
5. PC mode projects each side onto the top k directions of its own centered class
   means (as `procrustes_decoding_basic_incremental_pc.m`), so the rotation has
   k(k-1)/2 parameters. Non-PT in PC mode projects the source trials onto the target
   basis (same neurons only).

Summary metrics (`transfer_metrics.m`) are ratios of pooled means:
PT/control, PT − control, and headroom = (PT − control)/(self − control).

## Running

**Unattended (recommended):** the master script runs the drivers one after another and
logs everything printed to `results/logs/`:

```matlab
run_all('test')      % ~8 min, tiny configs in a temp folder: run once per machine first
run_all('all')       % every driver, ~20 h on one 8-worker machine
run_all('desktop')   % or split across two machines: ~11 h here ...
run_all('laptop')    % ... and ~9 h here
```

A job that fails is logged and skipped, the next one starts, and a summary of all jobs
is printed at the end. Calling `run_all` again with the same set resumes every driver
where it stopped. While running, each grid point prints a header, a progress line about
every 10% of its calls (with minutes left for that point), and a result line with the
hours left for that population.

**One driver at a time:** run any driver from MATLAB (the current folder does not
matter). **Each driver runs all populations by default.** To change a setting, set an
override first:

```matlab
cfg_override = struct('populations', {{'KO','V2'; 'FR','V1'}});   % or 'monkeys', {{'KO'}} for across-area
run_cue_transfer_incremental
clear cfg_override    % before running another driver
```

Any field of the driver's `cfg` block can be overridden; unknown fields are refused.
Each population writes its own file, saved after every grid point (N, half size,
draw or direction). Re-running the same configuration resumes where it stopped; a file
made with a different configuration stops the run (set `cfg.overwrite = true` or a new
`cfg.tag`). Seeds are fixed per (population, grid point, draw, item), so results do not
depend on which machine runs which population. Commit each machine's result files to
git; the two machines never write the same file.

Defaults: 20 neuron draws (1 when the whole population is used), 10 folds, 10 control
permutations per fold, PR from 20 half-splits, window [330 630] (30–330 ms after onset).

## Run time (measured 2026-10-11, R2025b, 8 process workers)

One `pt_transfer_core` call (10 SVM fits + predictions) takes 38–48 s on a worker,
almost independent of the number of neurons or PCs: the cost is the 1225 one-vs-one
learners of each fit. With 8 workers that is about 650 calls per hour. Call counts of
the default grids and the resulting time on one 8-worker machine:

| Driver | Calls | Time, one machine |
| --- | --- | --- |
| run_cue_transfer_incremental | 5,424 (FR V1 606, FR V2 1,566, KO V1 1,326, KO V2 1,926) | ~8.3 h |
| run_halves_incremental | 1,740 (240, 480, 420, 600) | ~2.7 h |
| run_across_area | 720 (360 per monkey) | ~1.1 h |
| run_pc_sweep, cue | 432 (108 per population) | ~0.7 h |
| run_pc_sweep, halves | 4,020 (FR V1 780, others 1,080 each) | ~6.2 h |
| run_pc_sweep, across | 648 (324 per monkey) | ~1.0 h |
| **Total** | **12,984** | **~20 h** |

The two-machine split in `run_all` gives the faster machine ('desktop') about 11 h of this
machine's speed (cue transfer vs N, across areas, PC sweep of cue transfer and across
areas) and the other ('laptop') about 9 h (halves vs size, PC sweep of halves). Each
grid point prints its run time, so the real speed ratio is visible after the first few
points.

Results are reproducible on process workers (single-threaded). Without a pool, parfor
runs on the client with multithreaded BLAS and the last digits can differ.
