function R = load_all_results(S)
%{
Filename: load_all_results.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Loads every neuron_decoding_rev result file into one struct:
      R.cue.<pop>, R.halves.<pop>, R.pc_cue.<pop>, R.pc_halves.<pop>   (pop = 'KO_V1', ...)
      R.across.<monkey>, R.pc_across.<monkey>                          (monkey = 'KO', 'FR')
    Each entry is the `res` struct written by the driver (see the driver headers),
    plus res.ctrl_m and res.ctrl_id_m: the control averaged over its permutations.
%}

if nargin < 1, S = figure_settings(); end
for p = 1:numel(S.pops)
    R.cue.(S.pops{p})       = load_one(S, fullfile('cue_transfer_incremental', [S.pops{p} '.mat']));
    R.halves.(S.pops{p})    = load_one(S, fullfile('halves_incremental',       [S.pops{p} '.mat']));
    R.pc_cue.(S.pops{p})    = load_one(S, fullfile('pc_sweep_cue',             [S.pops{p} '.mat']));
    R.pc_halves.(S.pops{p}) = load_one(S, fullfile('pc_sweep_halves',          [S.pops{p} '.mat']));
end
for m = 1:numel(S.monkeys)
    R.across.(S.monkeys{m})    = load_one(S, fullfile('across_area',     [S.monkeys{m} '.mat']));
    R.pc_across.(S.monkeys{m}) = load_one(S, fullfile('pc_sweep_across', [S.monkeys{m} '.mat']));
end
end

function res = load_one(S, rel)
f = fullfile(S.results_root, rel);
assert(exist(f, 'file') == 2, 'load_all_results: missing %s', f);
res = load(f, 'res').res;
assert(all(res.done(:)), 'load_all_results: %s is not finished', f);
res.ctrl_m    = mean(res.ctrl,    ndims(res.ctrl));       % mean over permutations
res.ctrl_id_m = mean(res.ctrl_id, ndims(res.ctrl_id));
end
