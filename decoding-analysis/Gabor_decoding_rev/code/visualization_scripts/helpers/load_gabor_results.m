function R = load_gabor_results(S)
%{
Filename: load_gabor_results.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-21

Description:
    Loads every Gabor_decoding_rev result file into one struct:
      R.cue.<model>   cue transfer vs number of cells (model = 'simple', 'complex')
      R.pc.<model>    PC sweep of cue transfer, whole model
    Each entry is the `res` struct written by the driver (see the driver
    headers), plus res.ctrl_m and res.ctrl_id_m: the control averaged over its
    permutations.
%}

if nargin < 1, S = gabor_figure_settings(); end
for i = 1:numel(S.models)
    R.cue.(S.models{i}) = load_one(S, fullfile('gabor_cue_transfer_incremental', [S.models{i} '.mat']));
    R.pc.(S.models{i})  = load_one(S, fullfile('gabor_pc_sweep',                 [S.models{i} '.mat']));
end
end

function res = load_one(S, rel)
f = fullfile(S.results_root, rel);
assert(exist(f, 'file') == 2, 'load_gabor_results: missing %s', f);
res = load(f, 'res').res;
assert(all(res.done(:)), 'load_gabor_results: %s is not finished', f);
res.ctrl_m    = mean(res.ctrl,    ndims(res.ctrl));       % mean over permutations
res.ctrl_id_m = mean(res.ctrl_id, ndims(res.ctrl_id));
end
