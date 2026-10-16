function R = load_ann_results(S)
%{
Filename: load_ann_results.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-13

Description:
    Loads every ANN_decoding_rev result file into one struct:
      R.cue.<key>, R.pc_cue.<key>     key = 'alexnet_pool1', ..., 'vgg_pool5'
      R.hier, R.pc_hier               AlexNet pool2 <-> pool5
    Each entry is the `res` struct written by the driver (see the driver
    headers), plus res.ctrl_m and res.ctrl_id_m: the control averaged over its
    permutations.
%}

if nargin < 1, S = ann_figure_settings(); end
for i = 1:numel(S.keys)
    R.cue.(S.keys{i})    = load_one(S, fullfile('ann_cue_transfer_incremental', [S.keys{i} '.mat']));
    R.pc_cue.(S.keys{i}) = load_one(S, fullfile('ann_pc_sweep_cue',             [S.keys{i} '.mat']));
end
R.hier    = load_one(S, fullfile('ann_hierarchy',          'alexnet_pool2_pool5.mat'));
R.pc_hier = load_one(S, fullfile('ann_pc_sweep_hierarchy', 'alexnet_pool2_pool5.mat'));
end

function res = load_one(S, rel)
f = fullfile(S.results_root, rel);
assert(exist(f, 'file') == 2, 'load_ann_results: missing %s', f);
res = load(f, 'res').res;
assert(all(res.done(:)), 'load_ann_results: %s is not finished', f);
res.ctrl_m    = mean(res.ctrl,    ndims(res.ctrl));       % mean over permutations
res.ctrl_id_m = mean(res.ctrl_id, ndims(res.ctrl_id));
end
