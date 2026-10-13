function out = bootstrap_metrics(Sa, Pa, Ca, unit_dim, S, seed)
%{
Filename: bootstrap_metrics.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Point estimates and percentile-bootstrap 95% confidence intervals of the
    normalized transfer metrics (transfer_metrics.m: ratios of POOLED means):
        ratio = PT / control,  diff = PT - control,
        headroom = (PT - control) / (self - control).
    The three arrays are first reshaped so that the resampling unit (neuron
    draw, or fold when a whole population is used once) is the row
    (to_resampling_units.m); self, PT and control of a unit are resampled together.

Inputs:
    Sa, Pa, Ca  self, PT and control accuracy arrays of the same shape (control
                already averaged over permutations); NaN entries ignored
    unit_dim    dimension of the arrays that holds the resampling unit
    S           figure_settings() (S.n_boot)
    seed        bootstrap seed

Output: struct with self, pt, ctrl, ratio, diff, headroom, *_ci [lo hi] and n_units
%}

Su = to_resampling_units(Sa, unit_dim);
Pu = to_resampling_units(Pa, unit_dim);
Cu = to_resampling_units(Ca, unit_dim);
m = transfer_metrics(Su, Pu, Cu);
out = struct('self', m.self, 'pt', m.pt, 'ctrl', m.ctrl, 'ratio', m.ratio, 'diff', m.diff, ...
    'headroom', m.headroom, 'n_units', size(Su, 1));
rs = RandStream('twister', 'Seed', seed);
nU = size(Su, 1);
B  = zeros(S.n_boot, 3);
for b = 1:S.n_boot
    i  = randi(rs, nU, nU, 1);
    mb = transfer_metrics(Su(i, :), Pu(i, :), Cu(i, :));
    B(b, :) = [mb.ratio, mb.diff, mb.headroom];
end
q = prctile(B, [2.5 97.5]);
out.ratio_ci    = q(:, 1)';
out.diff_ci     = q(:, 2)';
out.headroom_ci = q(:, 3)';
end
