%{
% FILENAME: figure_partial_sharing_results.m   (Simulations_v3_partial_sharing)
% AUTHOR:   Zitong Wang
% DATE:     2026-10-08
%
% DESCRIPTION:
%   Paper-style figures, and every number quoted, for the partially shared
%   family (alpha = K_s / K, K = 5) and its comparison with the recorded
%   populations. Counterpart of Simulations_v2/code/visualization_scripts/
%   figure_results_section.m.
%
%   headroom = (PT - control) / (self - control), from mean accuracies: the
%   share of the gap between the correspondence-shuffled floor and
%   self-decoding that PT recovers (0 = nothing shared, 1 = full transfer).
%   Ratios are ratios of mean accuracies. "Mean over N >= 20" averages a
%   grid quantity over N = 20..100.
%
%   FIGURE 1  fig1_partial_sharing_calibration.png   (synthetic)
%     A  headroom vs alpha, one line per number of classes
%     B  PT / control vs number of classes, one line per alpha
%     C  headroom vs N at 50 classes, one line per alpha
%   FIGURE 2  fig2_recorded_vs_partial_sharing.png
%     A  headroom vs N at 50 classes: alpha family and the recorded populations
%        (neuron subsampling, neuronal_decoding incremental_pr results)
%     B  recorded headroom vs number of classes (random stimulus subsets, all
%        neurons; real_data_stimulus_subsampling.m) over each alpha's model
%        headroom (mean over N >= 20)
%     C  floor (control / self) vs neurons per class: alpha family (mean over
%        alpha), Model 1 (Simulations_v2) and the recorded populations
%   TEXT      partial_sharing_numbers.txt
%
% INPUTS:
%   results/decoding_outputs/sweep_alpha_N_S_results.mat
%   results/decoding_outputs/real_stimsub_results.mat           (optional; partial runs fine)
%   ../Simulations_v2/results/decoding_outputs/sweep_N_S_results.mat   (Model 1)
%   ../neuronal_decoding/results/decoding_outputs/procrustes_decoding_basic_incremental_pr_results/<M>/<A>/
%
% OUTPUT:
%   results/figures/results_section/
%}

%% Configuration
clc; clear;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
sim_root = fullfile(script_dir, '..', '..');
out_dir  = fullfile(sim_root, 'results', 'figures', 'results_section');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
export_formats = {'png'};          % add 'tif', 'pdf' or 'eps' for submission
font = 'Arial';  fs = 8;
min_N_avg = 20;                    % "mean over N >= 20"
match_N   = [40 100];              % model N range used for 50-class comparisons with the recordings

alpha_hex = {'#ef9b74', '#e87c4c', '#dc5d28', '#b64919', '#8e3812', '#62260a'};   % alpha 0 -> 1 (validated ordinal ramp)
nstim_hex = {'#b9a7cf', '#9c84bd', '#7f62a9', '#634490', '#4a2f70'};              % n_stim 10 -> 50
c_m1  = [42 120 214] / 255;        % Model 1 (as Simulations_v2)
c_fam = [220 93 40] / 255;         % partially shared family as a whole (alpha ramp mid step)
c_ref = [0.55 0.55 0.53];
ink2  = [0.32 0.32 0.30];
real_pops = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};
real_col  = [0.10 0.10 0.10; 0.58 0.58 0.58; 0.10 0.10 0.10; 0.58 0.58 0.58];   % V1 dark, V2 light
real_ls   = {'-', '-', '--', '--'};                                               % FR solid, KO dashed
real_mk   = {'d', 'd', 's', 's'};
c = struct('self', 1, 'notf', 2, 'ctrl', 4, 'pt', 6, 'ctrl_id', 9);

%% Load
L = load(fullfile(sim_root, 'results', 'decoding_outputs', 'sweep_alpha_N_S_results.mat'), 'sweep');  sw = L.sweep;
A = sw.alpha_list;  N = sw.N_list;  S = sw.S_list;  nA = numel(A);  nS = numel(S);
kN = N >= min_N_avg;  kM = N >= match_N(1) & N <= match_N(2);
i50 = find(S == 50);
self = sw.mean(:,:,:,c.self);  pt = sw.mean(:,:,:,c.pt);  ctrl = sw.mean(:,:,:,c.ctrl);
hr = sw.headroom_of_means;  ratio = sw.ratio_of_means;  PR = sw.PR_emp;
hr_avg    = reshape(mean(hr(kN,:,:), 1, 'omitnan'), nS, nA);       % [nS x nA]
ratio_avg = reshape(mean(ratio(kN,:,:), 1, 'omitnan'), nS, nA);
cs_fam    = mean(ctrl, 3) ./ mean(self, 3);                          % floor of the family [nN x nS] (alpha-free)
A_col = hexes(alpha_hex);  S_col = hexes(nstim_hex);
if size(A_col, 1) ~= nA, A_col = interp1(linspace(0,1,size(A_col,1)), A_col, linspace(0,1,nA)); end

L1 = load(fullfile(sim_root, '..', 'Simulations_v2', 'results', 'decoding_outputs', 'sweep_N_S_results.mat'), 'sweep');
m1 = struct('N', L1.sweep.N_list, 'S', L1.sweep.S_list, 'self', L1.sweep.model1.mean(:,:,c.self), ...
            'ctrl', L1.sweep.model1.mean(:,:,c.ctrl), 'pt', L1.sweep.model1.mean(:,:,c.pt));
m1.cs = m1.ctrl ./ m1.self;

real_base = fullfile(sim_root, '..', 'neuronal_decoding', 'results', 'decoding_outputs', ...
                     'procrustes_decoding_basic_incremental_pr_results');
rec = [];
for p = 1:size(real_pops, 1)
    rec = [rec, load_real(fullfile(real_base, real_pops{p,1}, real_pops{p,2}), [real_pops{p,1} ' ' real_pops{p,2}], c)]; %#ok<AGROW>
end

RS = [];
f_rs = fullfile(sim_root, 'results', 'decoding_outputs', 'real_stimsub_results.mat');
if exist(f_rs, 'file'), Lr = load(f_rs, 'real'); RS = Lr.real; end

%% Numbers quoted in the write-up
fid = fopen(fullfile(out_dir, 'partial_sharing_numbers.txt'), 'w');
say(fid, 'Numbers for the partial-sharing write-up. Synthetic: sweep_alpha_N_S (K = %d latent dimensions per cue,\n', sw.n_latent);
say(fid, 'alpha = K_s/K, %d populations per point; alpha = 1 imported from Simulations_v2 sweep_N_S_model3_D5).\n', sw.cfg.R_pop);
say(fid, 'headroom = (PT - control)/(self - control) and ratio = PT/control, both from mean accuracies.\n\n');

say(fid, '--- design check: within-cue statistics across alpha (all %d grid points) ---\n', nnz(all(sw.done, 3)));
rng_self = max(self, [], 3) - min(self, [], 3);  rng_pr = max(PR, [], 3) - min(PR, [], 3);
say(fid, 'self-decoding: max range across alpha %.3f (median %.3f); mean SEM over populations %.3f\n', ...
    max(rng_self(:)), median(rng_self(:)), mean(sw.sem(:,:,:,c.self), 'all', 'omitnan'));
say(fid, 'participation ratio (z-scored means): max range across alpha %.2f (median %.2f)\n', max(rng_pr(:)), median(rng_pr(:)));
comp = reshape(sw.source(1,1,:) == 1, 1, []);              % alphas computed here (not imported)
fc = sw.frac_cv_capped(:,:,comp);
say(fid, 'fraction of depth-capped neurons (computed alphas): %.3f-%.3f\n', min(fc(:)), max(fc(:)));
xs = squeeze(mean(sw.per_pop.latent_cross_sv(:, :, 1, :, :), [1 2 4], 'omitnan')).';
say(fid, 'alpha = 0: cross-cue latent canonical correlations (mean) %s\n', mat2str(xs, 2));

say(fid, '\n--- alpha = 0 (no shared structure at the recorded dimensionality) ---\n');
r0 = ratio(kN,:,1);  h0 = hr(kN,:,1);
say(fid, 'N >= %d: PT/control median %.2f, range %.2f-%.2f; headroom %.2f to %.2f\n', min_N_avg, median(r0(:)), min(r0(:)), max(r0(:)), min(h0(:)), max(h0(:)));
cid = sw.mean(:,:,comp,c.ctrl_id) - 1 ./ S;
say(fid, 'control scored vs permuted identities (computed alphas): max |acc - 1/n_stim| = %.3f\n', max(abs(cid(:))));

say(fid, '\n--- headroom, mean over N >= %d (rows alpha, columns n_stim %s) ---\n', min_N_avg, mat2str(S));
for iA = 1:nA, say(fid, 'alpha %.1f: %s\n', A(iA), sprintf('%6.2f', hr_avg(:,iA))); end
say(fid, 'range over N >= %d at 50 classes, per alpha:\n', min_N_avg);
for iA = 1:nA, say(fid, '   alpha %.1f: %.2f-%.2f\n', A(iA), min(hr(kN,i50,iA)), max(hr(kN,i50,iA))); end
say(fid, '\n--- PT / control, mean over N >= %d (rows alpha, columns n_stim %s) ---\n', min_N_avg, mat2str(S));
for iA = 1:nA, say(fid, 'alpha %.1f: %s\n', A(iA), sprintf('%6.2f', ratio_avg(:,iA))); end

say(fid, '\n--- 50 classes, N = %d-%d (self and PR in the recorded range) ---\n', match_N);
for iA = 1:nA
    say(fid, 'alpha %.1f: self %.3f-%.3f | PT/control %.2f-%.2f | headroom %.2f-%.2f | control/self %.2f-%.2f | PR %.1f-%.1f\n', A(iA), ...
        min(self(kM,i50,iA)), max(self(kM,i50,iA)), min(ratio(kM,i50,iA)), max(ratio(kM,i50,iA)), min(hr(kM,i50,iA)), max(hr(kM,i50,iA)), ...
        min(ctrl(kM,i50,iA) ./ self(kM,i50,iA)), max(ctrl(kM,i50,iA) ./ self(kM,i50,iA)), min(PR(kM,i50,iA)), max(PR(kM,i50,iA)));
end
j50 = find(m1.S == 50);
say(fid, 'Model 1 (v2), 50 classes, N = %d-%d: control/self %.2f-%.2f\n', match_N, min(m1.cs(m1.N >= match_N(1), j50)), max(m1.cs(m1.N >= match_N(1), j50)));

% recorded at 50 classes: incremental_pr (full N) and, when available, the same-decoder C = 50 cell
h_curve = hr_avg(i50, :);                                    % headroom(alpha) at 50 classes, mean over N >= 20
r_curve = reshape(mean(ratio(kM, i50, :), 1), 1, nA);        % PT/control(alpha), 50 classes, N = 40-100
say(fid, '\n--- recorded populations, 50 classes, full population ---\n');
say(fid, 'alpha implied by headroom: interpolated on the 50-class model curve (mean over N >= %d)\n', min_N_avg);
say(fid, 'alpha implied by PT/control: interpolated on the 50-class model curve (mean over N = %d-%d)\n', match_N);
say(fid, '| Population | N | Self | PT | Control | PT/control | Headroom | Control/self | PR | alpha (headroom) | alpha (ratio) | source |\n');
say(fid, '| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |\n');
for p = 1:numel(rec)
    e = numel(rec(p).N);
    say(fid, '| %s | %d | %.3f | %.3f | %.3f | %.2f | %.2f | %.2f | %.1f | %.2f | %.2f | incremental_pr |\n', rec(p).name, rec(p).N(e), ...
        rec(p).self(e), rec(p).pt(e), rec(p).ctrl(e), rec(p).ratio(e), rec(p).headroom(e), rec(p).cs(e), rec(p).pr(e), ...
        inv_curve(h_curve, A, rec(p).headroom(e)), inv_curve(r_curve, A, rec(p).ratio(e)));
end
if ~isempty(RS)
    iC50 = find(RS.C_list == 50);
    for p = 1:numel(RS.populations)
        if ~RS.done(p, 1, iC50), continue; end
        a = squeeze(RS.mean(p, 1, iC50, :));
        say(fid, '| %s | %d | %.3f | %.3f | %.3f | %.2f | %.2f | %.2f | %.1f | %.2f | %.2f | stimsub C = 50 (same decoder as models) |\n', ...
            strrep(RS.populations{p}, '_', ' '), RS.N_used(p, 1, iC50), a(c.self), a(c.pt), a(c.ctrl), RS.ratio_of_means(p, 1, iC50), ...
            RS.headroom_of_means(p, 1, iC50), a(c.ctrl) / a(c.self), RS.PR_emp(p, 1, iC50), ...
            inv_curve(h_curve, A, RS.headroom_of_means(p, 1, iC50)), inv_curve(r_curve, A, RS.ratio_of_means(p, 1, iC50)));
    end
end
say(fid, 'recorded headroom over N (incremental_pr, 50 classes):\n');
for p = 1:numel(rec)
    k = rec(p).N >= min_N_avg;
    say(fid, '   %s: %.2f-%.2f (N = %d..%d); control/self %.2f-%.2f\n', rec(p).name, min(rec(p).headroom(k)), max(rec(p).headroom(k)), ...
        min(rec(p).N(k)), max(rec(p).N(k)), min(rec(p).cs(k)), max(rec(p).cs(k)));
end

if ~isempty(RS)
    say(fid, '\n--- recorded populations on random subsets of C stimuli (%d subsets per C < 50; %d repeats at C = 50) ---\n', ...
        RS.cfg.n_subsets, RS.cfg.n_repeat_full);
    for m = 1:numel(RS.modes)
        say(fid, 'mode %s: headroom (SEM over subsets) by C = %s\n', RS.modes{m}, mat2str(RS.C_list));
        for p = 1:numel(RS.populations)
            s_ = '';
            for iC = 1:numel(RS.C_list)
                if RS.done(p, m, iC)
                    s_ = [s_ sprintf('  %5.2f (%.2f)', RS.headroom_of_means(p, m, iC), RS.headroom_sem(p, m, iC))]; %#ok<AGROW>
                else
                    s_ = [s_ '      (n.d.)  ']; %#ok<AGROW>
                end
            end
            say(fid, '   %-6s %s\n', RS.populations{p}, s_);
        end
        say(fid, '   PT/control by C:\n');
        for p = 1:numel(RS.populations)
            say(fid, '   %-6s %s\n', RS.populations{p}, sprintf('  %5.2f', squeeze(RS.ratio_of_means(p, m, :))));
        end
        say(fid, '   control/self by C:\n');
        for p = 1:numel(RS.populations)
            say(fid, '   %-6s %s\n', RS.populations{p}, sprintf('  %5.2f', squeeze(RS.mean(p, m, :, c.ctrl) ./ RS.mean(p, m, :, c.self))));
        end
    end
    say(fid, 'model headroom by C (mean over N >= %d), for reference:\n', min_N_avg);
    for iA = 1:nA, say(fid, '   alpha %.1f: %s\n', A(iA), sprintf('  %5.2f', hr_avg(:, iA))); end
end
fclose(fid);

%% FIGURE 1: calibration of the partially shared family
fig1 = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 7.2 2.7]);
t1 = tiledlayout(fig1, 1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

% A: headroom vs alpha, one line per number of classes
ax = nexttile(t1);  hold(ax, 'on');
plot(ax, [0 1], [0 1], ':', 'Color', c_ref, 'LineWidth', 0.75);
yline(ax, 0, '-', 'Color', c_ref, 'LineWidth', 0.5);
hS = gobjects(1, nS);
for iS = 1:nS
    hS(iS) = plot(ax, A, hr_avg(iS,:), '-o', 'Color', S_col(iS,:), 'MarkerFaceColor', S_col(iS,:), 'MarkerSize', 3, 'LineWidth', 1.25);
end
text(ax, 0.80, 0.90, 'y = \alpha', 'Color', ink2, 'FontSize', fs - 1, 'FontName', font, 'Rotation', 33);
set(ax, 'XLim', [-0.03 1.03], 'YLim', [-0.08 1.02], 'XTick', 0:0.2:1);
style(ax, font, fs, 'Shared fraction \alpha', 'Headroom', 'Headroom vs. shared fraction', 'A');
legend(ax, hS, arrayfun(@(s) sprintf('%d classes', s), S, 'UniformOutput', false), 'Location', 'northwest', ...
    'Box', 'off', 'FontSize', fs - 1, 'FontName', font);

% B: PT / control vs number of classes, one line per alpha
ax = nexttile(t1);  hold(ax, 'on');
yline(ax, 1, '-', 'Color', c_ref, 'LineWidth', 0.75);
for iA = 1:nA
    plot(ax, S, ratio_avg(:,iA), '-o', 'Color', A_col(iA,:), 'MarkerFaceColor', A_col(iA,:), 'MarkerSize', 3, 'LineWidth', 1.25);
    text(ax, S(end) + 1.5, ratio_avg(end,iA), sprintf('%.1f', A(iA)), 'Color', ink2, 'FontSize', fs - 1, 'FontName', font);
end
set(ax, 'YScale', 'log', 'YLim', [0.85 5], 'YTick', [1 1.5 2 3 4], 'YTickLabel', {'1', '1.5', '2', '3', '4'}, ...
    'XLim', [5 58], 'XTick', S);
style(ax, font, fs, 'Number of stimulus classes', 'PT / control', 'Ratio depends on class count', 'B');

% C: headroom vs N at 50 classes, one line per alpha
ax = nexttile(t1);  hold(ax, 'on');
yline(ax, 0, '-', 'Color', c_ref, 'LineWidth', 0.5);
hA = gobjects(1, nA);
for iA = 1:nA
    hA(iA) = plot(ax, N, hr(:,i50,iA), '-o', 'Color', A_col(iA,:), 'MarkerFaceColor', A_col(iA,:), 'MarkerSize', 2.5, 'LineWidth', 1.25);
end
set(ax, 'XLim', [0 105], 'XTick', 0:25:100, 'YLim', [-0.08 1.02]);
style(ax, font, fs, 'Number of neurons N', 'Headroom', 'Headroom vs. N (50 classes)', 'C');
lg = legend(ax, hA, arrayfun(@(a) sprintf('\\alpha = %.1f', a), A, 'UniformOutput', false), 'Box', 'off', ...
    'FontSize', fs - 1, 'FontName', font);
lg.Layout.Tile = 'east';
save_fig(fig1, fullfile(out_dir, 'fig1_partial_sharing_calibration'), export_formats);

%% FIGURE 2: recorded populations vs the family
fig2 = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 7.2 2.7]);
t2 = tiledlayout(fig2, 1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
mk = {'MarkerSize', 3, 'LineWidth', 1.1};

% A: headroom vs N at 50 classes
ax = nexttile(t2);  hold(ax, 'on');
yline(ax, 0, '-', 'Color', c_ref, 'LineWidth', 0.5);
for iA = 1:nA
    h = plot(ax, N, hr(:,i50,iA), '-', 'Color', A_col(iA,:), 'LineWidth', 1.25);
    if iA == ceil(nA / 2), hFam = h; end
end
hR = gobjects(1, numel(rec));
for p = 1:numel(rec)
    hR(p) = plot(ax, rec(p).N, rec(p).headroom, [real_ls{p} real_mk{p}], 'Color', real_col(p,:), 'MarkerFaceColor', real_col(p,:), mk{:});
end
set(ax, 'XLim', [0 150], 'XTick', 0:50:150, 'YLim', [-0.08 1.02]);
style(ax, font, fs, 'Number of neurons N', 'Headroom', '50 classes', 'A');

% B: recorded headroom vs number of classes (stimulus subsampling, all neurons)
ax = nexttile(t2);  hold(ax, 'on');
yline(ax, 0, '-', 'Color', c_ref, 'LineWidth', 0.5);
for iA = 1:nA
    plot(ax, S, hr_avg(:,iA), '-', 'Color', A_col(iA,:), 'LineWidth', 1.25);
    text(ax, 52.5, hr_avg(end,iA), sprintf('%.1f', A(iA)), 'Color', ink2, 'FontSize', fs - 1, 'FontName', font);
end
if ~isempty(RS)
    jit = [-0.9 -0.3 0.3 0.9];
    for p = 1:numel(RS.populations)
        d = squeeze(RS.done(p, 1, :)).';
        if ~any(d), continue; end
        errorbar(ax, RS.C_list(d) + jit(p), as_row(RS.headroom_of_means(p, 1, d)), as_row(RS.headroom_sem(p, 1, d)), ...
            [real_ls{p} real_mk{p}], 'Color', real_col(p,:), 'MarkerFaceColor', real_col(p,:), 'CapSize', 0, mk{:});
    end
else
    text(ax, 0.5, 0.5, 'stimulus subsampling not run', 'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'Color', ink2, 'FontSize', fs, 'FontName', font);
end
set(ax, 'XLim', [5 57], 'XTick', S, 'YLim', [-0.3 1.02]);
style(ax, font, fs, 'Number of stimulus classes', 'Headroom', 'Stimulus subsets', 'B');

% C: floor vs neurons per class
ax = nexttile(t2);  hold(ax, 'on');
xline(ax, 1, ':', 'Color', c_ref, 'LineWidth', 0.75);
for iS = 1:numel(m1.S)
    h = plot(ax, m1.N / m1.S(iS), m1.cs(:,iS), '-', 'Color', [c_m1 0.55], 'LineWidth', 0.9);
    if iS == 1, hM1 = h; end
end
for iS = 1:nS
    h = plot(ax, N / S(iS), cs_fam(:,iS), '-', 'Color', [c_fam 0.7], 'LineWidth', 0.9);
    if iS == 1, hF = h; end
end
for p = 1:numel(rec)
    plot(ax, rec(p).N / 50, rec(p).cs, [real_ls{p} real_mk{p}], 'Color', real_col(p,:), 'MarkerFaceColor', real_col(p,:), mk{:});
end
if ~isempty(RS)
    for p = 1:numel(RS.populations)
        d = squeeze(RS.done(p, 1, :)).';
        if ~any(d), continue; end
        plot(ax, as_row(RS.N_used(p, 1, d)) ./ RS.C_list(d), as_row(RS.mean(p, 1, d, c.ctrl) ./ RS.mean(p, 1, d, c.self)), ...
            real_mk{p}, 'Color', real_col(p,:), 'MarkerFaceColor', 'w', 'MarkerSize', 3.5, 'LineWidth', 0.9);
    end
end
set(ax, 'XScale', 'log', 'XLim', [0.09 30], 'XTick', [0.1 1 10], 'XTickLabel', {'0.1', '1', '10'}, 'YLim', [0 1.05]);
style(ax, font, fs, 'Neurons per class (N / C)', 'Control / self-decoding', 'Overfitting floor', 'C');

hK = [hM1, hFam, hR, plot(ax, nan, nan, 'd', 'Color', ink2, 'MarkerFaceColor', 'w', 'MarkerSize', 3.5)];
lg = legend(ax, hK, [{'Model 1 (no shared)', sprintf('Shared, K = %d, \\alpha 0\\rightarrow1', sw.n_latent)}, ...
    arrayfun(@(r) r.name, rec, 'UniformOutput', false), {'recorded, C < 50'}], 'Box', 'off', ...
    'FontSize', fs - 1, 'FontName', font);
lg.Layout.Tile = 'east';
save_fig(fig2, fullfile(out_dir, 'fig2_recorded_vs_partial_sharing'), export_formats);
fprintf('Figures and numbers written to %s\n', out_dir);


%% ======================================================================
%% Local functions
%% ======================================================================
function say(fid, varargin)
% Print to the console and to the numbers file.
fprintf(varargin{:});
fprintf(fid, varargin{:});
end

function a = inv_curve(curve, alpha, y)
% alpha at which a monotone model curve reaches y (linear interpolation; NaN outside its range).
[cs, o] = sort(curve);
a = interp1(cs, alpha(o), y, 'linear', NaN);
end

function C = hexes(h)
C = cell2mat(cellfun(@(x) [hex2dec(x(2:3)) hex2dec(x(4:5)) hex2dec(x(6:7))] / 255, h(:), 'UniformOutput', false));
end

function v = as_row(A)
v = reshape(A, 1, []);
end

function style(ax, font, fs, xl, yl, ttl, letter)
% Journal-style axes: small type, ticks out, no box, panel letter top-left.
set(ax, 'FontName', font, 'FontSize', fs, 'TickDir', 'out', 'LineWidth', 0.75, 'Box', 'off', ...
    'XColor', [0.15 0.15 0.15], 'YColor', [0.15 0.15 0.15], 'TickLength', [0.02 0.02]);
xlabel(ax, xl);  ylabel(ax, yl);
title(ax, ttl, 'FontWeight', 'normal', 'FontSize', fs + 1);
text(ax, -0.24, 1.02, letter, 'Units', 'normalized', 'FontName', font, 'FontSize', fs + 4, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'bottom');
end

function save_fig(fig, stem, formats)
for k = 1:numel(formats)
    switch formats{k}
        case {'pdf', 'eps'}
            exportgraphics(fig, [stem '.' formats{k}], 'ContentType', 'vector');
        otherwise
            exportgraphics(fig, [stem '.' formats{k}], 'Resolution', 300);
    end
end
end

function R = load_real(pdir, name, c)
% Pool the six ordered-pair files of one recorded population (incremental_pr results):
% mean accuracies, PR, ratio, headroom and floor at each neuron count.
pairs = {'acec', 'ecac', 'ecex', 'exec', 'acex', 'exac'};
acc = containers.Map('KeyType', 'double', 'ValueType', 'any');
prs = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:numel(pairs)
    s = load(fullfile(pdir, [pairs{k} '_results.mat']));
    fn = fieldnames(s);  res = s.(fn{1});
    for q = 1:numel(res)
        n = res{q}.neuron_num;
        a = res{q}.acc;  v = [res{q}.pr_stim1(:); res{q}.pr_stim2(:)];
        if isKey(acc, n), acc(n) = [acc(n); a]; prs(n) = [prs(n); v]; else, acc(n) = a; prs(n) = v; end
    end
end
Ns = sort(cell2mat(keys(acc)));
R = struct('name', name, 'N', Ns, 'self', nan(size(Ns)), 'pt', nan(size(Ns)), 'ctrl', nan(size(Ns)), ...
           'ratio', nan(size(Ns)), 'headroom', nan(size(Ns)), 'cs', nan(size(Ns)), 'pr', nan(size(Ns)));
for i = 1:numel(Ns)
    a = mean(acc(Ns(i)), 1);
    R.self(i) = a(c.self);  R.pt(i) = a(c.pt);  R.ctrl(i) = a(c.ctrl);
    R.ratio(i)    = R.pt(i) / R.ctrl(i);
    R.headroom(i) = (R.pt(i) - R.ctrl(i)) / (R.self(i) - R.ctrl(i));
    R.cs(i)       = R.ctrl(i) / R.self(i);
    R.pr(i)       = mean(prs(Ns(i)));
end
end
