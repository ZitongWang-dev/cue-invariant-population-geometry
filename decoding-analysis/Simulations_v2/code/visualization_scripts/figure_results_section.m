%{
% FILENAME: figure_results_section.m   (Simulations_v2)
% AUTHOR:   Zitong Wang
% DATE:     2026-09-24 (2026-10-07: each model read from its own results file;
%           Model 3 = 5-D shared latent)
%
% DESCRIPTION:
%   Paper-style figures, and every number quoted, for the results write-up of
%   the simulation benchmark (Model 1 = independent random means, Model 3 =
%   shared low-dimensional geometry; N x n_stim sweep) and its comparison with
%   the recorded populations.
%
%   Each model is read from its own results file (`sources` below):
%     Model 1 : sweep_N_S             n_stim 5-50
%     Model 3 : sweep_N_S_model3_D5   5-D shared latent, n_stim 10-50 (a 5-D
%               latent needs at least 6 stimuli). Point the Model 3 entry at
%               'sweep_N_S' to plot the earlier 2-D version instead.
%   The Model 3 run reuses Model 1's seeds, so each Model 3 population has the
%   same neurons as the Model 1 population at the same grid point.
%
%   FIGURE 1  fig1_synthetic_benchmark.png
%     A  PT vs control at every grid point, both models (log axes), with the
%        PT = control, 2x and 4x lines.
%     B  control / self-decoding vs N / n_stim, one line per n_stim (darker =
%        more classes): the overfitting floor as a fraction of the ceiling.
%     C  PT / control vs n_stim (dots = individual N, line = mean over N).
%   FIGURE 2  fig2_recorded_vs_models.png   (50 stimulus classes)
%     A  self-decoding, B  PT / control, C  participation ratio, all vs N, for
%     both models and the four recorded populations (incremental_pr results:
%     same decoder, same correspondence-shuffled control, N by subsampling).
%   TEXT      results_section_numbers.txt   numbers quoted in the write-up,
%             including Table 1 (50 classes, full population size).
%
%   Ratios are ratios of mean accuracies (mean PT / mean control), for models
%   and recorded data alike. PR = participation ratio of the z-scored,
%   trial-averaged condition means (same definition for models and data).
%
% INPUTS:
%   results/decoding_outputs/<tag>_results.mat   (one per entry of `sources`)
%   ../../../neuronal_decoding/results/decoding_outputs/
%       procrustes_decoding_basic_incremental_pr_results/<M>/<A>/*_results.mat
%
% OUTPUT:
%   results/figures/results_section/
%}

%% Configuration
clc; clear;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
sim_root = fullfile(script_dir, '..', '..');
sources  = {'model1', 'sweep_N_S';                 % {model, results tag}
            'model3', 'sweep_N_S_model3_D5'};
out_dir  = fullfile(sim_root, 'results', 'figures', 'results_section');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
export_formats = {'png'};          % add 'tif', 'pdf' or 'eps' for submission
font = 'Arial';  fs = 8;

% colours: categorical slots 1-2 (blue, orange); recorded data in neutral ink
c_m1  = [42 120 214] / 255;        % Model 1
c_m3  = [235 104 52] / 255;        % Model 3
c_ref = [0.55 0.55 0.53];          % reference lines
ink2  = [0.32 0.32 0.30];          % secondary text
real_pops = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};
real_col  = [0.10 0.10 0.10; 0.58 0.58 0.58; 0.10 0.10 0.10; 0.58 0.58 0.58];   % V1 dark, V2 light
real_ls   = {'-', '-', '--', '--'};                                               % FR solid, KO dashed
real_mk   = {'d', 'd', 's', 's'};
c = struct('self', 1, 'notf', 2, 'ctrl', 4, 'pt', 6, 'ctrl_id', 9);

%% Load
m1 = load_model(sim_root, sources(1,:), c);
m3 = load_model(sim_root, sources(2,:), c);
lab1 = 'Model 1: no shared structure';
lab3 = sprintf('Model 3: shared %d-D geometry', m3.D);

real_base = fullfile(sim_root, '..', 'neuronal_decoding', 'results', 'decoding_outputs', ...
                     'procrustes_decoding_basic_incremental_pr_results');
real = [];
for p = 1:size(real_pops, 1)
    real = [real, load_real(fullfile(real_base, real_pops{p,1}, real_pops{p,2}), ...
                            [real_pops{p,1} ' ' real_pops{p,2}])]; %#ok<AGROW>
end

%% Numbers quoted in the write-up
fid = fopen(fullfile(out_dir, 'results_section_numbers.txt'), 'w');
say(fid, 'Numbers for the simulation results write-up. Model 1 from %s (%d populations per point);\n', m1.tag, m1.R);
say(fid, 'Model 3 from %s (%d-D latent, %d populations per point). Ratios = mean PT / mean control.\n\n', m3.tag, m3.D, m3.R);

say(fid, '--- Model 1 (no shared structure) ---\n');
say(fid, 'PT/control, all %d points: median %.2f, range %.2f-%.2f; N >= 20: %.2f-%.2f\n', numel(m1.ratio), ...
    median(m1.ratio(:)), min(m1.ratio(:)), max(m1.ratio(:)), min(m1.ratio(m1.NN >= 20)), max(m1.ratio(m1.NN >= 20)));
say(fid, 'control scored vs permuted identities: max |acc - 1/n_stim| = %.3f\n', max(abs(m1.ctrl_id - 1 ./ m1.SS), [], 'all'));
for N = [5 50 100]
    iN = find(m1.N == N);
    say(fid, '50 classes, N = %3d: self %.3f  PT %.3f  control %.3f\n', N, m1.self(iN,m1.i50), m1.pt(iN,m1.i50), m1.ctrl(iN,m1.i50));
end
say(fid, '50 classes: self range %.3f-%.3f, PR %.1f-%.1f over N\n', min(m1.self(:,m1.i50)), max(m1.self(:,m1.i50)), ...
    min(m1.PR(:,m1.i50)), max(m1.PR(:,m1.i50)));
q = m1.NN ./ m1.SS;
bins = {q <= 0.25, q == 1 & m1.SS >= 20, q == 1 & m1.SS < 20, q == 2, q >= 6};
bin_names = {'N/n_stim <= 0.25', 'N/n_stim = 1 (n_stim >= 20)', 'N/n_stim = 1 (n_stim < 20)', 'N/n_stim = 2', 'N/n_stim >= 6'};
for b = 1:numel(bins)
    say(fid, 'control/self, %-28s: %.2f-%.2f (%d points)\n', bin_names{b}, min(m1.cs(bins{b})), max(m1.cs(bins{b})), nnz(bins{b}));
end

say(fid, '\n--- both models ---\n');
say(fid, 'no-transform: max |acc - 1/n_stim| = %.3f (Model 1), %.3f (Model 3)\n', ...
    max(abs(m1.notf - 1 ./ m1.SS), [], 'all'), max(abs(m3.notf - 1 ./ m3.SS), [], 'all'));

say(fid, '\n--- Model 3 (shared %d-D geometry) ---\n', m3.D);
say(fid, 'PT/self over the grid: %.2f-%.2f\n', min(m3.pt(:) ./ m3.self(:)), max(m3.pt(:) ./ m3.self(:)));
for iS = 1:numel(m3.S)
    sel = m3.N >= 20;
    say(fid, 'n_stim = %2d: PT/control mean over N %.2f (all N %.2f-%.2f; N >= 20 %.2f-%.2f) | control/self %.2f-%.2f (N = 100: %.2f) | self %.3f-%.3f\n', ...
        m3.S(iS), mean(m3.ratio(:,iS)), min(m3.ratio(:,iS)), max(m3.ratio(:,iS)), min(m3.ratio(sel,iS)), max(m3.ratio(sel,iS)), ...
        min(m3.cs(:,iS)), max(m3.cs(:,iS)), m3.cs(end,iS), min(m3.self(:,iS)), max(m3.self(:,iS)));
end
say(fid, '50 classes, N = 100: self %.3f  PT %.3f  control %.3f\n', m3.self(end,m3.i50), m3.pt(end,m3.i50), m3.ctrl(end,m3.i50));
say(fid, '50 classes: self range %.3f-%.3f, PR %.1f-%.1f over N\n', min(m3.self(:,m3.i50)), max(m3.self(:,m3.i50)), ...
    min(m3.PR(:,m3.i50)), max(m3.PR(:,m3.i50)));

say(fid, '\n--- recorded populations (50 classes) ---\n');
for p = 1:numel(real)
    say(fid, '%s: N = %d..%d | PT/control %.2f-%.2f over N | PR %.1f-%.1f | self %.3f-%.3f\n', real(p).name, ...
        real(p).N(1), real(p).N(end), min(real(p).ratio), max(real(p).ratio), min(real(p).pr), max(real(p).pr), ...
        min(real(p).self), max(real(p).self));
end
lo = min(arrayfun(@(r) r.self(end), real));  hi = max(arrayfun(@(r) r.self(end), real));
for mm = {m1, m3}
    m = mm{1};
    match = find(m.self(:,m.i50) >= lo - 0.05 & m.self(:,m.i50) <= hi + 0.05);
    say(fid, '%s at 50 classes with self in the recorded range (%.2f-%.2f, +/-0.05):\n', m.name, lo, hi);
    for i = match.'
        say(fid, '   N = %3d: self %.3f  PT/control %.2f  PR %.1f\n', m.N(i), m.self(i,m.i50), m.ratio(i,m.i50), m.PR(i,m.i50));
    end
end

say(fid, '\nTable 1 (50 classes, full population size)\n');
say(fid, '| Population | N | Self-decoding | PT | Control | PT / control | PR |\n| --- | --- | --- | --- | --- | --- | --- |\n');
for p = 1:numel(real)
    say(fid, '| %s | %d | %.2f | %.2f | %.2f | %.2f | %.1f |\n', real(p).name, real(p).N(end), real(p).self(end), ...
        real(p).pt(end), real(p).ctrl(end), real(p).ratio(end), real(p).pr(end));
end
for mm = {m1, m3}
    m = mm{1};
    say(fid, '| %s | %d | %.2f | %.2f | %.2f | %.2f | %.1f |\n', m.name, m.N(end), m.self(end,m.i50), m.pt(end,m.i50), ...
        m.ctrl(end,m.i50), m.ratio(end,m.i50), m.PR(end,m.i50));
end
fclose(fid);

%% FIGURE 1: synthetic benchmark
fig1 = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 7.2 2.6]);
t1 = tiledlayout(fig1, 1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');

% A: PT vs control, every grid point
ax = nexttile(t1);  hold(ax, 'on');
lims = [0.02 1];
xx = logspace(log10(lims(1)), log10(lims(2)), 50);
plot(ax, xx, xx, '-', 'Color', c_ref, 'LineWidth', 0.75);
plot(ax, xx, 2 * xx, ':', 'Color', c_ref, 'LineWidth', 0.75);
plot(ax, xx, 4 * xx, ':', 'Color', c_ref, 'LineWidth', 0.75);
text(ax, 0.75, 0.92, '\times1', 'Color', ink2, 'FontSize', fs - 1, 'FontName', font, 'HorizontalAlignment', 'right');
text(ax, 0.40, 0.92, '\times2', 'Color', ink2, 'FontSize', fs - 1, 'FontName', font, 'HorizontalAlignment', 'right');
text(ax, 0.20, 0.92, '\times4', 'Color', ink2, 'FontSize', fs - 1, 'FontName', font, 'HorizontalAlignment', 'right');
h1 = scatter(ax, m1.ctrl(:), m1.pt(:), 16, c_m1, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.3);
h3 = scatter(ax, m3.ctrl(:), m3.pt(:), 16, c_m3, 'filled', 'MarkerEdgeColor', 'w', 'LineWidth', 0.3);
set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', lims, 'YLim', lims, ...
    'XTick', [0.02 0.05 0.1 0.2 0.5 1], 'YTick', [0.02 0.05 0.1 0.2 0.5 1]);
style(ax, font, fs, 'Control accuracy', 'PT accuracy', 'PT vs. control', 'A');
lg1 = legend(ax, [h1 h3], {lab1, lab3}, 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', fs, 'FontName', font);
lg1.Layout.Tile = 'north';

% B: overfitting floor vs neurons per class
ax = nexttile(t1);  hold(ax, 'on');
ramp1 = one_hue_ramp(c_m1, numel(m1.S));  ramp3 = one_hue_ramp(c_m3, numel(m3.S));
xline(ax, 1, ':', 'Color', c_ref, 'LineWidth', 0.75);
for iS = 1:numel(m1.S)
    plot(ax, m1.N / m1.S(iS), m1.cs(:,iS), '-o', 'Color', ramp1(iS,:), 'MarkerFaceColor', ramp1(iS,:), 'MarkerSize', 2.5, 'LineWidth', 1);
end
for iS = 1:numel(m3.S)
    plot(ax, m3.N / m3.S(iS), m3.cs(:,iS), '-o', 'Color', ramp3(iS,:), 'MarkerFaceColor', ramp3(iS,:), 'MarkerSize', 2.5, 'LineWidth', 1);
end
text(ax, m3.N(end) / m3.S(1), m3.cs(end,1) - 0.13, sprintf('%d classes', m3.S(1)), 'Color', ink2, 'FontSize', fs - 1, ...
    'FontName', font, 'HorizontalAlignment', 'right');
text(ax, m3.N(end) / m3.S(end) * 1.12, m3.cs(end,end) - 0.06, sprintf('%d classes', m3.S(end)), 'Color', ink2, ...
    'FontSize', fs - 1, 'FontName', font);
text(ax, 1.08, 0.06, 'N = n_{stim}', 'Color', ink2, 'FontSize', fs - 1, 'FontName', font);
set(ax, 'XScale', 'log', 'XLim', [0.09 30], 'YLim', [0 1.05], 'XTick', [0.1 1 10], 'XTickLabel', {'0.1', '1', '10'});
style(ax, font, fs, 'Neurons per class (N / n_{stim})', 'Control / self-decoding', 'Overfitting floor', 'B');

% C: PT relative to the floor vs number of classes
ax = nexttile(t1);  hold(ax, 'on');
yline(ax, 1, '-', 'Color', c_ref, 'LineWidth', 0.75);
for mm = {{m1, c_m1}, {m3, c_m3}}
    m = mm{1}{1};  col = mm{1}{2};
    jit = linspace(-1.8, 1.8, numel(m.N));
    for iS = 1:numel(m.S)
        scatter(ax, m.S(iS) + jit, m.ratio(:,iS), 8, col, 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');
    end
    plot(ax, m.S, mean(m.ratio, 1), '-o', 'Color', col, 'MarkerFaceColor', col, 'MarkerSize', 3.5, 'LineWidth', 1.5);
end
set(ax, 'YScale', 'log', 'YLim', [0.8 8], 'YTick', [1 2 4 8], 'XLim', [0 55], 'XTick', union(m1.S, m3.S));
style(ax, font, fs, 'Number of stimulus classes', 'PT / control', 'PT relative to the floor', 'C');

save_fig(fig1, fullfile(out_dir, 'fig1_synthetic_benchmark'), export_formats);

%% FIGURE 2: recorded populations vs the two models at 50 classes
fig2 = figure('Color', 'w', 'Units', 'inches', 'Position', [0.5 0.5 7.2 2.6]);
t2 = tiledlayout(fig2, 1, 3, 'TileSpacing', 'loose', 'Padding', 'compact');
mk = {'MarkerSize', 3, 'LineWidth', 1.25};
hm = gobjects(1, 2 + numel(real));

% A: self-decoding
ax = nexttile(t2);  hold(ax, 'on');
yline(ax, 1/50, ':', 'Color', c_ref, 'LineWidth', 0.75);
errorbar(ax, m1.N, m1.self(:,m1.i50), m1.self_sem(:,m1.i50), '-o', 'Color', c_m1, 'MarkerFaceColor', c_m1, 'CapSize', 0, mk{:});
errorbar(ax, m3.N, m3.self(:,m3.i50), m3.self_sem(:,m3.i50), '-o', 'Color', c_m3, 'MarkerFaceColor', c_m3, 'CapSize', 0, mk{:});
for p = 1:numel(real)
    errorbar(ax, real(p).N, real(p).self, real(p).self_sem, [real_ls{p} real_mk{p}], ...
        'Color', real_col(p,:), 'MarkerFaceColor', real_col(p,:), 'CapSize', 0, mk{:});
end
text(ax, 148, 0.05, 'chance', 'Color', ink2, 'FontSize', fs - 1, 'FontName', font, 'HorizontalAlignment', 'right');
set(ax, 'XLim', [0 150], 'YLim', [0 1], 'XTick', 0:50:150);
style(ax, font, fs, 'Number of neurons N', 'Self-decoding accuracy', 'Self-decoding', 'A');

% B: PT / control (these line handles also feed the legend: no error-bar glyphs)
ax = nexttile(t2);  hold(ax, 'on');
yline(ax, 1, '-', 'Color', c_ref, 'LineWidth', 0.75);
hm(1) = plot(ax, m1.N, m1.ratio(:,m1.i50), '-o', 'Color', c_m1, 'MarkerFaceColor', c_m1, mk{:});
hm(2) = plot(ax, m3.N, m3.ratio(:,m3.i50), '-o', 'Color', c_m3, 'MarkerFaceColor', c_m3, mk{:});
for p = 1:numel(real)
    hm(2 + p) = plot(ax, real(p).N, real(p).ratio, [real_ls{p} real_mk{p}], 'Color', real_col(p,:), ...
        'MarkerFaceColor', real_col(p,:), mk{:});
end
set(ax, 'XLim', [0 150], 'XTick', 0:50:150, 'YScale', 'log', 'YLim', [0.8 8], 'YTick', [1 1.5 2 4 8]);
style(ax, font, fs, 'Number of neurons N', 'PT / control', 'PT relative to the floor', 'B');

% C: participation ratio (models drawn last: Model 3 lies on the recorded curves)
ax = nexttile(t2);  hold(ax, 'on');
for p = 1:numel(real)
    plot(ax, real(p).N, real(p).pr, [real_ls{p} real_mk{p}], 'Color', real_col(p,:), 'MarkerFaceColor', real_col(p,:), mk{:});
end
plot(ax, m1.N, m1.PR(:,m1.i50), '-o', 'Color', c_m1, 'MarkerFaceColor', c_m1, mk{:});
plot(ax, m3.N, m3.PR(:,m3.i50), '-o', 'Color', c_m3, 'MarkerFaceColor', c_m3, mk{:});
set(ax, 'XLim', [0 150], 'XTick', 0:50:150, 'YLim', [0 32]);
style(ax, font, fs, 'Number of neurons N', 'Participation ratio', 'Dimensionality', 'C');

lg = legend(ax, hm, [{'Model 1', sprintf('Model 3 (%d-D)', m3.D)}, arrayfun(@(r) r.name, real, 'UniformOutput', false)], ...
    'Box', 'off', 'FontSize', fs - 1, 'FontName', font);
lg.Layout.Tile = 'east';

save_fig(fig2, fullfile(out_dir, 'fig2_recorded_vs_models'), export_formats);
fprintf('Figures and numbers written to %s\n', out_dir);


%% ======================================================================
%% Local functions
%% ======================================================================
function m = load_model(sim_root, src, c)
% One model's grid arrays from its own results file.
L  = load(fullfile(sim_root, 'results', 'decoding_outputs', [src{2} '_results.mat']), 'sweep');
sw = L.sweep;  d = sw.(src{1});
m = struct('name', src{1}, 'tag', src{2}, 'N', sw.N_list, 'S', sw.S_list, 'R', sw.cfg.R_pop, 'D', NaN);
if isfield(sw.cfg.gen, 'n_latent')
    m.D = sw.cfg.gen.n_latent;
elseif strcmp(src{1}, 'model3')
    m.D = 2;                                       % generator default when not overridden
end
m.self     = d.mean(:,:,c.self);   m.self_sem = d.sem(:,:,c.self);
m.pt       = d.mean(:,:,c.pt);     m.ctrl     = d.mean(:,:,c.ctrl);
m.notf     = d.mean(:,:,c.notf);   m.ctrl_id  = d.mean(:,:,c.ctrl_id);
m.ratio    = d.ratio_of_means;     m.cs       = m.ctrl ./ m.self;
m.PR       = d.PR_emp;
[m.NN, m.SS] = ndgrid(m.N, m.S);
m.i50 = find(m.S == 50);
end

function say(fid, varargin)
% Print to the console and to the numbers file.
fprintf(varargin{:});
fprintf(fid, varargin{:});
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

function cmap = one_hue_ramp(base, n)
% One-hue ramp, light -> full strength, for n ordered lines (never near-white).
t = linspace(0.35, 1, n).';
cmap = (1 - t) .* [1 1 1] + t .* base;
cmap = cmap .* (1 - 0.25 * (t - 0.35) / 0.65);      % darken the far end slightly
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

function R = load_real(pdir, name)
% Pool the six ordered-pair files of one recorded population (incremental_pr
% results): mean/SEM over all resampling rows and mean PR at each neuron count.
pairs = {'acec', 'ecac', 'ecex', 'exec', 'acex', 'exac'};
acc = containers.Map('KeyType', 'double', 'ValueType', 'any');
prs = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:numel(pairs)
    s = load(fullfile(pdir, [pairs{k} '_results.mat']));
    fn = fieldnames(s);  res = s.(fn{1});
    for q = 1:numel(res)
        n = res{q}.neuron_num;
        a = res{q}.acc;  v = [res{q}.pr_stim1(:); res{q}.pr_stim2(:)];
        if isKey(acc, n)
            acc(n) = [acc(n); a];  prs(n) = [prs(n); v];
        else
            acc(n) = a;  prs(n) = v;
        end
    end
end
Ns = sort(cell2mat(keys(acc)));
R = struct('name', name, 'N', Ns, 'self', nan(size(Ns)), 'self_sem', nan(size(Ns)), ...
           'pt', nan(size(Ns)), 'ctrl', nan(size(Ns)), 'ratio', nan(size(Ns)), 'pr', nan(size(Ns)));
for i = 1:numel(Ns)
    A = acc(Ns(i));  n = size(A, 1);
    R.self(i)     = mean(A(:,1));
    R.self_sem(i) = std(A(:,1)) / sqrt(n);
    R.pt(i)       = mean(A(:,6));
    R.ctrl(i)     = mean(A(:,4));
    R.ratio(i)    = R.pt(i) / R.ctrl(i);
    R.pr(i)       = mean(prs(Ns(i)));
end
end
