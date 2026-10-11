%{
% FILENAME: visualization_simulation_sweep.m   (Simulations_v2)
% AUTHOR:   Zitong Wang
% DATE:     2026-09-17
%
% DESCRIPTION:
%   Figures and a text summary for the N x n_stim simulation sweep produced by
%   procrustes_decoding_simulation_sweep.m (Model 1 = independent random means,
%   the null; Model 3 = linear mixed selectivity, shared rank-2 geometry).
%   Everything is read from the one structured results file the driver writes
%   (results/decoding_outputs/<tag>_results.mat). The file is checkpointed
%   after every grid point, so this also works on a partial run (points not
%   yet done are blank).
%
%   Quantities (columns of pt_transfer_decode_population.m):
%     self     within-cue self-decoding (ceiling)                       col 1
%     PT       rotation-only Procrustes transfer decoding               col 6
%     control  correspondence-shuffled, rotation-only, vs target labels col 4
%              (the overfitting floor; matched to PT)
%     ratio    PT / control per population, then averaged (error = SEM over populations)
%
%   FIGURE 1  sweep_heatmaps_<model>.png      (one per model)
%     2 x 3 heatmaps over N (x) and n_stim (y): self, PT, control (shared 0-1
%     colour scale); PT/control ratio (diverging around 1, log scale);
%     PT - control margin (diverging around 0); empirical PR of the z-scored
%     condition means (the real-data PR definition). Every cell is annotated.
%   FIGURE 2  sweep_curves_vs_N_<model>.png   (one per model)
%     three panels vs N, one curve per n_stim (single-hue ramp, light = few
%     classes, dark = many): (a) self; (b) PT (solid) and control (dashed);
%     (c) PT/control with the 1.0 reference. Error bars = SEM over populations.
%   FIGURE 3  model_comparison_neuron_to_class_ratio.png
%     both models against the reviewer's decisive quantity N/n_stim (log x):
%     (a) control/self, the overfitting floor as a fraction of the ceiling;
%     (b) PT/self; (c) PT/control. Colour = model, one thin line per n_stim.
%   FIGURE 4  real_vs_models_nstim50.png
%     the n_stim = 50 rows of both models beside the four recorded populations
%     (neuronal_decoding incremental_pr results: same decoder, same corrected
%     control, N swept by neuron subsampling): self, PT and control, ratio vs N.
%   TEXT      summary_<tag>.txt   per-model tables of every grid point.
%
% INPUTS:
%   results/decoding_outputs/<tag>_results.mat   (variable `sweep`)
%   optional: ../../../neuronal_decoding/results/decoding_outputs/
%             procrustes_decoding_basic_incremental_pr_results/<M>/<A>/*_results.mat
%
% OUTPUT:
%   results/figures/<tag>/  (all PNGs and the text summary)
%}

%% Configuration
clc; clearvars -except viz_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
sim_root = fullfile(script_dir, '..', '..');

tag       = 'sweep_N_S';   % results file to plot: results/decoding_outputs/<tag>_results.mat
ref_nstim = 50;            % n_stim row compared with the recorded populations in Figure 4

% display names and colours per model (models not listed fall back to their field name)
model_display = struct( ...
    'model1', struct('long',  'Model 1: independent random means (null)', ...
                     'short', 'Model 1 (null)', 'color', [0.00 0.45 0.70]), ...          % blue
    'model3', struct('long',  'Model 3: linear mixed selectivity (shared rank-2 geometry)', ...
                     'short', 'Model 3 (shared geometry)', 'color', [0.84 0.37 0.00]));  % vermillion
col = struct('self',1, 'notf',2, 'pt_full',3, 'ctrl',4, 'scale',5, 'pt',6, 'trans',7, 'chance',8, 'ctrl_id',9);

show_real = true;    % overlay the recorded populations (n_stim = 50) in Figure 4
real_base = fullfile(sim_root, '..', 'neuronal_decoding', 'results', 'decoding_outputs', ...
                     'procrustes_decoding_basic_incremental_pr_results');
real_pops = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};
real_gray = [0.15 0.15 0.15; 0.45 0.45 0.45; 0.15 0.15 0.15; 0.45 0.45 0.45];
real_style = {'-', '-', '--', '--'};                                 % FR solid, KO dashed

% optional overrides from the workspace, e.g. viz_override = struct('tag','sweep_N_S_pc5','ref_nstim',50)
if exist('viz_override', 'var') && isstruct(viz_override)
    fo = fieldnames(viz_override);
    for i = 1:numel(fo), eval(sprintf('%s = viz_override.%s;', fo{i}, fo{i})); end
end

fig_dir = fullfile(sim_root, 'results', 'figures', tag);
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Load the structured results file (points not done yet are NaN)
results_file = fullfile(sim_root, 'results', 'decoding_outputs', sprintf('%s_results.mat', tag));
if ~exist(results_file, 'file')
    error('No results file %s. Run procrustes_decoding_simulation_sweep.m first.', results_file);
end
loaded = load(results_file, 'sweep');
sweep  = loaded.sweep;
D3 = 2;                                             % Model 3 latent dimension (generator default)
if isfield(sweep.cfg, 'gen') && isfield(sweep.cfg.gen, 'n_latent'), D3 = sweep.cfg.gen.n_latent; end
model_display.model3.long  = sprintf('Model 3: linear mixed selectivity (shared %d-D geometry)', D3);
model_display.model3.short = sprintf('Model 3 (shared %d-D geometry)', D3);
models = sweep.models(cellfun(@(m) any(sweep.(m).done(:)), sweep.models));   % models with data
assert(~isempty(models), 'No grid point is done yet in %s.', results_file);
nM = numel(models);
SW = cellfun(@(m) sweep.(m), models, 'UniformOutput', false);
model_names  = cell(1, nM);
model_short  = cell(1, nM);
model_colors = zeros(nM, 3);
fallback     = lines(nM);
for m = 1:nM
    if isfield(model_display, models{m})
        d = model_display.(models{m});
    else
        d = struct('long', models{m}, 'short', models{m}, 'color', fallback(m, :));
    end
    model_names{m} = d.long;  model_short{m} = d.short;  model_colors(m, :) = d.color;
    fprintf('%s: %d / %d grid points done\n', models{m}, nnz(SW{m}.done), numel(SW{m}.done));
end
N_list = sweep.N_list;  S_list = sweep.S_list;
nN = numel(N_list);  nS = numel(S_list);
R_pop  = sweep.cfg.R_pop;
S_colors = seq_ramp([0.35 0.20 0.55], nS);                            % one hue, light -> dark with n_stim

%% FIGURE 1: heatmaps per model
for m = 1:nM
    sw = SW{m};
    self_m  = sw.mean(:,:,col.self);
    pt_m    = sw.mean(:,:,col.pt);
    ctrl_m  = sw.mean(:,:,col.ctrl);
    ratio_m = sw.ratio_mean;
    marg_m  = pt_m - ctrl_m;
    pr_m    = sw.PR_emp;

    figure('Color', 'w', 'Position', [60 60 1500 760]);
    heat_panel(subplot(2,3,1), self_m,  N_list, S_list, 'self-decoding',            [0 1],  seq_ramp([0.10 0.30 0.55], 64), '%.2f');
    heat_panel(subplot(2,3,2), pt_m,    N_list, S_list, 'PT (rotation-only)',       [0 1],  seq_ramp([0.10 0.30 0.55], 64), '%.2f');
    heat_panel(subplot(2,3,3), ctrl_m,  N_list, S_list, 'control (shuffled correspondence)', [0 1], seq_ramp([0.10 0.30 0.55], 64), '%.2f');
    L = max(abs(log2(ratio_m(:))), [], 'omitnan');  L = max(L, 1);
    ax = heat_panel(subplot(2,3,4), log2(ratio_m), N_list, S_list, 'PT / control  (1 = pure overfitting)', [-L L], div_ramp(64), '');
    annotate_cells(ax, ratio_m, '%.2f');
    cb = colorbar(ax); tk = [-L 0 L]; cb.Ticks = tk; cb.TickLabels = arrayfun(@(t) sprintf('%.2g', 2^t), tk, 'UniformOutput', false);
    Lm = max(abs(marg_m(:)), [], 'omitnan');  Lm = max(Lm, 0.05);
    heat_panel(subplot(2,3,5), marg_m,  N_list, S_list, 'PT - control',             [-Lm Lm], div_ramp(64), '%.2f');
    heat_panel(subplot(2,3,6), pr_m,    N_list, S_list, 'participation ratio (z-scored means)', [0 max(pr_m(:),[],'omitnan')], seq_ramp([0.10 0.30 0.55], 64), '%.1f');
    sgtitle(sprintf('%s   |   %d populations per point', model_names{m}, R_pop), 'FontWeight', 'bold');
    saveas(gcf, fullfile(fig_dir, sprintf('sweep_heatmaps_%s.png', models{m})));
end

%% FIGURE 2: curves vs N per model, one line per n_stim
for m = 1:nM
    sw = SW{m};
    figure('Color', 'w', 'Position', [60 60 1500 420]);
    h_leg = gobjects(1, nS);

    ax1 = subplot(1,3,1); hold(ax1, 'on');
    for iS = 1:nS
        h_leg(iS) = errorbar(ax1, N_list, sw.mean(:,iS,col.self), sw.sem(:,iS,col.self), '-o', ...
            'Color', S_colors(iS,:), 'MarkerFaceColor', S_colors(iS,:), 'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
        yline(ax1, 1/S_list(iS), ':', 'Color', S_colors(iS,:), 'HandleVisibility', 'off');
    end
    style_axes(ax1, 'number of neurons N', 'accuracy', 'self-decoding (dotted = chance 1/n_{stim})');
    ylim(ax1, [0 1]);
    legend(ax1, h_leg, arrayfun(@(s) sprintf('n_{stim} = %d', s), S_list, 'UniformOutput', false), 'Location', 'southeast', 'Box', 'off');

    ax2 = subplot(1,3,2); hold(ax2, 'on');
    for iS = 1:nS
        errorbar(ax2, N_list, sw.mean(:,iS,col.pt),   sw.sem(:,iS,col.pt),   '-o',  'Color', S_colors(iS,:), 'MarkerFaceColor', S_colors(iS,:), 'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
        errorbar(ax2, N_list, sw.mean(:,iS,col.ctrl), sw.sem(:,iS,col.ctrl), '--s', 'Color', S_colors(iS,:), 'MarkerFaceColor', 'w',            'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
    end
    hk = [plot(ax2, nan, nan, '-o',  'Color', 'k', 'MarkerFaceColor', 'k', 'LineWidth', 1.5), ...
          plot(ax2, nan, nan, '--s', 'Color', 'k', 'MarkerFaceColor', 'w', 'LineWidth', 1.5)];
    style_axes(ax2, 'number of neurons N', 'accuracy', 'PT (solid) vs control (dashed)');
    ylim(ax2, [0 1]);
    legend(ax2, hk, {'PT (rotation-only)', 'control (shuffled correspondence)'}, 'Location', 'northwest', 'Box', 'off');

    ax3 = subplot(1,3,3); hold(ax3, 'on');
    for iS = 1:nS
        errorbar(ax3, N_list, sw.ratio_mean(:,iS), sw.ratio_sem(:,iS), '-o', ...
            'Color', S_colors(iS,:), 'MarkerFaceColor', S_colors(iS,:), 'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
    end
    yline(ax3, 1, 'k-', 'HandleVisibility', 'off');
    style_axes(ax3, 'number of neurons N', 'PT / control', 'PT relative to its overfitting floor');
    set(ax3, 'YScale', 'log');
    yl = ylim(ax3); ylim(ax3, [min(0.5, yl(1)), max(2, yl(2))]);

    sgtitle(model_names{m}, 'FontWeight', 'bold');
    saveas(gcf, fullfile(fig_dir, sprintf('sweep_curves_vs_N_%s.png', models{m})));
end

%% FIGURE 3: both models against the neuron-to-class ratio N / n_stim
figure('Color', 'w', 'Position', [60 60 1500 420]);
panels = {'control / self', 'PT / self', 'PT / control'};
for k = 1:3
    ax = subplot(1,3,k); hold(ax, 'on');
    hm = gobjects(1, nM);
    for m = 1:nM
        sw = SW{m};
        for iS = 1:nS
            x = N_list / S_list(iS);
            switch k
                case 1, y = sw.mean(:,iS,col.ctrl) ./ sw.mean(:,iS,col.self);
                case 2, y = sw.mean(:,iS,col.pt)   ./ sw.mean(:,iS,col.self);
                case 3, y = sw.ratio_mean(:,iS);
            end
            h = plot(ax, x, y, '-o', 'Color', model_colors(m,:), 'MarkerFaceColor', model_colors(m,:), ...
                     'MarkerSize', 4, 'LineWidth', 1);
            if iS == 1, hm(m) = h; end
        end
    end
    set(ax, 'XScale', 'log');
    if k == 3
        set(ax, 'YScale', 'log');
        yline(ax, 1, 'k-', 'HandleVisibility', 'off');
    else
        ylim(ax, [0 1.05]);
    end
    xline(ax, 1, 'k:', 'N = n_{stim}', 'LabelVerticalAlignment', 'bottom', ...
          'LabelOrientation', 'horizontal', 'HandleVisibility', 'off');
    style_axes(ax, 'N / n_{stim}  (neuron-to-class ratio)', panels{k}, panels{k});
    if k == 1, legend(ax, hm, model_short, 'Location', 'southeast', 'Box', 'off'); end
end
sgtitle('Overfitting floor and PT margin against the neuron-to-class ratio (one line per n_{stim})', 'FontWeight', 'bold');
saveas(gcf, fullfile(fig_dir, 'model_comparison_neuron_to_class_ratio.png'));

%% FIGURE 4: the n_stim = ref_nstim rows of both models beside the recorded populations (50 stimuli)
iS50 = find(S_list == ref_nstim, 1);
if ~isempty(iS50)
    real = struct('name', {}, 'N', {}, 'self', {}, 'self_sem', {}, 'pt', {}, 'pt_sem', {}, 'ctrl', {}, 'ctrl_sem', {}, 'ratio', {});
    if show_real
        for p = 1:size(real_pops,1)
            try
                real(end+1) = load_real_population(fullfile(real_base, real_pops{p,1}, real_pops{p,2}), ...
                    [real_pops{p,1} ' ' real_pops{p,2}], col); %#ok<SAGROW>
            catch ME
                warning('Real data for %s %s not loaded (%s); skipped.', real_pops{p,1}, real_pops{p,2}, ME.message);
            end
        end
    end

    figure('Color', 'w', 'Position', [60 60 1500 420]);
    ax1 = subplot(1,3,1); hold(ax1, 'on');
    ax2 = subplot(1,3,2); hold(ax2, 'on');
    ax3 = subplot(1,3,3); hold(ax3, 'on');
    hl = []; hn = {};
    for m = 1:nM
        sw = SW{m};
        c  = model_colors(m,:);
        h = errorbar(ax1, N_list, sw.mean(:,iS50,col.self), sw.sem(:,iS50,col.self), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
        errorbar(ax2, N_list, sw.mean(:,iS50,col.pt),   sw.sem(:,iS50,col.pt),   '-o',  'Color', c, 'MarkerFaceColor', c,   'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
        errorbar(ax2, N_list, sw.mean(:,iS50,col.ctrl), sw.sem(:,iS50,col.ctrl), '--s', 'Color', c, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
        errorbar(ax3, N_list, sw.ratio_mean(:,iS50), sw.ratio_sem(:,iS50), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 3);
        hl(end+1) = h; hn{end+1} = model_short{m}; %#ok<SAGROW>
    end
    for p = 1:numel(real)
        c = real_gray(p,:); ls = real_style{p};
        h = errorbar(ax1, real(p).N, real(p).self, real(p).self_sem, [ls 'd'], 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 5, 'LineWidth', 1.2, 'CapSize', 3);
        errorbar(ax2, real(p).N, real(p).pt,   real(p).pt_sem,   [ls 'd'], 'Color', c, 'MarkerFaceColor', c,   'MarkerSize', 5, 'LineWidth', 1.2, 'CapSize', 3);
        errorbar(ax2, real(p).N, real(p).ctrl, real(p).ctrl_sem, [ls 's'], 'Color', c, 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'LineWidth', 1.2, 'CapSize', 3);
        plot(ax3, real(p).N, real(p).ratio, [ls 'd'], 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 5, 'LineWidth', 1.2);
        hl(end+1) = h; hn{end+1} = ['recorded ' real(p).name]; %#ok<SAGROW>
    end
    style_axes(ax1, 'number of neurons N', 'accuracy', sprintf('self-decoding, n_{stim} = %d', ref_nstim));
    ylim(ax1, [0 1]); yline(ax1, 1/ref_nstim, 'k:', 'chance', 'HandleVisibility', 'off');
    legend(ax1, hl, hn, 'Location', 'northwest', 'Box', 'off');
    hk = [plot(ax2, nan, nan, '-', 'Color', 'k', 'LineWidth', 1.5), plot(ax2, nan, nan, '--', 'Color', 'k', 'LineWidth', 1.5)];
    style_axes(ax2, 'number of neurons N', 'accuracy', sprintf('PT (solid, filled) vs control (dashed, open), n_{stim} = %d', ref_nstim));
    ylim(ax2, [0 1]);
    legend(ax2, hk, {'PT (rotation-only)', 'control'}, 'Location', 'northwest', 'Box', 'off');
    yline(ax3, 1, 'k-', 'HandleVisibility', 'off');
    style_axes(ax3, 'number of neurons N', 'PT / control', sprintf('PT relative to its overfitting floor, n_{stim} = %d', ref_nstim));
    set(ax3, 'YScale', 'log');
    sgtitle(sprintf('Recorded populations (50 stimuli) vs the two synthetic models at n_{stim} = %d (same decoder, same corrected control)', ref_nstim), 'FontWeight', 'bold');
    saveas(gcf, fullfile(fig_dir, sprintf('real_vs_models_nstim%d.png', ref_nstim)));
end

%% TEXT SUMMARY
fid = fopen(fullfile(fig_dir, sprintf('summary_%s.txt', tag)), 'w');
for m = 1:nM
    sw = SW{m};
    fprintf(fid, '===== %s  (%s; %d populations per point; mean over populations, SEM in brackets) =====\n', model_names{m}, tag, R_pop);
    for iS = 1:nS
        fprintf(fid, '\n-- n_stim = %d (chance %.3f) --\n', S_list(iS), 1/S_list(iS));
        fprintf(fid, '%5s  %13s  %13s  %13s  %13s  %13s  %8s  %7s\n', 'N', 'self', 'no-transform', 'PT(rot)', 'control', 'PT/control', 'ctrl-id', 'PR_emp');
        for iN = 1:nN
            if ~sw.done(iN,iS), fprintf(fid, '%5d  (not done)\n', N_list(iN)); continue; end
            fprintf(fid, '%5d  %5.3f (%.3f)  %5.3f (%.3f)  %5.3f (%.3f)  %5.3f (%.3f)  %5.2f (%.2f)   %6.3f  %7.1f\n', N_list(iN), ...
                sw.mean(iN,iS,col.self), sw.sem(iN,iS,col.self), sw.mean(iN,iS,col.notf), sw.sem(iN,iS,col.notf), ...
                sw.mean(iN,iS,col.pt), sw.sem(iN,iS,col.pt), sw.mean(iN,iS,col.ctrl), sw.sem(iN,iS,col.ctrl), ...
                sw.ratio_mean(iN,iS), sw.ratio_sem(iN,iS), sw.mean(iN,iS,col.ctrl_id), sw.PR_emp(iN,iS));
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);
fprintf('Figures and summary written to %s\n', fig_dir);
type(fullfile(fig_dir, sprintf('summary_%s.txt', tag)));


%% ======================================================================
%% Local functions
%% ======================================================================
function ax = heat_panel(ax, M, N_list, S_list, ttl, clim, cmap, fmt)
% Heatmap of a [nN x nS] grid quantity with N on x and n_stim on y, annotated.
im = imagesc(ax, 1:numel(N_list), 1:numel(S_list), M.');
im.AlphaData = double(~isnan(M.'));                  % missing grid points stay blank
set(ax, 'YDir', 'normal', 'XTick', 1:numel(N_list), 'XTickLabel', N_list, ...
        'YTick', 1:numel(S_list), 'YTickLabel', S_list, 'TickLength', [0 0], 'FontSize', 9);
colormap(ax, cmap);
if ~isempty(clim) && clim(2) > clim(1), caxis(ax, clim); end
colorbar(ax);
xlabel(ax, 'number of neurons N'); ylabel(ax, 'number of stimuli n_{stim}');
title(ax, ttl);
if ~isempty(fmt), annotate_cells(ax, M, fmt); end
end

function annotate_cells(ax, M, fmt)
% Write each cell's value in the cell; ink colour follows the cell's brightness.
im = findobj(ax, 'Type', 'image');
cd = im.CData;  cl = caxis(ax);  cmap = colormap(ax);
for i = 1:size(M,1)
    for j = 1:size(M,2)
        if isnan(M(i,j)), continue; end
        v  = (cd(j,i) - cl(1)) / (cl(2) - cl(1));
        ci = max(1, min(size(cmap,1), round(1 + v * (size(cmap,1) - 1))));
        lum = cmap(ci,:) * [0.299 0.587 0.114].';
        if lum < 0.5, ink = [1 1 1]; else, ink = [0 0 0]; end
        text(ax, i, j, sprintf(fmt, M(i,j)), 'HorizontalAlignment', 'center', ...
             'FontSize', 8, 'Color', ink);
    end
end
end

function cmap = seq_ramp(dark, n)
% Single-hue sequential ramp: near-white -> `dark`.
light = 1 - 0.12 * (1 - dark);
t = linspace(0, 1, n).';
cmap = (1 - t) .* light + t .* dark;
end

function cmap = div_ramp(n)
% Diverging ramp: blue (low) -> neutral gray (mid) -> vermillion (high).
lo = [0.00 0.45 0.70];  mid = [0.93 0.93 0.93];  hi = [0.84 0.37 0.00];
h = floor(n/2);
t = linspace(0, 1, h).';
cmap = [(1 - t) .* lo + t .* mid; (1 - t) .* mid + t .* hi];
end

function style_axes(ax, xl, yl, ttl)
xlabel(ax, xl); ylabel(ax, yl); title(ax, ttl);
box(ax, 'off'); grid(ax, 'on'); set(ax, 'GridAlpha', 0.15, 'FontSize', 10);
end

function real = load_real_population(pdir, name, col)
% Pool the six ordered-pair files of one recorded population (incremental_pr
% results): mean and SEM over all resampling rows at each neuron count.
pairs = {'acec','ecac','ecex','exec','acex','exac'};
acc_by_N = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:numel(pairs)
    f = fullfile(pdir, [pairs{k} '_results.mat']);
    if ~exist(f, 'file'), continue; end
    r = load(f); fn = fieldnames(r); res = r.(fn{1});
    for q = 1:numel(res)
        n = res{q}.neuron_num;
        if isKey(acc_by_N, n), acc_by_N(n) = [acc_by_N(n); res{q}.acc]; else, acc_by_N(n) = res{q}.acc; end
    end
end
assert(acc_by_N.Count > 0, 'no pair files in %s', pdir);
Ns = sort(cell2mat(keys(acc_by_N)));
real = struct('name', name, 'N', Ns, 'self', nan(size(Ns)), 'self_sem', nan(size(Ns)), ...
              'pt', nan(size(Ns)), 'pt_sem', nan(size(Ns)), 'ctrl', nan(size(Ns)), 'ctrl_sem', nan(size(Ns)), 'ratio', nan(size(Ns)));
for i = 1:numel(Ns)
    A = acc_by_N(Ns(i));  n = size(A,1);
    real.self(i) = mean(A(:,col.self)); real.self_sem(i) = std(A(:,col.self))/sqrt(n);
    real.pt(i)   = mean(A(:,col.pt));   real.pt_sem(i)   = std(A(:,col.pt))/sqrt(n);
    real.ctrl(i) = mean(A(:,col.ctrl)); real.ctrl_sem(i) = std(A(:,col.ctrl))/sqrt(n);
    real.ratio(i) = real.pt(i) / real.ctrl(i);
end
end
