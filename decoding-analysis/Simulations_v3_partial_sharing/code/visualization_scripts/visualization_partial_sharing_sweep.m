%{
% FILENAME: visualization_partial_sharing_sweep.m   (Simulations_v3_partial_sharing)
% AUTHOR:   Zitong Wang
% DATE:     2026-10-08
%
% DESCRIPTION:
%   Figures and a text summary for the partially shared sweep
%   (procrustes_decoding_partial_sharing_sweep.m: alpha = K_s/K x N x n_stim),
%   set beside the recorded populations. Counterpart of
%   Simulations_v2/code/visualization_scripts/visualization_simulation_sweep.m.
%   Works on partial results (points not done yet are blank).
%
%   Quantities (columns of pt_transfer_decode_population.m):
%     self      within-cue self-decoding (ceiling)                        col 1
%     PT        rotation-only Procrustes transfer decoding                col 6
%     control   correspondence-shuffled, rotation-only, vs target labels  col 4
%     ratio     PT / control
%     headroom  (PT - control) / (self - control): share of the gap between the
%               shuffled floor and self-decoding that PT recovers (0 = nothing
%               shared, 1 = full transfer). Lines use the mean accuracies
%               (headroom_of_means); error bars are the SEM over populations of
%               the per-population headroom.
%
%   FIGURE 1  sweep_heatmaps_alpha<100*alpha>.png   (one per alpha)
%     2 x 3 heatmaps over N (x) and n_stim (y): self, PT, control (0-1),
%     PT/control (log, diverging around 1), headroom (0-1), PR of the z-scored
%     condition means. Every cell annotated.
%   FIGURE 2  sweep_curves_vs_N.png
%     columns = n_stim; rows: (1) PT per alpha (solid), self (black) and control
%     (gray dashed; both averaged over alpha, which they do not depend on);
%     (2) PT/control per alpha; (3) headroom per alpha.
%   FIGURE 3  alpha_calibration.png
%     against alpha, one line per n_stim: (a) headroom and (b) PT/control
%     (means over N >= min_N_avg); (c) self-decoding and (d) PR at N = 100,
%     which must not depend on alpha (fixed-K design).
%   FIGURE 4  real_vs_alpha_nstim50.png
%     50 classes, against N: headroom, PT/control and control/self for every
%     alpha, Model 1 (Simulations_v2, floor panel) and the four recorded
%     populations (neuronal_decoding incremental_pr results, neuron subsampling).
%   FIGURE 5  real_stimsub_vs_alpha.png   (when real_stimsub_results.mat exists)
%     recorded headroom against the number of classes C (stimulus subsampling,
%     real_data_stimulus_subsampling.m) for both modes (all neurons; N / C held
%     at the recorded ratio), with each alpha's model headroom (mean over
%     N >= min_N_avg) as reference lines; and the floor (control / self) against
%     N / C for the alpha family, Model 1 and the recorded subsets.
%   TEXT      summary_<tag>.txt   every grid point per alpha, and the recorded
%             stimulus-subsampling table.
%
%   Colours: alpha = one-hue orange ramp, light (alpha = 0) -> dark (alpha = 1,
%   the colour of Model 3 in Simulations_v2); n_stim = one-hue purple ramp (as v2);
%   both ramps pass the dataviz ordinal checks (monotone lightness, adjacent
%   delta-L >= 0.06, light end >= 2:1 on white, one hue). Model 1 = v2 blue.
%   Recorded populations in neutral ink: V1 dark, V2 light, FR solid, KO dashed.
%
% INPUTS:
%   results/decoding_outputs/<tag>_results.mat                 (variable `sweep`)
%   results/decoding_outputs/<real_tag>_results.mat            (variable `real`, optional)
%   ../Simulations_v2/results/decoding_outputs/sweep_N_S_results.mat   (Model 1)
%   ../neuronal_decoding/results/decoding_outputs/procrustes_decoding_basic_incremental_pr_results/<M>/<A>/
%
% OUTPUT:
%   results/figures/<tag>/
%}

%% Configuration
clc; clearvars -except viz_override;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
sim_root = fullfile(script_dir, '..', '..');

tag       = 'sweep_alpha_N_S';   % results/decoding_outputs/<tag>_results.mat
real_tag  = 'real_stimsub';      % recorded stimulus subsampling (optional)
ref_nstim = 50;                  % n_stim row compared with the recorded populations
min_N_avg = 20;                  % averages over N use N >= min_N_avg
ref_N     = 100;                 % N used for the within-cue invariance panels
col = struct('self',1, 'notf',2, 'pt_full',3, 'ctrl',4, 'scale',5, 'pt',6, 'trans',7, 'chance',8, 'ctrl_id',9);

alpha_hex = {'#ef9b74', '#e87c4c', '#dc5d28', '#b64919', '#8e3812', '#62260a'};   % alpha 0 -> 1
nstim_hex = {'#b9a7cf', '#9c84bd', '#7f62a9', '#634490', '#4a2f70'};              % n_stim 10 -> 50
c_m1      = [42 120 214] / 255;                                                  % Model 1 (v2)
c_ref     = [0.55 0.55 0.53];
real_pops  = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};
real_gray  = [0.10 0.10 0.10; 0.58 0.58 0.58; 0.10 0.10 0.10; 0.58 0.58 0.58];
real_style = {'-', '-', '--', '--'};
real_mk    = {'d', 'd', 's', 's'};

v2_root   = fullfile(sim_root, '..', 'Simulations_v2');
real_base = fullfile(sim_root, '..', 'neuronal_decoding', 'results', 'decoding_outputs', ...
                     'procrustes_decoding_basic_incremental_pr_results');

if exist('viz_override', 'var') && isstruct(viz_override)
    fo = fieldnames(viz_override);
    for i = 1:numel(fo), eval(sprintf('%s = viz_override.%s;', fo{i}, fo{i})); end
end
fig_dir = fullfile(sim_root, 'results', 'figures', tag);
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

%% Load
results_file = fullfile(sim_root, 'results', 'decoding_outputs', sprintf('%s_results.mat', tag));
assert(exist(results_file, 'file') == 2, 'No results file %s.', results_file);
L = load(results_file, 'sweep');  sw = L.sweep;
A_list = sw.alpha_list;  N_list = sw.N_list;  S_list = sw.S_list;
nA = numel(A_list);  nN = numel(N_list);  nS = numel(S_list);
R_pop = sw.cfg.R_pop;  K = sw.n_latent;
fprintf('%s: %d / %d points done\n', tag, nnz(sw.done), numel(sw.done));
A_colors = ramp_colors(alpha_hex, nA);
S_colors = ramp_colors(nstim_hex, nS);
A_names  = arrayfun(@(a) sprintf('\\alpha = %.1f', a), A_list, 'UniformOutput', false);
S_names  = arrayfun(@(s) sprintf('n_{stim} = %d', s), S_list, 'UniformOutput', false);
kN = N_list >= min_N_avg;
iS50 = find(S_list == ref_nstim, 1);
iNref = find(N_list == ref_N, 1);

M1 = [];                                            % Model 1 (no shared structure, PR ~ 30)
f1 = fullfile(v2_root, 'results', 'decoding_outputs', 'sweep_N_S_results.mat');
if exist(f1, 'file')
    L1 = load(f1, 'sweep');
    M1 = struct('N', L1.sweep.N_list, 'S', L1.sweep.S_list, 'self', L1.sweep.model1.mean(:,:,col.self), ...
                'ctrl', L1.sweep.model1.mean(:,:,col.ctrl));
end

rec = [];                                           % recorded populations, incremental_pr (50 stimuli)
for p = 1:size(real_pops, 1)
    try
        rec = [rec, load_real(fullfile(real_base, real_pops{p,1}, real_pops{p,2}), ...
                              [real_pops{p,1} ' ' real_pops{p,2}], col)]; %#ok<AGROW>
    catch ME
        warning('Recorded %s %s not loaded (%s).', real_pops{p,1}, real_pops{p,2}, ME.message);
    end
end

RS = [];                                            % recorded stimulus subsampling
f_rs = fullfile(sim_root, 'results', 'decoding_outputs', sprintf('%s_results.mat', real_tag));
if exist(f_rs, 'file')
    Lr = load(f_rs, 'real');  RS = Lr.real;
    fprintf('%s: %d / %d cells done\n', real_tag, nnz(RS.done), numel(RS.done));
end

% model quantities averaged over N >= min_N_avg: [nS x nA]
hr_avg    = squeeze(mean(sw.headroom_of_means(kN, :, :), 1, 'omitnan'));
ratio_avg = squeeze(mean(sw.ratio_of_means(kN, :, :), 1, 'omitnan'));
hr_avg = reshape(hr_avg, nS, nA);  ratio_avg = reshape(ratio_avg, nS, nA);

%% FIGURE 1: heatmaps per alpha
for iA = 1:nA
    self_m = sw.mean(:,:,iA,col.self);  pt_m = sw.mean(:,:,iA,col.pt);  ctrl_m = sw.mean(:,:,iA,col.ctrl);
    ratio_m = sw.ratio_of_means(:,:,iA);  hr_m = sw.headroom_of_means(:,:,iA);  pr_m = sw.PR_emp(:,:,iA);
    figure('Color', 'w', 'Position', [60 60 1500 760]);
    heat_panel(subplot(2,3,1), self_m, N_list, S_list, 'self-decoding', [0 1], seq_ramp([0.10 0.30 0.55], 64), '%.2f');
    heat_panel(subplot(2,3,2), pt_m,   N_list, S_list, 'PT (rotation-only)', [0 1], seq_ramp([0.10 0.30 0.55], 64), '%.2f');
    heat_panel(subplot(2,3,3), ctrl_m, N_list, S_list, 'control (shuffled correspondence)', [0 1], seq_ramp([0.10 0.30 0.55], 64), '%.2f');
    Lr2 = max(abs(log2(ratio_m(:))), [], 'omitnan');  Lr2 = max(Lr2, 1);
    ax = heat_panel(subplot(2,3,4), log2(ratio_m), N_list, S_list, 'PT / control  (1 = no margin)', [-Lr2 Lr2], div_ramp(64), '');
    annotate_cells(ax, ratio_m, '%.2f');
    cb = colorbar(ax);  tk = [-Lr2 0 Lr2];  cb.Ticks = tk;
    cb.TickLabels = arrayfun(@(t) sprintf('%.2g', 2^t), tk, 'UniformOutput', false);
    heat_panel(subplot(2,3,5), hr_m, N_list, S_list, 'headroom (PT - control) / (self - control)', [0 1], seq_ramp([0.62 0.24 0.07], 64), '%.2f');
    heat_panel(subplot(2,3,6), pr_m, N_list, S_list, 'participation ratio (z-scored means)', [0 max(sw.PR_emp(:), [], 'omitnan')], ...
               seq_ramp([0.10 0.30 0.55], 64), '%.1f');
    sgtitle(sprintf('Partial sharing, \\alpha = %.1f  (%d of K = %d latent dimensions shared)   |   %d populations per point', ...
        A_list(iA), sw.n_shared_list(iA), K, R_pop), 'FontWeight', 'bold');
    saveas(gcf, fullfile(fig_dir, sprintf('sweep_heatmaps_alpha%03d.png', round(100 * A_list(iA)))));
end

%% FIGURE 2: curves vs N, columns = n_stim, lines = alpha
figure('Color', 'w', 'Position', [40 40 1700 900]);
tl = tiledlayout(3, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
hA = gobjects(1, nA);
for row = 1:3
    for iS = 1:nS
        ax = nexttile(tl, (row - 1) * nS + iS);  hold(ax, 'on');
        switch row
            case 1
                self_a = mean(sw.mean(:,iS,:,col.self), 3, 'omitnan');
                ctrl_a = mean(sw.mean(:,iS,:,col.ctrl), 3, 'omitnan');
                hs = plot(ax, N_list, self_a, '-', 'Color', 'k', 'LineWidth', 1.5);
                hc = plot(ax, N_list, ctrl_a, '--', 'Color', c_ref, 'LineWidth', 1.5);
                for iA = 1:nA
                    hA(iA) = errorbar(ax, N_list, sw.mean(:,iS,iA,col.pt), sw.sem(:,iS,iA,col.pt), '-o', 'Color', A_colors(iA,:), ...
                        'MarkerFaceColor', A_colors(iA,:), 'MarkerSize', 3, 'LineWidth', 1.2, 'CapSize', 0);
                end
                ylim(ax, [0 1]);
                style_axes(ax, 'N', 'accuracy', sprintf('n_{stim} = %d: PT, self, control', S_list(iS)));
                if iS == 1
                    lg2 = legend(ax, [hs hc hA], [{'self (mean over \alpha)', 'control (mean over \alpha)'}, A_names], ...
                        'Box', 'off', 'FontSize', 8);
                    lg2.Layout.Tile = 'east';
                end
            case 2
                for iA = 1:nA
                    plot(ax, N_list, sw.ratio_of_means(:,iS,iA), '-o', 'Color', A_colors(iA,:), ...
                        'MarkerFaceColor', A_colors(iA,:), 'MarkerSize', 3, 'LineWidth', 1.2);
                end
                yline(ax, 1, 'k-', 'HandleVisibility', 'off');
                set(ax, 'YScale', 'log', 'YTick', [1 1.5 2 3 4 6], 'YTickLabel', {'1', '1.5', '2', '3', '4', '6'});  ylim(ax, [0.8 7]);
                style_axes(ax, 'N', 'PT / control', 'PT relative to the floor');
            case 3
                for iA = 1:nA
                    errorbar(ax, N_list, sw.headroom_of_means(:,iS,iA), sw.headroom_sem(:,iS,iA), '-o', 'Color', A_colors(iA,:), ...
                        'MarkerFaceColor', A_colors(iA,:), 'MarkerSize', 3, 'LineWidth', 1.2, 'CapSize', 0);
                end
                yline(ax, 0, 'k-', 'HandleVisibility', 'off');
                ylim(ax, [-0.25 1.1]);
                style_axes(ax, 'N', 'headroom', 'share of floor-to-self gap recovered');
        end
    end
end
title(tl, sprintf('Partial sharing (K = %d): curves against N for every \\alpha and n_{stim}  |  %d populations per point', K, R_pop), 'FontWeight', 'bold');
saveas(gcf, fullfile(fig_dir, 'sweep_curves_vs_N.png'));

%% FIGURE 3: calibration against alpha
figure('Color', 'w', 'Position', [60 60 1200 820]);
tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
hSl = gobjects(1, nS);
ax = nexttile(tl);  hold(ax, 'on');
plot(ax, [0 1], [0 1], ':', 'Color', c_ref, 'HandleVisibility', 'off');
for iS = 1:nS
    hSl(iS) = plot(ax, A_list, hr_avg(iS,:), '-o', 'Color', S_colors(iS,:), 'MarkerFaceColor', S_colors(iS,:), 'MarkerSize', 5, 'LineWidth', 1.5);
end
yline(ax, 0, 'k-', 'HandleVisibility', 'off');
xlim(ax, [-0.03 1.03]);  ylim(ax, [-0.1 1.05]);
style_axes(ax, '\alpha (shared fraction of latent dimensions)', 'headroom', sprintf('(a) headroom, mean over N \\geq %d (dotted: y = \\alpha)', min_N_avg));
legend(ax, hSl, S_names, 'Location', 'northwest', 'Box', 'off');
ax = nexttile(tl);  hold(ax, 'on');
for iS = 1:nS
    plot(ax, A_list, ratio_avg(iS,:), '-o', 'Color', S_colors(iS,:), 'MarkerFaceColor', S_colors(iS,:), 'MarkerSize', 5, 'LineWidth', 1.5);
end
yline(ax, 1, 'k-', 'HandleVisibility', 'off');
set(ax, 'YScale', 'log', 'YTick', [1 1.5 2 3 4 6], 'YTickLabel', {'1', '1.5', '2', '3', '4', '6'});  xlim(ax, [-0.03 1.03]);  ylim(ax, [0.8 7]);
style_axes(ax, '\alpha', 'PT / control', sprintf('(b) PT / control, mean over N \\geq %d', min_N_avg));
if ~isempty(iNref)
    ax = nexttile(tl);  hold(ax, 'on');
    for iS = 1:nS
        errorbar(ax, A_list, as_row(sw.mean(iNref,iS,:,col.self)), as_row(sw.sem(iNref,iS,:,col.self)), '-o', 'Color', S_colors(iS,:), ...
            'MarkerFaceColor', S_colors(iS,:), 'MarkerSize', 5, 'LineWidth', 1.5, 'CapSize', 0);
    end
    xlim(ax, [-0.03 1.03]);  ylim(ax, [0 1]);
    style_axes(ax, '\alpha', 'self-decoding', sprintf('(c) self-decoding at N = %d (must not depend on \\alpha)', ref_N));
    ax = nexttile(tl);  hold(ax, 'on');
    for iS = 1:nS
        plot(ax, A_list, as_row(sw.PR_emp(iNref,iS,:)), '-o', 'Color', S_colors(iS,:), 'MarkerFaceColor', S_colors(iS,:), 'MarkerSize', 5, 'LineWidth', 1.5);
    end
    xlim(ax, [-0.03 1.03]);  ylim(ax, [0 12]);
    style_axes(ax, '\alpha', 'participation ratio', sprintf('(d) PR at N = %d (must not depend on \\alpha)', ref_N));
end
title(tl, sprintf('Calibration of the partially shared family (K = %d latent dimensions per cue, K_s shared, \\alpha = K_s / K)', K), 'FontWeight', 'bold');
saveas(gcf, fullfile(fig_dir, 'alpha_calibration.png'));

%% FIGURE 4: 50 classes against N, recorded populations overlaid
if ~isempty(iS50)
    figure('Color', 'w', 'Position', [60 60 1600 470]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax1 = nexttile(tl);  hold(ax1, 'on');
    ax2 = nexttile(tl);  hold(ax2, 'on');
    ax3 = nexttile(tl);  hold(ax3, 'on');
    hl = gobjects(0);  hn = {};
    for iA = 1:nA
        c = A_colors(iA,:);
        h = plot(ax1, N_list, sw.headroom_of_means(:,iS50,iA), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
        plot(ax2, N_list, sw.ratio_of_means(:,iS50,iA), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
        plot(ax3, N_list, sw.mean(:,iS50,iA,col.ctrl) ./ sw.mean(:,iS50,iA,col.self), '-o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.5);
        hl(end+1) = h;  hn{end+1} = A_names{iA}; %#ok<SAGROW>
    end
    if ~isempty(M1)
        j50 = find(M1.S == ref_nstim, 1);
        hl(end+1) = plot(ax3, M1.N, M1.ctrl(:,j50) ./ M1.self(:,j50), '-o', 'Color', c_m1, 'MarkerFaceColor', c_m1, 'MarkerSize', 4, 'LineWidth', 1.5);
        hn{end+1} = 'Model 1 (v2, no shared structure)';
    end
    for p = 1:numel(rec)
        c = real_gray(p,:);  st = [real_style{p} real_mk{p}];
        h = plot(ax1, rec(p).N, rec(p).headroom, st, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.2);
        plot(ax2, rec(p).N, rec(p).ratio, st, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.2);
        plot(ax3, rec(p).N, rec(p).cs, st, 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4, 'LineWidth', 1.2);
        hl(end+1) = h;  hn{end+1} = ['recorded ' rec(p).name]; %#ok<SAGROW>
    end
    yline(ax1, 0, 'k-', 'HandleVisibility', 'off');  ylim(ax1, [-0.15 1.05]);
    style_axes(ax1, 'number of neurons N', 'headroom', sprintf('headroom, n_{stim} = %d', ref_nstim));
    yline(ax2, 1, 'k-', 'HandleVisibility', 'off');  set(ax2, 'YScale', 'log', 'YTick', [1 1.5 2 3 4 6], 'YTickLabel', {'1', '1.5', '2', '3', '4', '6'});  ylim(ax2, [0.8 7]);
    style_axes(ax2, 'number of neurons N', 'PT / control', sprintf('PT / control, n_{stim} = %d', ref_nstim));
    ylim(ax3, [0 1.05]);
    style_axes(ax3, 'number of neurons N', 'control / self', sprintf('floor (control / self), n_{stim} = %d', ref_nstim));
    lg = legend(ax3, hl, hn, 'Box', 'off', 'FontSize', 8);  lg.Layout.Tile = 'east';
    title(tl, sprintf('Recorded populations (incremental\\_pr, neuron subsampling) vs the partially shared family at n_{stim} = %d', ref_nstim), ...
        'FontWeight', 'bold');
    saveas(gcf, fullfile(fig_dir, sprintf('real_vs_alpha_nstim%d.png', ref_nstim)));
end

%% FIGURE 5: recorded stimulus subsampling vs alpha
if ~isempty(RS) && any(RS.done(:))
    figure('Color', 'w', 'Position', [60 60 1600 470]);
    tl = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    hl = gobjects(0);  hn = {};
    for m = 1:numel(RS.modes)
        ax = nexttile(tl);  hold(ax, 'on');
        for iA = 1:nA
            h = plot(ax, S_list, hr_avg(:,iA), '-', 'Color', A_colors(iA,:), 'LineWidth', 1.5);
            if m == 1, hl(end+1) = h; hn{end+1} = sprintf('model %s', A_names{iA}); end %#ok<SAGROW>
            text(ax, S_list(end) + 1, hr_avg(end,iA), sprintf('%.1f', A_list(iA)), 'Color', [0.32 0.32 0.30], 'FontSize', 8);
        end
        for p = 1:numel(RS.populations)
            d = squeeze(RS.done(p, m, :)).';
            if ~any(d), continue; end
            c = real_gray(p,:);  st = [real_style{p} real_mk{p}];
            h = errorbar(ax, RS.C_list(d), as_row(RS.headroom_of_means(p, m, d)), as_row(RS.headroom_sem(p, m, d)), st, ...
                'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 5, 'LineWidth', 1.2, 'CapSize', 0);
            if m == 1, hl(end+1) = h; hn{end+1} = ['recorded ' strrep(RS.populations{p}, '_', ' ')]; end %#ok<SAGROW>
        end
        yline(ax, 0, 'k-', 'HandleVisibility', 'off');
        xlim(ax, [5 57]);  ylim(ax, [-0.3 1.05]);  set(ax, 'XTick', S_list);
        if strcmp(RS.modes{m}, 'full_N'), mt = 'all recorded neurons'; else, mt = 'N / C held at the recorded ratio'; end
        style_axes(ax, 'number of stimulus classes C', 'headroom', sprintf('recorded headroom vs C (%s)', mt));
    end
    % floor against neurons per class
    ax = nexttile(tl);  hold(ax, 'on');
    cs_a = mean(sw.mean(:,:,:,col.ctrl), 3, 'omitnan') ./ mean(sw.mean(:,:,:,col.self), 3, 'omitnan');   % [nN x nS], alpha-free
    for iS = 1:nS
        h = plot(ax, N_list / S_list(iS), cs_a(:,iS), '-', 'Color', A_colors(4,:), 'LineWidth', 1);
        if iS == 1, hf1 = h; end
    end
    hf = hf1;  hfn = {sprintf('partially shared family (K = %d, mean over \\alpha)', K)};
    if ~isempty(M1)
        for iS = 1:numel(M1.S)
            h = plot(ax, M1.N / M1.S(iS), M1.ctrl(:,iS) ./ M1.self(:,iS), '-', 'Color', c_m1, 'LineWidth', 1);
            if iS == 1, hf(end+1) = h; hfn{end+1} = 'Model 1 (v2)'; end %#ok<SAGROW>
        end
    end
    for p = 1:numel(RS.populations)
        for m = 1:numel(RS.modes)
            d = squeeze(RS.done(p, m, :)).';
            if ~any(d), continue; end
            x = as_row(RS.N_used(p, m, d)) ./ RS.C_list(d);
            y = as_row(RS.mean(p, m, d, col.ctrl) ./ RS.mean(p, m, d, col.self));
            mfc = real_gray(p,:);  if m == 2, mfc = 'w'; end
            plot(ax, x, y, real_mk{p}, 'Color', real_gray(p,:), 'MarkerFaceColor', mfc, 'MarkerSize', 5, 'LineWidth', 1);
        end
    end
    hf(end+1) = plot(ax, nan, nan, 'd', 'Color', 'k', 'MarkerFaceColor', 'k');  hfn{end+1} = 'recorded, all neurons (filled)';
    hf(end+1) = plot(ax, nan, nan, 'd', 'Color', 'k', 'MarkerFaceColor', 'w');  hfn{end+1} = 'recorded, N/C held (open)';
    set(ax, 'XScale', 'log');  ylim(ax, [0 1.05]);
    xline(ax, 1, 'k:', 'HandleVisibility', 'off');
    style_axes(ax, 'neurons per class N / C', 'control / self', 'floor against neurons per class');
    legend(ax, hf, hfn, 'Location', 'southeast', 'Box', 'off', 'FontSize', 7);
    lg = legend(nexttile(tl, 1), hl, hn, 'Box', 'off', 'FontSize', 7);  lg.Layout.Tile = 'east';
    title(tl, sprintf(['Recorded populations on random subsets of C stimuli (same decoder and control) vs the model headroom ' ...
        '(mean over N \\geq %d)'], min_N_avg), 'FontWeight', 'bold');
    saveas(gcf, fullfile(fig_dir, 'real_stimsub_vs_alpha.png'));
end

%% TEXT SUMMARY
fid = fopen(fullfile(fig_dir, sprintf('summary_%s.txt', tag)), 'w');
for iA = 1:nA
    fprintf(fid, '===== alpha = %.1f (%d of K = %d latent dims shared; %s; %d populations per point; mean, SEM in brackets) =====\n', ...
        A_list(iA), sw.n_shared_list(iA), K, tag, R_pop);
    for iS = 1:nS
        fprintf(fid, '\n-- n_stim = %d (chance %.3f) --\n', S_list(iS), 1/S_list(iS));
        fprintf(fid, '%5s  %13s  %13s  %13s  %13s  %10s  %15s  %8s  %7s  %6s\n', 'N', 'self', 'no-transform', 'PT(rot)', ...
            'control', 'PT/control', 'headroom', 'ctrl-id', 'PR_emp', 'capped');
        for iN = 1:nN
            if ~sw.done(iN,iS,iA), fprintf(fid, '%5d  (not done)\n', N_list(iN)); continue; end
            a = squeeze(sw.mean(iN,iS,iA,:));  e = squeeze(sw.sem(iN,iS,iA,:));
            fprintf(fid, '%5d  %5.3f (%.3f)  %5.3f (%.3f)  %5.3f (%.3f)  %5.3f (%.3f)  %10.2f  %6.2f (%.2f)    %6.3f  %7.1f  %6.2f\n', ...
                N_list(iN), a(col.self), e(col.self), a(col.notf), e(col.notf), a(col.pt), e(col.pt), a(col.ctrl), e(col.ctrl), ...
                sw.ratio_of_means(iN,iS,iA), sw.headroom_of_means(iN,iS,iA), sw.headroom_sem(iN,iS,iA), a(col.ctrl_id), ...
                sw.PR_emp(iN,iS,iA), sw.frac_cv_capped(iN,iS,iA));
        end
    end
    fprintf(fid, '\n');
end
if ~isempty(RS)
    fprintf(fid, '===== recorded populations, random subsets of C stimuli (%s; mean over subsets, SEM in brackets) =====\n', real_tag);
    fprintf(fid, '%-6s %-14s %3s %4s  %13s  %13s  %13s  %10s  %15s  %7s\n', 'pop', 'mode', 'C', 'N', 'self', 'PT(rot)', 'control', ...
        'PT/control', 'headroom', 'PR_emp');
    for p = 1:numel(RS.populations)
        for m = 1:numel(RS.modes)
            for iC = 1:numel(RS.C_list)
                if ~RS.done(p,m,iC), continue; end
                a = squeeze(RS.mean(p,m,iC,:));  e = squeeze(RS.sem(p,m,iC,:));
                fprintf(fid, '%-6s %-14s %3d %4.0f  %5.3f (%.3f)  %5.3f (%.3f)  %5.3f (%.3f)  %10.2f  %6.2f (%.2f)  %7.1f\n', ...
                    RS.populations{p}, RS.modes{m}, RS.C_list(iC), RS.N_used(p,m,iC), a(col.self), e(col.self), a(col.pt), e(col.pt), ...
                    a(col.ctrl), e(col.ctrl), RS.ratio_of_means(p,m,iC), RS.headroom_of_means(p,m,iC), RS.headroom_sem(p,m,iC), RS.PR_emp(p,m,iC));
            end
        end
    end
end
fclose(fid);
fprintf('Figures and summary written to %s\n', fig_dir);


%% ======================================================================
%% Local functions
%% ======================================================================
function ax = heat_panel(ax, M, N_list, S_list, ttl, clim, cmap, fmt)
% Heatmap of a [nN x nS] grid quantity with N on x and n_stim on y, annotated.
im = imagesc(ax, 1:numel(N_list), 1:numel(S_list), M.');
im.AlphaData = double(~isnan(M.'));
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
        text(ax, i, j, sprintf(fmt, M(i,j)), 'HorizontalAlignment', 'center', 'FontSize', 8, 'Color', ink);
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

function C = ramp_colors(hexes, n)
% The validated ordinal ramp when n matches its length; otherwise interpolated along it.
R = cell2mat(cellfun(@(h) [hex2dec(h(2:3)) hex2dec(h(4:5)) hex2dec(h(6:7))] / 255, hexes(:), 'UniformOutput', false));
if n == size(R, 1), C = R; return; end
C = interp1(linspace(0, 1, size(R, 1)), R, linspace(0, 1, n));
end

function v = as_row(A)
% Any array as a row vector.
v = reshape(A, 1, []);
end

function style_axes(ax, xl, yl, ttl)
xlabel(ax, xl); ylabel(ax, yl); title(ax, ttl);
box(ax, 'off'); grid(ax, 'on'); set(ax, 'GridAlpha', 0.15, 'FontSize', 9);
end

function R = load_real(pdir, name, col)
% Pool the six ordered-pair files of one recorded population (incremental_pr
% results): mean accuracies, ratio, headroom and floor at each neuron count.
pairs = {'acec','ecac','ecex','exec','acex','exac'};
acc = containers.Map('KeyType', 'double', 'ValueType', 'any');
for k = 1:numel(pairs)
    s = load(fullfile(pdir, [pairs{k} '_results.mat']));
    fn = fieldnames(s);  res = s.(fn{1});
    for q = 1:numel(res)
        n = res{q}.neuron_num;
        if isKey(acc, n), acc(n) = [acc(n); res{q}.acc]; else, acc(n) = res{q}.acc; end
    end
end
Ns = sort(cell2mat(keys(acc)));
R = struct('name', name, 'N', Ns, 'self', nan(size(Ns)), 'pt', nan(size(Ns)), 'ctrl', nan(size(Ns)), ...
           'ratio', nan(size(Ns)), 'headroom', nan(size(Ns)), 'cs', nan(size(Ns)));
for i = 1:numel(Ns)
    a = mean(acc(Ns(i)), 1);
    R.self(i) = a(col.self);  R.pt(i) = a(col.pt);  R.ctrl(i) = a(col.ctrl);
    R.ratio(i)    = R.pt(i) / R.ctrl(i);
    R.headroom(i) = (R.pt(i) - R.ctrl(i)) / (R.self(i) - R.ctrl(i));
    R.cs(i)       = R.ctrl(i) / R.self(i);
end
end
