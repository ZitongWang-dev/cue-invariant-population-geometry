%{
% FILENAME: calibrate_from_real_data.m   (Simulations_v2)
% AUTHOR:   Zitong Wang
% DATE:     2026-09-17
%
% DESCRIPTION:
%   Measures the single-neuron statistics of the recorded populations that the
%   v2 generators reproduce, so the synthetic operating point is set from the
%   data rather than by hand. For every neuron of every population (FR/KO x
%   V1/V2) and every cue, using the same loader and window as the decoding
%   pipeline (multiclass_svmloader_PT, 330-630 ms, 50 stimuli x 10 trials):
%
%     mean count    r    = mean over stimuli of the trial-averaged count
%     tuning depth  cv   = sqrt(max(var_s(mean_s) - mean_s(var_t)/10, 0)) / r
%                          = SD of the tuning curve across the 50 stimuli, with
%                          the trial-noise contribution to each 10-trial mean
%                          removed, divided by the mean count (a noise-corrected
%                          coefficient of variation of the tuning curve)
%     Fano factor        = mean over stimuli of (trial variance / trial mean)
%
%   and fits a log-normal (median, log-SD; robust: median and IQR/1.349 of the
%   log values) to r and to cv across neurons. The pooled fits are the defaults
%   of generate_model1_trial_data.m / generate_model3_trial_data.m
%   (rate_median, rate_log_sd, cv_median, cv_log_sd). A Fano factor ~ 1 is the
%   justification for drawing Poisson trials.
%
%   Why this replaces the v1 calibration (Simulations/.../calibrate_weight_scale.m):
%   v1 set one signal SD ("weight_scale") and a large positivity offset for all
%   neurons, then searched weight_scale for a target self-decoding, so the
%   tuning depth was a by-product of the offset. v2 sets depth and rate per
%   neuron from these measured distributions and lets self-decoding be whatever
%   the geometry (Model 1 vs Model 3), N and n_stim make it.
%
% INPUTS:
%   ../../../neuronal_decoding_control/neuronal_data/<M>/<A>/<M>_<A>_allstim.mat
%   ../../../neuronal_decoding_control/code/decoding_scripts/multiclass_svmloader_PT.m
%
% OUTPUT:
%   results/calibration/empirical_operating_point.mat  (per-neuron values + fits)
%   results/figures/calibration_single_neuron_stats.png
%}

%% INITIALIZATION
clc; clear;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
sim_root  = fullfile(script_dir, '..', '..');                       % Simulations_v2/
da_root   = fullfile(sim_root, '..');                               % decoding-analysis/
addpath(fullfile(da_root, 'neuronal_decoding_control', 'code', 'decoding_scripts'));
data_base = fullfile(da_root, 'neuronal_decoding_control', 'neuronal_data');

pops       = {'FR','V1'; 'FR','V2'; 'KO','V1'; 'KO','V2'};
cue_names  = {'ac','ec','ex'};
timewindow = [330 630];
n_stim     = 50;
n_trials   = 10;

%% MEASURE
rows = [];   % [pop cue neuron r cv fano cv_raw maxdev]
fprintf('%-6s %-3s %4s | %-24s | %-26s | %-16s\n', 'pop', 'cue', 'N', 'count: med (IQR)', 'tuning cv: med (IQR)', 'Fano: med (IQR)');
for p = 1:size(pops,1)
    m = pops{p,1}; a = pops{p,2};
    d = load(fullfile(data_base, m, a, sprintf('%s_%s_allstim.mat', m, a)), 'three_stim_array');
    for c = 1:numel(d.three_stim_array)
        data = multiclass_svmloader_PT(d.three_stim_array{c}, timewindow);   % [500 x N]
        N = size(data, 2);
        means = zeros(n_stim, N); vars = zeros(n_stim, N);
        for s = 1:n_stim
            r_idx = (s-1)*n_trials + (1:n_trials);
            means(s,:) = mean(data(r_idx,:), 1);
            vars(s,:)  = var(data(r_idx,:), 0, 1);
        end
        r       = mean(means, 1);
        noise_v = mean(vars, 1);
        sig_v   = max(var(means, 0, 1) - noise_v / n_trials, 0);   % noise-corrected tuning variance
        cv      = sqrt(sig_v) ./ max(r, eps);
        cv_raw  = std(means, 0, 1) ./ max(r, eps);                 % uncorrected (for reference)
        fano    = noise_v ./ max(r, eps);
        maxdev  = max(abs(means - r), [], 1) ./ max(r, eps);
        rows = [rows; repmat([p c], N, 1), (1:N).', r.', cv.', fano.', cv_raw.', maxdev.']; %#ok<AGROW>
        fprintf('%-6s %-3s %4d | %5.2f (%5.2f-%5.2f)      | %5.3f (%5.3f-%5.3f)      | %4.2f (%4.2f-%4.2f)\n', ...
            [m ' ' a], cue_names{c}, N, median(r), prctile(r,25), prctile(r,75), ...
            median(cv), prctile(cv,25), prctile(cv,75), median(fano), prctile(fano,25), prctile(fano,75));
    end
end

%% POOLED LOG-NORMAL FITS (robust)
r_all    = rows(:,4);
cv_all   = rows(:,5);
fano_all = rows(:,6);
ok_r  = r_all > 0;
ok_cv = cv_all > 0;                        % entries whose corrected tuning variance floored at 0 are excluded
fit = struct();
fit.rate_median  = median(r_all(ok_r));
fit.rate_log_sd  = iqr(log(r_all(ok_r))) / 1.349;
fit.cv_median    = median(cv_all(ok_cv));
fit.cv_log_sd    = iqr(log(cv_all(ok_cv))) / 1.349;
fit.frac_cv_zero = mean(~ok_cv);
fit.fano_median  = median(fano_all);
fit.cv_pctile    = prctile(cv_all(ok_cv), [5 25 50 75 95]);

fprintf('\n=== pooled over %d neuron-cue entries (%d neurons x 3 cues) ===\n', size(rows,1), round(size(rows,1)/3));
fprintf('mean count  : median %.2f, log-SD %.2f  (5-95%%: %.1f-%.1f)\n', fit.rate_median, fit.rate_log_sd, prctile(r_all,5), prctile(r_all,95));
fprintf('tuning cv   : median %.3f, log-SD %.2f  (5/25/50/75/95%%: %s); %.1f%% of entries floored at 0 (excluded)\n', ...
    fit.cv_median, fit.cv_log_sd, num2str(fit.cv_pctile, '%.3f '), 100*fit.frac_cv_zero);
fprintf('Fano factor : median %.2f (IQR %.2f-%.2f)  -> Poisson trials are adequate\n', fit.fano_median, prctile(fano_all,25), prctile(fano_all,75));
fprintf('\nSuggested generator defaults:\n  rate_median = %.0f;  rate_log_sd = %.2f;  cv_median = %.2f;  cv_log_sd = %.2f\n', ...
    fit.rate_median, fit.rate_log_sd, fit.cv_median, fit.cv_log_sd);

%% SAVE
out_dir = fullfile(sim_root, 'results', 'calibration');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
columns = {'pop_idx','cue_idx','neuron','mean_count','tuning_cv','fano','tuning_cv_uncorrected','max_abs_dev_over_mean'};
save(fullfile(out_dir, 'empirical_operating_point.mat'), 'rows', 'columns', 'pops', 'cue_names', 'fit', 'timewindow');

%% FIGURE: measured distributions vs the log-normal fits used by the generators
fig_dir = fullfile(sim_root, 'results', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
figure('Color', 'w', 'Position', [100 100 1200 340]);

subplot(1,3,1); hold on;
edges = linspace(-1, 2.2, 33);
histogram(log10(r_all(ok_r)), edges, 'Normalization', 'pdf', 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
xx = linspace(edges(1), edges(end), 200);
plot(xx, normpdf(xx, log10(fit.rate_median), fit.rate_log_sd / log(10)), 'k-', 'LineWidth', 1.5);
xlabel('log_{10} mean count (330-630 ms)'); ylabel('density'); box off;
title(sprintf('mean count: median %.1f, log-SD %.2f', fit.rate_median, fit.rate_log_sd));

subplot(1,3,2); hold on;
edges = linspace(-2.2, 0.3, 33);
histogram(log10(cv_all(ok_cv)), edges, 'Normalization', 'pdf', 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
xx = linspace(edges(1), edges(end), 200);
plot(xx, normpdf(xx, log10(fit.cv_median), fit.cv_log_sd / log(10)), 'k-', 'LineWidth', 1.5);
xlabel('log_{10} tuning depth (noise-corrected CV)'); ylabel('density'); box off;
title(sprintf('tuning depth: median %.2f, log-SD %.2f', fit.cv_median, fit.cv_log_sd));

subplot(1,3,3); hold on;
histogram(fano_all, linspace(0, 3, 31), 'Normalization', 'pdf', 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none');
xline(1, 'k--', 'Poisson');
xlabel('Fano factor (trial var / mean)'); ylabel('density'); box off;
title(sprintf('Fano: median %.2f', fit.fano_median));

sgtitle('Single-neuron statistics of the recorded populations (all 4 populations x 3 cues)');
saveas(gcf, fullfile(fig_dir, 'calibration_single_neuron_stats.png'));
fprintf('Saved %s\n', fullfile(fig_dir, 'calibration_single_neuron_stats.png'));
