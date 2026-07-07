%{
% FILENAME: calibrate_weight_scale.m
% AUTHOR:   Zitong Wang
% DATE:     2026-07-01
%
% DESCRIPTION:
%   Calibration sweep for a simulation-family benchmark. The model is selected
%   by cfg.model_fn (contract: [data_trial, labels, diag] = model_fn(gen_cfg)),
%   so the same script calibrates Models 1/2/3/5/6. Finds the signal magnitude
%   that puts self-decoding in a target band, and exposes the auto-vs-fixed
%   baseline tradeoff (and the modulation-depth ceiling) in one table.
%
%   NOTE: the useful ws_grid range is model-specific. Model 3 needs LARGE values
%   (self climbs slowly, capped by rank-2 + positivity); Model 1 self-decodes
%   easily and needs SMALL values to bring self DOWN into the band. Set ws_grid
%   accordingly per model.
%
%   Signal magnitude is controlled by weight_scale. For Model 3 this is the SD
%   of the readout weights; rescaling the stimulus grid is algebraically
%   identical. For Model 1 it is the SD of the random means directly. Seeds are
%   averaged per point so the trend is not buried in single-seed noise.
%
%   Two baseline modes are swept side by side:
%     'auto'   baseline = ceil(5.5 * weight_scale * M_grid), tracks the signal.
%              Positivity guaranteed. Signal and Poisson noise grow together, so
%              SNR ~ sqrt(weight_scale) and the mean count climbs with the sweep.
%     'fixed'  baseline held constant (= cfg.fixed_baseline). Over its VALID
%              range (weight_scale <~ baseline/7.6) noise is fixed and signal
%              grows, so SNR ~ weight_scale there -- but past that edge rates go
%              negative, frac_floored leaves zero, and the model is rectified
%              (no longer clean rank-2 linear). The usable point is the largest
%              self with frac_floored == 0.
%
%   Expected reading: neither mode reaches a data-like self band (0.4-0.6) at a
%   realistic mean count. 50 classes packed in a 2D plane, with modulation depth
%   capped near 13-19% by positivity, cannot be resolved without either absurd
%   counts (auto) or rectification (fixed pushed past its edge). That ceiling is
%   the informative result.
%
%   Reported per (mode, weight_scale), averaged over seeds:
%     self         self-decoding accuracy (ECOC on cue 1, k-fold CV), mean +/- SD
%     mean_count   mean spike count
%     frac_floored max fraction of rates floored across seeds (want exactly 0)
%
% USAGE:
%   Run the script. Edit the CONFIGURATION block to change the sweep.
%}

%% INITIALIZATION
clc; clear;

%% CONFIGURATION
cfg = struct();

% --- model selection (swap these two to calibrate a different model) ---
cfg.model_name     = 'model1';
cfg.model_fn       = @generate_model1_trial_data;   % generator handle, contract: model_fn(gen_cfg)

% ---- locked operating points (for reference / reruns) ----
%   Model 3: model_fn=@generate_model3_trial_data, baseline='auto', rho=0, N=100
%            weight_scale=12 -> self~0.49, PR~1.98, count~138, frac_floored=0
%   Model 1: model_fn=@generate_model1_trial_data, baseline='auto', rho=0, N=100
%            matched-self:  weight_scale=0.70 -> self~0.48, PR~33, count~4,  frac_floored=0
%            near-ceiling:  weight_scale=2    -> self~0.97, PR~33, count~11, frac_floored=0

% --- sweep: ws_grid is model-specific ---
%   Model 3: [0.7 1 1.5 2 3 4 6 8 10 12 14]   self climbs slowly (rank-2 + positivity cap); in band ~12
%   Model 1: [0.2 0.3 0.5 0.7 1 1.5 2 3]      self near ceiling; calibrate DOWN into the band
cfg.N              = 100;
cfg.ws_grid        = [0.2 0.3 0.5 0.7 1 1.5 2 3];   % Model 1 range
cfg.n_seeds        = 6;                        % seeds averaged per point
cfg.base_seed      = 2000;
cfg.fixed_baseline = 8;                        % baseline for the 'fixed' mode (data-like count)
cfg.n_folds        = 5;                        % k-fold for the self-decoding proxy
cfg.target_self    = [0.4 0.6];                % band to flag
cfg.make_figure    = true;
cfg.save           = true;

% Baseline modes: {label, baseline value passed to the generator}
modes      = {'auto', cfg.fixed_baseline};
mode_names = {'auto', sprintf('fixed(%g)', cfg.fixed_baseline)};

%% BUILD JOB LIST (mode x weight_scale x seed) AND RUN
n_mode = numel(modes);
n_ws   = numel(cfg.ws_grid);
n_seed = cfg.n_seeds;

[MI, WI, SI] = ndgrid(1:n_mode, 1:n_ws, 1:n_seed);
MI = MI(:); WI = WI(:); SI = SI(:);
n_job = numel(MI);

self_j  = zeros(n_job, 1);
mc_j    = zeros(n_job, 1);
floor_j = zeros(n_job, 1);
pr_j    = zeros(n_job, 1);

ws_grid  = cfg.ws_grid;          % local copies for clean parfor slicing
base_seed = cfg.base_seed;
N = cfg.N;  n_folds = cfg.n_folds;
model_fn = cfg.model_fn;         % local handle for parfor broadcast

parfor j = 1:n_job
    gc = struct('N', N, 'weight_scale', ws_grid(WI(j)), 'rho', 0, ...
                'baseline', modes{MI(j)}, 'seed', base_seed + SI(j));

    % Suppress the generator's floor warning; frac_floored is read as a column.
    warn_state = warning('off', 'all');
    [d, lab, g] = model_fn(gc);
    warning(warn_state);

    self_j(j)  = self_decode(d{1}, lab, n_folds);
    mc_j(j)    = g.mean_count;
    floor_j(j) = g.frac_floored;
    pr_j(j)    = mean(g.signal_PR);   % mean over cues; PR is a model-class property
end

%% AGGREGATE OVER SEEDS
self_arr  = reshape(self_j,  [n_mode, n_ws, n_seed]);
mc_arr    = reshape(mc_j,    [n_mode, n_ws, n_seed]);
floor_arr = reshape(floor_j, [n_mode, n_ws, n_seed]);
pr_arr    = reshape(pr_j,    [n_mode, n_ws, n_seed]);

self_mean = mean(self_arr, 3);        % [n_mode x n_ws]
self_sd   = std(self_arr, 0, 3);
mc_mean   = mean(mc_arr, 3);
floor_max = max(floor_arr, [], 3);    % worst seed (want 0)
pr_mean   = mean(pr_arr, 3);          % ~2 when clean; inflates once rectified

%% PRINT TABLES
for mi = 1:n_mode
    fprintf('\n=== baseline mode: %s ===\n', mode_names{mi});
    fprintf('%8s  %14s  %11s  %8s  %12s\n', 'weight', 'self (mean±sd)', 'mean_count', 'PR', 'frac_floor');
    for wi = 1:n_ws
        flag = '';
        if floor_max(mi,wi) == 0 && ...
           self_mean(mi,wi) >= cfg.target_self(1) && self_mean(mi,wi) <= cfg.target_self(2)
            flag = '  <-- in band';
        end
        if floor_max(mi,wi) > 0
            flag = [flag, '  [FLOORED: rectifying]'];
        end
        fprintf('%8.2f  %6.3f ± %.3f  %11.1f  %8.3f  %12.4f%s\n', ...
            cfg.ws_grid(wi), self_mean(mi,wi), self_sd(mi,wi), ...
            mc_mean(mi,wi), pr_mean(mi,wi), floor_max(mi,wi), flag);
    end
end

%% RECOMMENDATION PER MODE (best clean, i.e. unfloored, operating point)
fprintf('\n=== best clean operating point (frac_floored == 0) ===\n');
mid = mean(cfg.target_self);
for mi = 1:n_mode
    clean = find(floor_max(mi,:) == 0);
    if isempty(clean)
        fprintf('%-10s : no unfloored point in the swept range.\n', mode_names{mi});
        continue;
    end
    [~, ix] = min(abs(self_mean(mi,clean) - mid));
    wi = clean(ix);
    in_band = self_mean(mi,wi) >= cfg.target_self(1) && self_mean(mi,wi) <= cfg.target_self(2);
    status  = 'IN BAND'; if ~in_band, status = 'below band (ceiling)'; end
    fprintf('%-10s : weight_scale=%.2f -> self=%.3f, count=%.1f  [%s]\n', ...
        mode_names{mi}, cfg.ws_grid(wi), self_mean(mi,wi), mc_mean(mi,wi), status);
end
fprintf('(chance = %.3f)\n', 1/50);

%% SAVE
if cfg.save
    save_path = fullfile('..','..','results','decoding_outputs', cfg.model_name, 'calibration');
    if ~exist(save_path, 'dir'), mkdir(save_path); end
    sweep = struct('cfg', cfg, 'mode_names', {mode_names}, 'ws_grid', cfg.ws_grid, ...
                   'self_mean', self_mean, 'self_sd', self_sd, ...
                   'mean_count', mc_mean, 'signal_PR', pr_mean, 'frac_floored', floor_max);
    save(fullfile(save_path, 'calibration_sweep.mat'), 'sweep');
    fprintf('\nSaved sweep to %s\n', save_path);
end

%% FIGURE
if cfg.make_figure
    colors = {[0 0 0], [0.85 0.2 0.2]};   % auto = black, fixed = red
    figure('Color', 'w', 'Position', [100 100 1300 380]);

    % Panel 1: self-decoding vs signal magnitude
    subplot(1,3,1); hold on;
    for mi = 1:n_mode
        errorbar(cfg.ws_grid, self_mean(mi,:), self_sd(mi,:), '-o', ...
            'Color', colors{mi}, 'MarkerFaceColor', colors{mi}, 'DisplayName', mode_names{mi});
        fl = floor_max(mi,:) > 0;                      % mark rectified points
        if any(fl)
            plot(cfg.ws_grid(fl), self_mean(mi,fl), 'x', 'Color', colors{mi}, ...
                'MarkerSize', 12, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    yline(cfg.target_self(1), '--', 'target', 'HandleVisibility', 'off');
    yline(cfg.target_self(2), '--', '', 'HandleVisibility', 'off');
    yline(1/50, ':', 'chance', 'HandleVisibility', 'off');
    xlabel('weight\_scale (signal magnitude)'); ylabel('self-decoding');
    legend('Location', 'northwest'); box off;
    title('self-decoding ( x = rectified )');

    % Panel 2: mean count vs signal magnitude
    subplot(1,3,2); hold on;
    for mi = 1:n_mode
        plot(cfg.ws_grid, mc_mean(mi,:), '-o', 'Color', colors{mi}, ...
            'MarkerFaceColor', colors{mi}, 'DisplayName', mode_names{mi});
    end
    xlabel('weight\_scale (signal magnitude)'); ylabel('mean spike count');
    legend('Location', 'northwest'); box off;
    title('mean count (the auto tradeoff)');

    % Panel 3: signal PR vs signal magnitude (rank-2 => 2; inflates once rectified)
    subplot(1,3,3); hold on;
    for mi = 1:n_mode
        plot(cfg.ws_grid, pr_mean(mi,:), '-o', 'Color', colors{mi}, ...
            'MarkerFaceColor', colors{mi}, 'DisplayName', mode_names{mi});
        fl = floor_max(mi,:) > 0;
        if any(fl)
            plot(cfg.ws_grid(fl), pr_mean(mi,fl), 'x', 'Color', colors{mi}, ...
                'MarkerSize', 12, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    yline(2, ':', 'rank-2', 'HandleVisibility', 'off');
    xlabel('weight\_scale (signal magnitude)'); ylabel('signal PR');
    legend('Location', 'best'); box off;
    title('participation ratio ( x = rectified )');

    sgtitle(sprintf('%s calibration: self-decoding, count, and PR vs signal magnitude', cfg.model_name));
    saveas(gcf, fullfile('..','..','results','decoding_outputs', cfg.model_name, ...
        'calibration', 'calibration_sweep.png'));   % uncomment for production
end


%% ======================================================================
%% Local functions
%% ======================================================================
function acc = self_decode(counts, labels, n_folds)
% Self-decoding proxy: linear ECOC on one cue's z-scored trials, k-fold CV.
model = fitcecoc(zscore(counts), labels);
acc   = 1 - kfoldLoss(crossval(model, 'KFold', n_folds));
end
