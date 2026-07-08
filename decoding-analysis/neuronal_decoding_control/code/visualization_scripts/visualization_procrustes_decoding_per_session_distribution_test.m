%{
Filename: visualization_distribution_test_per_session.m
Author: Zitong Wang
Date: 2026-05-23

Description:
    Visualization of per-cell Mann-Whitney U results produced by
    distribution_test_per_session.m. Three subplots, one per decoding type.

    Within each subplot:
      - x-axis: session columns. V1 sessions on the left, V2 sessions on
                the right, separated by a visual gap.
      - y-axis: centered AUC = P(affine > none) - 0.5
                (positive: affine higher;  negative: affine lower)
      - Each session column carries 3 markers (one per stim pair) with
        a small horizontal jitter:
          o = EC-AC   s = EC-EX   ^ = AC-EX
      - Marker color: V1 = blue, V2 = orange
      - Marker fill:
          filled  -> within-cell Mann-Whitney p < sig_thresh
          unfilled (outline only) -> n.s.
      - Dashed line at y = 0 marks the null (no distributional shift).

    Also prints to the command window:
      - The granular per-(area, session, pair, decoding type) table
        (AUC, centered AUC, raw p, sig marker)
      - A directional summary: fraction of sessions with AUC < 0.5 per
        (area, pair, decoding type).
%}

clc; clear;

%% Configuration
monkey        = 'KO';   % 'FR' or 'KO'
file_location = ['..\..\results\decoding_outputs\','procrustes_decoding_per_session_results'];

areas       = {'V1', 'V2'};
dec_idx     = [1, 2, 6];
dec_labels  = {'self', 'non-PT', 'PT (rot)'};
pair_labels = {'EC-AC', 'EC-EX', 'AC-EX'};
pair_groups = { {'acec','ecac'}, {'ecex','exec'}, {'acex','exac'} };
pair_shapes = {'o', 's', '^'};

sig_thresh  = 0.05;     % within-cell raw p threshold for filled vs unfilled
gap_width   = 1.5;      % visual gap between V1 and V2 session groups
pair_jitter = [-0.22, 0, +0.22];

n_areas = numel(areas);
n_dec   = numel(dec_labels);
n_pairs = numel(pair_labels);

%% Load per-session data and compute AUC + p per cell
pair_names_all = {'acec','ecac','ecex','exec','acex','exac'};
data = struct();
for mode_cell = {'none', 'affine'}
    mode = mode_cell{1};
    for a = 1:n_areas
        area = areas{a};
        for k = 1:numel(pair_names_all)
            pn    = pair_names_all{k};
            fname = fullfile(pwd, file_location, mode, monkey, area, [pn '_results.mat']);
            S     = load(fname);
            data.(mode).(area).(pn) = S.([pn '_results']);
        end
    end
end

results = struct();
for a = 1:n_areas
    area = areas{a};
    n_s  = numel(data.none.(area).acec);
    assert(numel(data.affine.(area).acec) == n_s, ...
        '%s %s: session count mismatch between modes', monkey, area);

    results.(area).n_sessions = n_s;
    results.(area).auc        = nan(n_s, n_pairs, n_dec);
    results.(area).centered   = nan(n_s, n_pairs, n_dec);
    results.(area).p          = nan(n_s, n_pairs, n_dec);

    for s = 1:n_s
        for p = 1:n_pairs
            s1 = pair_groups{p}{1};
            s2 = pair_groups{p}{2};
            none_samples = [data.none.(area).(s1)(s).accuracy; ...
                            data.none.(area).(s2)(s).accuracy];
            aff_samples  = [data.affine.(area).(s1)(s).accuracy; ...
                            data.affine.(area).(s2)(s).accuracy];

            for d_idx = 1:n_dec
                d        = dec_idx(d_idx);
                aff_col  = aff_samples(:, d);
                none_col = none_samples(:, d);

                [pval, ~, stats] = ranksum(aff_col, none_col);
                n_a = numel(aff_col); n_n = numel(none_col);
                U_a = stats.ranksum - n_a * (n_a + 1) / 2;
                auc = U_a / (n_a * n_n);

                results.(area).auc(s, p, d_idx)      = auc;
                results.(area).centered(s, p, d_idx) = auc - 0.5;
                results.(area).p(s, p, d_idx)        = pval;
            end
        end
    end
end

%% Print granular per-(area, session, pair, decoding) table
fprintf('\n=== Mann-Whitney U per (session, pair, decoding type) | monkey %s ===\n', monkey);
fprintf('AUC = P(affine > none) + 0.5*P(affine == none);  centered = AUC - 0.5\n');
fprintf('Two-sided p from ranksum. No multiple-comparisons correction at this level.\n\n');
fprintf('Area | Session | Pair  | Decoding   |  AUC   | centered |   p      | sig\n');
fprintf('%s\n', repmat('-', 1, 78));

for a = 1:n_areas
    area = areas{a};
    n_s  = numel(data.none.(area).acec);
    for s = 1:n_s
        for p = 1:n_pairs
            for d_idx = 1:n_dec
                fprintf('%-4s |   %3d   | %-5s | %-9s  | %5.3f  | %+6.3f   | %.4f   | %s\n', ...
                    area, s, pair_labels{p}, dec_labels{d_idx}, ...
                    results.(area).auc(s, p, d_idx), ...
                    results.(area).centered(s, p, d_idx), ...
                    results.(area).p(s, p, d_idx), ...
                    sig_marker(results.(area).p(s, p, d_idx)));
            end
        end
        fprintf('%s\n', repmat('-', 1, 78));
    end
end

%% Directional summary
fprintf('\n=== Directional summary: fraction of sessions with AUC < 0.5 (affine < none) ===\n');
fprintf('Area | Pair  | Decoding   | n_below / n_total\n');
fprintf('%s\n', repmat('-', 1, 44));
for a = 1:n_areas
    area = areas{a};
    for p = 1:n_pairs
        for d_idx = 1:n_dec
            v       = results.(area).auc(:, p, d_idx);
            n_below = sum(v < 0.5);
            n_total = numel(v);
            fprintf('%-4s | %-5s | %-9s  |   %d / %d\n', ...
                area, pair_labels{p}, dec_labels{d_idx}, n_below, n_total);
        end
    end
end
fprintf('\n');

%% Layout: x-positions for sessions
n_v1    = results.V1.n_sessions;
n_v2    = results.V2.n_sessions;
v1_x    = 1:n_v1;
v2_x    = (n_v1 + gap_width + 1) : (n_v1 + gap_width + n_v2);
sep_x   = n_v1 + gap_width/2 + 0.5;
x_left  = 0.3;
x_right = v2_x(end) + 0.7;

%% Symmetric y-range across the three subplots
all_centered = [results.V1.centered(:); results.V2.centered(:)];
all_centered = all_centered(~isnan(all_centered) & isfinite(all_centered));
y_max        = max(0.5, max(abs(all_centered)) * 1.15);

%% Plot
v1_base = [0      0.4470 0.7410];
v2_base = [0.8500 0.3250 0.0980];

figure('Position', [100 100 1500 500]);
ax_handles = gobjects(1, n_dec);

for d_idx = 1:n_dec
    ax_handles(d_idx) = subplot(1, n_dec, d_idx);
    hold on;

    % V1 sessions
    for s = 1:n_v1
        x_center = v1_x(s);
        for p = 1:n_pairs
            x_pos = x_center + pair_jitter(p);
            y_pos = results.V1.centered(s, p, d_idx);
            p_val = results.V1.p(s, p, d_idx);
            plot_cell(x_pos, y_pos, pair_shapes{p}, v1_base, p_val < sig_thresh);
        end
    end

    % V2 sessions
    for s = 1:n_v2
        x_center = v2_x(s);
        for p = 1:n_pairs
            x_pos = x_center + pair_jitter(p);
            y_pos = results.V2.centered(s, p, d_idx);
            p_val = results.V2.p(s, p, d_idx);
            plot_cell(x_pos, y_pos, pair_shapes{p}, v2_base, p_val < sig_thresh);
        end
    end

    % Zero line + V1/V2 separator
    yline(0,     '--k', 'LineWidth', 1, 'HandleVisibility','off');
    xline(sep_x, ':',   'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'HandleVisibility','off');

    xlim([x_left, x_right]);
    ylim([-y_max, y_max]);

    % x-tick labels: S1..S5 then S1..S7
    xticks([v1_x, v2_x]);
    xticklabels(cellfun(@(i) sprintf('S%d', i), num2cell([1:n_v1, 1:n_v2]), ...
        'UniformOutput', false));
    xtickangle(0);

    title(dec_labels{d_idx});
    if d_idx == 1
        ylabel('Centered AUC: P(affine > none) - 0.5');
    end

    % Group labels at the top
    text(mean(v1_x), y_max * 0.95, 'V1', ...
        'HorizontalAlignment','center', 'FontWeight','bold', ...
        'FontSize', 12, 'Color', v1_base);
    text(mean(v2_x), y_max * 0.95, 'V2', ...
        'HorizontalAlignment','center', 'FontWeight','bold', ...
        'FontSize', 12, 'Color', v2_base);

    grid on; box on;

    % Legend in the rightmost subplot
    if d_idx == n_dec
        hl_v1 = plot(NaN, NaN, 'o', 'MarkerFaceColor', v1_base, ...
            'MarkerEdgeColor', v1_base, 'MarkerSize', 8);
        hl_v2 = plot(NaN, NaN, 'o', 'MarkerFaceColor', v2_base, ...
            'MarkerEdgeColor', v2_base, 'MarkerSize', 8);
        hl_p1 = plot(NaN, NaN, pair_shapes{1}, 'MarkerFaceColor', [0.4 0.4 0.4], ...
            'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerSize', 8);
        hl_p2 = plot(NaN, NaN, pair_shapes{2}, 'MarkerFaceColor', [0.4 0.4 0.4], ...
            'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerSize', 8);
        hl_p3 = plot(NaN, NaN, pair_shapes{3}, 'MarkerFaceColor', [0.4 0.4 0.4], ...
            'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerSize', 8);
        hl_sig  = plot(NaN, NaN, 'o', 'MarkerFaceColor', [0.4 0.4 0.4], ...
            'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerSize', 8);
        hl_nsig = plot(NaN, NaN, 'o', 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', [0.4 0.4 0.4], 'MarkerSize', 8, 'LineWidth', 1.4);

        legend([hl_v1, hl_v2, hl_p1, hl_p2, hl_p3, hl_sig, hl_nsig], ...
            {'V1', 'V2', pair_labels{1}, pair_labels{2}, pair_labels{3}, ...
             sprintf('p < %.2f', sig_thresh), 'n.s.'}, ...
            'Location', 'eastoutside', 'NumColumns', 1);
    end

    hold off;
end

linkaxes(ax_handles, 'y');
sgtitle(sprintf(['%s  |  per-cell Mann-Whitney U (affine vs none)  |  ' ...
                 'shape = stim pair, color = area, filled = p < %.2f'], ...
                monkey, sig_thresh));

%% =================== LOCAL FUNCTIONS ===================
function plot_cell(x, y, shape, area_color, is_sig)
% Draw one marker. Filled if significant, outline-only if not.
    if isnan(y), return; end
    if is_sig
        plot(x, y, shape, ...
            'MarkerFaceColor', area_color, ...
            'MarkerEdgeColor', area_color, ...
            'MarkerSize', 8, 'LineWidth', 1.0, ...
            'HandleVisibility','off');
    else
        plot(x, y, shape, ...
            'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', area_color, ...
            'MarkerSize', 8, 'LineWidth', 1.4, ...
            'HandleVisibility','off');
    end
end

function s = sig_marker(p)
    if isnan(p),       s = '';
    elseif p < 0.001,  s = '***';
    elseif p < 0.01,   s = '**';
    elseif p < 0.05,   s = '*';
    else,              s = 'n.s.';
    end
end