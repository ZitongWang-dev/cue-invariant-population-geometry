function b = plot_grouped_bars(ax, M, E, group_labels, series, S, y_max)
%{
Filename: plot_grouped_bars.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Grouped bar chart in the style of the published figures: one group per
    stimulus pair / rendering, one bar per series, red error bars, dashed 2%
    chance line. Each series = struct with fields
        name, color, edge ('self' | 'nopt' | 'pt' | 'ctrl'), alpha (optional)

Inputs:
    ax            target axes
    M, E          [n_groups x n_series] means and error half-widths
    group_labels  {1 x n_groups}
    series        1 x n_series struct array
    S             figure_settings()
    y_max         upper y limit
%}

b = bar(ax, M, 'grouped');
for i = 1:numel(series)
    st = S.edge.(series(i).edge);
    b(i).FaceColor = series(i).color;
    b(i).EdgeColor = [0 0 0];
    b(i).LineStyle = st{1};
    b(i).LineWidth = st{2};
    if isfield(series, 'alpha') && ~isempty(series(i).alpha), b(i).FaceAlpha = series(i).alpha; end
end
hold(ax, 'on');
for i = 1:numel(series)
    errorbar(ax, b(i).XEndPoints, M(:, i), E(:, i), 'Color', S.c.err, 'LineStyle', 'none', ...
        'LineWidth', 1, 'HandleVisibility', 'off');
end
yline(ax, S.chance, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1);
hold(ax, 'off');
set(ax, 'XTick', 1:numel(group_labels), 'XTickLabel', group_labels, 'FontSize', 10);
ylim(ax, [0 y_max]);
ylabel(ax, 'Decoding accuracy');
box(ax, 'on');
end
