function save_figure(fig, S, folder, name)
%{
Filename: save_figure.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Saves a figure as results/figures/<folder>/<name>.png (150 dpi) and .fig,
    then closes it.
%}

d = fullfile(S.fig_root, folder);
if ~exist(d, 'dir'), mkdir(d); end
exportgraphics(fig, fullfile(d, [name '.png']), 'Resolution', 150);
savefig(fig, fullfile(d, [name '.fig']));
close(fig);
fprintf('  saved %s/%s.png\n', folder, name);
end
