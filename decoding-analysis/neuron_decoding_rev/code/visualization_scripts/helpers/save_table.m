function save_table(T, S, folder, name)
%{
Filename: save_table.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Writes the plotted numbers of a figure series as
    results/figures/<folder>/<name>.csv, next to its figures.
%}

d = fullfile(S.fig_root, folder);
if ~exist(d, 'dir'), mkdir(d); end
writetable(T, fullfile(d, [name '.csv']));
fprintf('  wrote %s/%s.csv (%d rows)\n', folder, name, height(T));
end
