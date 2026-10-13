function [m, e] = mean_sem(x)
%{
Filename: mean_sem.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Mean and SEM over all non-NaN entries of an array of any shape: the error
    bars of the published figures (SEM over all repeat estimates).
%}

x = x(~isnan(x));
m = mean(x);
e = std(x) / sqrt(numel(x));
end
