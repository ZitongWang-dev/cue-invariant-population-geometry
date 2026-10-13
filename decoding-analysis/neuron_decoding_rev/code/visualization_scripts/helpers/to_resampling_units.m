function U = to_resampling_units(A, unit_dim)
%{
Filename: to_resampling_units.m   (visualization helper)
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Reshapes an accuracy array to [n_units x everything else], with the
    resampling unit (neuron draw, or fold when there is a single draw) as rows.
    Rows that are entirely NaN (unused draw slots) are dropped. Feeds
    bootstrap_metrics.m.
%}

nd = max(ndims(A), unit_dim);
A  = permute(A, [unit_dim, setdiff(1:nd, unit_dim)]);
U  = reshape(A, size(A, 1), []);
U  = U(~all(isnan(U), 2), :);
end
