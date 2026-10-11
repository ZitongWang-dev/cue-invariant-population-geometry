function run_all(run_set)
%{
Filename: run_all.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Master script: runs the decoding drivers one after another, unattended, and
    writes everything printed to a log file in results/logs/.

        run_all('all')       every driver, main-figure analyses first and the long
                             PC halves sweep last (about 20 h on one 8-worker
                             machine at the 2026-10-11 speed)
        run_all('desktop')   one half of a two-machine split (about 11 h)
        run_all('laptop')    the other half (about 9 h)
        run_all('test')      tiny configurations (about 8 min) written to a
                             temporary folder; run it once on each machine first
                             to check MATLAB, the toolboxes and the parallel pool

    run_all with no argument uses DEFAULT_SET below, so pressing Run in the
    editor also works.

    Each job runs in its own workspace (run_driver.m). A job that stops with an
    error is logged and skipped, and the next job starts. Every driver resumes
    where it stopped, so after fixing a problem just call run_all again with the
    same set.
%}

DEFAULT_SET = 'all';
if nargin < 1 || isempty(run_set), run_set = DEFAULT_SET; end

code_dir = fileparts(mfilename('fullpath'));
addpath(code_dir);

%% Jobs: {driver, cfg_override ([] = defaults), description, hours at the measured speed}
switch lower(run_set)
    case 'all'
        jobs = {
            'run_cue_transfer_incremental', [],                                  'cue transfer vs N (Fig 3b, 4c)',            8.3
            'run_across_area',              [],                                  'across areas (Fig 4b, 12)',                 1.1
            'run_halves_incremental',       [],                                  'within-area halves vs size (Fig 4a, 13)',   2.7
            'run_pc_sweep', struct('comparisons', {{'cue', 'across'}}),          'PC sweep: cue transfer and across areas',   1.7
            'run_pc_sweep', struct('comparisons', {{'halves'}}),                 'PC sweep: within-area halves',              6.2};
    case 'desktop'
        jobs = {
            'run_cue_transfer_incremental', [],                                  'cue transfer vs N (Fig 3b, 4c)',            8.3
            'run_across_area',              [],                                  'across areas (Fig 4b, 12)',                 1.1
            'run_pc_sweep', struct('comparisons', {{'cue', 'across'}}),          'PC sweep: cue transfer and across areas',   1.7};
    case 'laptop'
        jobs = {
            'run_halves_incremental',       [],                                  'within-area halves vs size (Fig 4a, 13)',   2.7
            'run_pc_sweep', struct('comparisons', {{'halves'}}),                 'PC sweep: within-area halves',              6.2};
    case 'test'
        t = fullfile(tempdir, 'neuron_decoding_rev_test');
        jobs = {
            'run_cue_transfer_incremental', struct('populations', {{'FR','V1'}}, 'N_list', 48, 'n_perm', 2, ...
                'out_root', t, 'overwrite', true), 'test: cue transfer, FR V1 at full size', 0.02
            'run_halves_incremental', struct('populations', {{'FR','V1'}}, 'h_list', 24, 'n_draws', 1, 'n_perm', 2, ...
                'out_root', t, 'overwrite', true), 'test: halves, FR V1', 0.01
            'run_across_area', struct('monkeys', {{'FR'}}, 'n_draws', 1, 'n_perm', 2, ...
                'out_root', t, 'overwrite', true), 'test: across areas, FR', 0.03
            'run_pc_sweep', struct('populations', {{'FR','V1'}}, 'monkeys', {{'FR'}}, 'k_base', 5, ...
                'n_draws_halves', 1, 'n_perm', 2, 'out_root', t, 'overwrite', true), 'test: PC sweep, FR', 0.08};
    otherwise
        error('run_all: unknown set ''%s'' (use all, desktop, laptop or test).', run_set);
end
n_jobs = size(jobs, 1);

%% Log file
if strcmpi(run_set, 'test')
    log_dir = fullfile(tempdir, 'neuron_decoding_rev_test');
else
    log_dir = fullfile(code_dir, '..', '..', 'results', 'logs');
end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end
log_file = fullfile(log_dir, sprintf('run_all_%s_%s.log', lower(run_set), char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
diary(log_file);
cleanup = onCleanup(@() diary('off'));

fprintf('run_all(''%s'') started %s on %s\nLog: %s\n\n', run_set, char(datetime('now')), getenv('COMPUTERNAME'), log_file);
for j = 1:n_jobs
    fprintf('  job %d/%d  %-30s %-45s ~%.1f h\n', j, n_jobs, jobs{j, 1}, jobs{j, 3}, jobs{j, 4});
end
fprintf('  total ~%.1f h at the 2026-10-11 speed (8 workers)\n', sum([jobs{:, 4}]));

if isempty(gcp('nocreate')), parpool('Processes'); end
p = gcp;
fprintf('Parallel pool: %d workers\n', p.NumWorkers);

%% Run the jobs one after another
status = repmat({'not run'}, n_jobs, 1);
t_all = tic;
for j = 1:n_jobs
    fprintf('\n################ job %d/%d: %s (%s) | started %s ################\n', ...
        j, n_jobs, jobs{j, 3}, jobs{j, 1}, char(datetime('now')));
    t0 = tic;
    try
        run_driver(jobs{j, 1}, jobs{j, 2});
        status{j} = sprintf('done in %.2f h', toc(t0) / 3600);
    catch ME
        status{j} = sprintf('FAILED after %.2f h: %s', toc(t0) / 3600, ME.message);
        fprintf(2, '\nJob %d failed: %s\n%s\n', j, ME.message, getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

fprintf('\n################ run_all(''%s'') finished %s, %.2f h ################\n', run_set, char(datetime('now')), toc(t_all) / 3600);
for j = 1:n_jobs
    fprintf('  job %d  %-45s %s\n', j, jobs{j, 3}, status{j});
end
fprintf('Log: %s\n', log_file);
end
