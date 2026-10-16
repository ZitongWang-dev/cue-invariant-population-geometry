function run_all_gabor(run_set)
%{
Filename: run_all_gabor.m
Author:   Zitong Wang
Date:     2026-10-15

Description:
    Master script of Gabor_decoding_rev: runs the Gabor-model decoding drivers
    one after another, unattended, and writes everything printed to a log file
    in results/logs/. Same behaviour as neuron_decoding_rev/run_all.m and
    ANN_decoding_rev/run_all_ann.m.

        run_all_gabor('all')    both drivers (about 5 h on one 8-worker machine
                                at the CNN speed measured on 2026-10-12)
        run_all_gabor('test')   tiny configurations (about 5 min) written to a
                                temporary folder; run it once on the machine first

    run_all_gabor with no argument uses DEFAULT_SET below, so pressing Run in
    the editor also works.

    Each job runs in its own workspace (run_gabor_driver.m). A job that stops
    with an error is logged and skipped, and the next job starts. Every driver
    resumes where it stopped, so after fixing a problem just call run_all_gabor
    again with the same set.
%}

DEFAULT_SET = 'all';
if nargin < 1 || isempty(run_set), run_set = DEFAULT_SET; end

code_dir = fileparts(mfilename('fullpath'));
addpath(code_dir);

%% Jobs: {driver, cfg_override ([] = defaults), description, hours at the measured speed}
switch lower(run_set)
    case 'all'
        jobs = {
            'run_gabor_cue_transfer_incremental', [],   'cue transfer vs cells, simple and complex (Fig 15, point 9)', 4.5
            'run_gabor_pc_sweep',                 [],   'PC sweep: simple and complex',                                0.3};
    case 'test'
        t = fullfile(tempdir, 'Gabor_decoding_rev_test');
        jobs = {
            'run_gabor_cue_transfer_incremental', struct('N_list', [48 192], 'n_draws', 1, 'n_perm', 2, ...
                'out_root', t, 'overwrite', true), 'test: cue transfer at 48 (Fig 15a) and whole models (Fig 15b)', 0.03
            'run_gabor_pc_sweep', struct('k_base', 5, 'n_perm', 2, 'out_root', t, 'overwrite', true), ...
                'test: PC sweep at k = 5 and the maximum k', 0.05};
    otherwise
        error('run_all_gabor: unknown set ''%s'' (use all or test).', run_set);
end
n_jobs = size(jobs, 1);

%% Log file
if strcmpi(run_set, 'test')
    log_dir = fullfile(tempdir, 'Gabor_decoding_rev_test');
else
    log_dir = fullfile(code_dir, '..', '..', 'results', 'logs');
end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end
log_file = fullfile(log_dir, sprintf('run_all_gabor_%s_%s.log', lower(run_set), char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
diary(log_file);
cleanup = onCleanup(@() diary('off'));

fprintf('run_all_gabor(''%s'') started %s on %s\nLog: %s\n\n', run_set, char(datetime('now')), getenv('COMPUTERNAME'), log_file);
for j = 1:n_jobs
    fprintf('  job %d/%d  %-36s %-62s ~%.1f h\n', j, n_jobs, jobs{j, 1}, jobs{j, 3}, jobs{j, 4});
end
fprintf('  total ~%.1f h (8 workers)\n', sum([jobs{:, 4}]));

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
        run_gabor_driver(jobs{j, 1}, jobs{j, 2});
        status{j} = sprintf('done in %.2f h', toc(t0) / 3600);
    catch ME
        status{j} = sprintf('FAILED after %.2f h: %s', toc(t0) / 3600, ME.message);
        fprintf(2, '\nJob %d failed: %s\n%s\n', j, ME.message, getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

fprintf('\n################ run_all_gabor(''%s'') finished %s, %.2f h ################\n', run_set, char(datetime('now')), toc(t_all) / 3600);
for j = 1:n_jobs
    fprintf('  job %d  %-62s %s\n', j, jobs{j, 3}, status{j});
end
fprintf('Log: %s\n', log_file);
end
