function run_all_ann(run_set)
%{
Filename: run_all_ann.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Master script of ANN_decoding_rev: runs the CNN decoding drivers one after
    another, unattended, and writes everything printed to a log file in
    results/logs/. Same behaviour as neuron_decoding_rev/run_all.m.

        run_all_ann('all')       every driver (about 11 h on one 8-worker machine
                                 at the speed measured on 2026-10-12)
        run_all_ann('desktop')   one half of a two-machine split: the four large
                                 layers vs units, the hierarchy and its PC sweep
                                 (about 6.5 h at that speed)
        run_all_ann('laptop')    the other half: the four small layers vs units and
                                 the cue PC sweep of all layers (about 5 h)
        run_all_ann('test')      tiny configurations (about 6 min) written to a
                                 temporary folder; run it once on each machine first

    run_all_ann with no argument uses DEFAULT_SET below, so pressing Run in the
    editor also works.

    Each job runs in its own workspace (run_ann_driver.m). A job that stops with
    an error is logged and skipped, and the next job starts. Every driver
    resumes where it stopped, so after fixing a problem just call run_all_ann
    again with the same set.
%}

DEFAULT_SET = 'all';
if nargin < 1 || isempty(run_set), run_set = DEFAULT_SET; end

code_dir = fileparts(mfilename('fullpath'));
addpath(code_dir);

big   = {'alexnet',5; 'vgg',3; 'vgg',4; 'vgg',5};       % 249-511 units: grid up to 160 / 320
small = {'alexnet',1; 'alexnet',2; 'vgg',1; 'vgg',2};    % 64-250 units

%% Jobs: {driver, cfg_override ([] = defaults), description, hours at the measured speed}
switch lower(run_set)
    case 'all'
        jobs = {
            'run_ann_cue_transfer_incremental', [],                                 'cue transfer vs units, 8 layers (Fig 3d, 4d)', 9.5
            'run_ann_hierarchy',                [],                                 'AlexNet pool2 <-> pool5',                      0.2
            'run_ann_pc_sweep',                 [],                                 'PC sweep: cue transfer and hierarchy',         1.7};
    case 'desktop'
        jobs = {
            'run_ann_cue_transfer_incremental', struct('layers', {big}),            'cue transfer vs units: AlexNet pool5, VGG pool3-5', 5.8
            'run_ann_hierarchy',                [],                                 'AlexNet pool2 <-> pool5',                      0.2
            'run_ann_pc_sweep', struct('comparisons', {{'hierarchy'}}),             'PC sweep: hierarchy',                          0.5};
    case 'laptop'
        jobs = {
            'run_ann_cue_transfer_incremental', struct('layers', {small}),          'cue transfer vs units: AlexNet pool1-2, VGG pool1-2', 3.8
            'run_ann_pc_sweep', struct('comparisons', {{'cue'}}),                   'PC sweep: cue transfer, 8 layers',             1.2};
    case 'test'
        t = fullfile(tempdir, 'ANN_decoding_rev_test');
        jobs = {
            'run_ann_cue_transfer_incremental', struct('layers', {{'vgg',1}}, 'N_list', 64, 'n_perm', 2, ...
                'out_root', t, 'overwrite', true), 'test: cue transfer, VGG pool1 whole layer (Fig 3d point)', 0.01
            'run_ann_hierarchy', struct('n_draws', 1, 'n_perm', 2, 'out_root', t, 'overwrite', true), ...
                'test: AlexNet pool2 <-> pool5, 1 draw', 0.03
            'run_ann_pc_sweep', struct('layers', {{'vgg',1}}, 'k_base', 5, 'n_perm', 2, ...
                'out_root', t, 'overwrite', true), 'test: PC sweep, VGG pool1 and hierarchy', 0.06};
    otherwise
        error('run_all_ann: unknown set ''%s'' (use all, desktop, laptop or test).', run_set);
end
n_jobs = size(jobs, 1);

%% Log file
if strcmpi(run_set, 'test')
    log_dir = fullfile(tempdir, 'ANN_decoding_rev_test');
else
    log_dir = fullfile(code_dir, '..', '..', 'results', 'logs');
end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end
log_file = fullfile(log_dir, sprintf('run_all_ann_%s_%s.log', lower(run_set), char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
diary(log_file);
cleanup = onCleanup(@() diary('off'));

fprintf('run_all_ann(''%s'') started %s on %s\nLog: %s\n\n', run_set, char(datetime('now')), getenv('COMPUTERNAME'), log_file);
for j = 1:n_jobs
    fprintf('  job %d/%d  %-34s %-52s ~%.1f h\n', j, n_jobs, jobs{j, 1}, jobs{j, 3}, jobs{j, 4});
end
fprintf('  total ~%.1f h at the 2026-10-12 speed (8 workers)\n', sum([jobs{:, 4}]));

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
        run_ann_driver(jobs{j, 1}, jobs{j, 2});
        status{j} = sprintf('done in %.2f h', toc(t0) / 3600);
    catch ME
        status{j} = sprintf('FAILED after %.2f h: %s', toc(t0) / 3600, ME.message);
        fprintf(2, '\nJob %d failed: %s\n%s\n', j, ME.message, getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end

fprintf('\n################ run_all_ann(''%s'') finished %s, %.2f h ################\n', run_set, char(datetime('now')), toc(t_all) / 3600);
for j = 1:n_jobs
    fprintf('  job %d  %-52s %s\n', j, jobs{j, 3}, status{j});
end
fprintf('Log: %s\n', log_file);
end
