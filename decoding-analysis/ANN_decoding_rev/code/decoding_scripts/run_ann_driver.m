function run_ann_driver(driver, override)
%{
Filename: run_ann_driver.m
Author:   Zitong Wang
Date:     2026-10-12

Description:
    Runs one driver script of this folder (e.g. 'run_ann_pc_sweep') inside this
    function's own workspace, with an optional cfg_override. Used by
    run_all_ann.m: the drivers start with clearvars, which here clears only this
    function's workspace, so the master script's job list survives, and a
    cfg_override left in the base workspace is never picked up by accident.
    (Same as neuron_decoding_rev/run_driver.m, which looks for drivers in its
    own folder.)

Inputs:
    driver   : name of a driver script in this folder, without '.m'
    override : struct of cfg fields to override, or [] for the defaults
%}

if nargin > 1 && ~isempty(override)
    cfg_override = override; %#ok<NASGU> read by the driver script
end
run(fullfile(fileparts(mfilename('fullpath')), [driver '.m']));
end
