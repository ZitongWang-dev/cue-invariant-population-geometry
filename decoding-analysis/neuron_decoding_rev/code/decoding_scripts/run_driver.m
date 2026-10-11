function run_driver(driver, override)
%{
Filename: run_driver.m
Author:   Zitong Wang
Date:     2026-10-11

Description:
    Runs one driver script (e.g. 'run_pc_sweep') inside this function's own
    workspace, with an optional cfg_override. Used by run_all.m: the drivers
    start with clearvars, which here clears only this function's workspace, so
    the master script's job list survives, and a cfg_override left in the base
    workspace is never picked up by accident.

Inputs:
    driver   : name of a driver script in this folder, without '.m'
    override : struct of cfg fields to override, or [] for the defaults
%}

if nargin > 1 && ~isempty(override)
    cfg_override = override; %#ok<NASGU> read by the driver script
end
run(fullfile(fileparts(mfilename('fullpath')), [driver '.m']));
end
