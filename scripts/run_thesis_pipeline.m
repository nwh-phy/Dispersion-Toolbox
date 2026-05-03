%% run_thesis_pipeline.m — baseline thesis q-EELS analysis wrapper
% Run from the repository root with:
%   run('scripts/run_thesis_pipeline.m')
%
% The implementation lives in src/thesis/run_thesis_pipeline.m so it can be
% called from tests and from other scripts.

clearvars; clc;
project_dir = fileparts(fileparts(mfilename('fullpath')));
cd(project_dir);
run(fullfile(project_dir, 'startup.m'));

cfg = thesis_config();
run_thesis_pipeline(cfg);
fprintf('\nThesis pipeline complete. Output: %s\n', cfg.output.baseline_dir);
