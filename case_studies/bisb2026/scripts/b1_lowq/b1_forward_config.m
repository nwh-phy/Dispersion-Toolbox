function paths = b1_forward_config(session, run_root)
%B1_FORWARD_CONFIG  Dataset and file locations shared by the P1b/P2a/P2b stages.
% Attach this struct as cfg.paths. For a new dataset, set input_dir, response_dir,
% operator_dir and fit_dir from the outputs of preceding stages as they become available.
if nargin < 1 || isempty(session), session = '590_PL2_10w'; end
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
d = bisb_datasets(session);
if nargin < 2 || isempty(run_root)
    run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z');
end
paths = struct('session', d.id, 'run_root', char(run_root));
paths.stage1_file = fullfile(run_root, 'stage1', sprintf('channels_%s.mat', d.id));
paths.diag_dir = fullfile(run_root, 'p1b', ['diag_' d.id]);
paths.p1b_dir = fullfile(run_root, 'p1b');
paths.p2a_dir = fullfile(run_root, 'p2a');
paths.p2b_dir = fullfile(run_root, 'p2b');
paths.l1_file = fullfile(root, 'paper_results', 'b1_components_v2', ...
    '20260912T021813081Z_31fcc543_18af79b6', d.id, 'L1_minimal.mat');
paths.highq_file = fullfile(run_root, 'stage5', sprintf('e0_stage2_2d_h5_%s.csv', d.id));
paths.fit_prefix = ['fit_' d.id '_'];
paths.legacy_fit_prefix = '';
paths.input_dir = '';
paths.response_dir = '';
paths.operator_dir = '';
paths.fit_dir = '';
if strcmp(d.id, '590_PL2_10w')
    paths.input_dir = fullfile(paths.p2a_dir, 'inputs_590_PL2_10w_66385d96');
    paths.response_dir = fullfile(paths.input_dir, 'response_42ffc295');
    paths.operator_dir = fullfile(paths.response_dir, 'operator_e9052bc0');
    paths.fit_dir = fullfile(paths.p2b_dir, 'fit_590_691c36ad');
    paths.legacy_fit_prefix = 'fit_590_';
end
end
