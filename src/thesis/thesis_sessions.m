function sessions = thesis_sessions(cfg)
%THESIS_SESSIONS  Dataset registry for the thesis q-EELS pipeline.
%   The registry keeps sample labels, folder paths, and dq calibration in one
%   place so scripts do not hard-code local dataset names in multiple files.

if nargin < 1 || isempty(cfg)
    cfg = thesis_config();
end
root_dir = cfg.project_root;
base_dir = fullfile(root_dir, '20260120 BiSb');

sessions = repmat(struct( ...
    'name', '', ...
    'role', '', ...
    'path', '', ...
    'dq_Ainv', NaN, ...
    'enabled', true, ...
    'sample_label', cfg.notes.sample_label, ...
    'naming_caution', cfg.notes.naming_caution), 1, 3);

sessions(1).name = '590_PL2_10w';
sessions(1).role = 'primary_10w_defocus';
sessions(1).path = fullfile(base_dir, '590 PL2 10w 0.004 10sx300');
sessions(1).dq_Ainv = 0.005;

sessions(2).name = 'no_PL2_20w_2film';
sessions(2).role = 'primary_20w_defocus';
sessions(2).path = fullfile(base_dir, 'no pl2 20w 0.004 10sx300 2film');
sessions(2).dq_Ainv = 0.0025;

sessions(3).name = 'n0_PL2_10w_repeat';
sessions(3).role = 'repeat_10w_defocus_sanity_check';
sessions(3).path = fullfile(base_dir, 'n0 pl2 10w 0.004 10s x300');
sessions(3).dq_Ainv = 0.005;
end
