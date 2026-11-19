function cfg = b1_forward_prep_config(session, run_root)
%B1_FORWARD_PREP_CONFIG  Stage1/P1b acquisition settings for the four Bi runs.
if nargin<2 || isempty(run_root)
    root=bisb_find_project_root(fileparts(mfilename('fullpath')));
    run_root=fullfile(root,'paper_results','b1_forward_bi_20261119');
end
cfg.version='P0-2026-11-19';
cfg.paths=b1_forward_config(session,run_root);
d=bisb_datasets(cfg.paths.session);
switch d.id
    case 'Bi1107_10w',  cfg.frames_per_block=60; cfg.qmax=0.035;
    case 'Bi1107_20w',  cfg.frames_per_block=30; cfg.qmax=0.0125;
    case 'Bi1117_10w',  cfg.frames_per_block=30; cfg.qmax=0.035;
    case 'Bi1117_20w',  cfg.frames_per_block=60; cfg.qmax=0.0125;
    otherwise, error('b1_forward_prep_config:Session','Unsupported P0 acquisition %s.',d.id);
end
cfg.n_blocks=6;
p=b1_fullq_config();
cfg.reference_abs_q=p.reference_abs_q;
cfg.zlp_window=p.zlp_window;
cfg.raw_profile_window=12;
cfg.holdout=[150 250];
cfg.noise_q=[0.0005 0.001 0.002 0.005];
assert(d.n_frames==cfg.n_blocks*cfg.frames_per_block,'Unexpected registry frame count.');
end
