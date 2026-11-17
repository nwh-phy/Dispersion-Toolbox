function cfg = b1_low_branch_bi_config(data_root, id)
%B1_LOW_BRANCH_BI_CONFIG  b1_low_branch_matlab on one 2026-11 Bi sequence (layer check).
%   cfg = b1_low_branch_bi_config(data_root, id) is a modified copy of
%   b1_low_branch_matlab_config(data_root, {id}): the last sixth of the frames is
%   held out (300 + 60 of 360, 150 + 30 of 180; the 20260120 runs held out 50 of
%   300, b1_low_branch_bi1117_config 60 of 360). Everything else (binning, joint
%   ZLP model, kinematic factor, bootstrap) is unchanged; h comes from bisb_datasets.
%
%   See also: b1_low_branch_matlab, b1_low_branch_matlab_config, b1_low_branch_bi1117_config
arguments
    data_root (1,:) char {mustBeFolder}
    id (1,:) char
end
d = bisb_datasets(id);
cfg = b1_low_branch_matlab_config(data_root, {id});
cfg.version = 'B1-low-MATLAB-Bi-20261117';
cfg.preparation.train_frames = d.n_frames - round(d.n_frames / 6);
% single-channel limit on a whole channel (0.006 = 47.0 channels at dq 0.000128 would cut the
% first merged bin; unchanged at dq 0.0005 and 0.00025)
cfg.preparation.native_max = round(cfg.preparation.native_max / d.dq) * d.dq;
end
