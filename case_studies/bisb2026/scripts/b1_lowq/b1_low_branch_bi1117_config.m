function cfg = b1_low_branch_bi1117_config(data_root, h)
%B1_LOW_BRANCH_BI1117_CONFIG  b1_low_branch_matlab on the 2026-11-17 Bi 20w sequence.
%   cfg = b1_low_branch_bi1117_config(data_root) is a modified copy of
%   b1_low_branch_matlab_config(data_root, {'Bi1117_20w'}): the first 300 of the
%   360 frames train the one- and two-B1 fits and the last 60 are held out (the
%   20260120 runs used 250 of 300 frames). Everything else (binning, joint ZLP
%   model, kinematic factor, bootstrap) is unchanged.
%   cfg = b1_low_branch_bi1117_config(data_root, h) replaces the slot acceptance
%   h (1/A) of the registry, e.g. h = 0 for the group convention; use a separate
%   output folder for every h.
%
%   See also: b1_low_branch_matlab, b1_low_branch_matlab_config, bisb_datasets
arguments
    data_root (1,:) char {mustBeFolder}
    h (1,1) double = NaN
end
cfg = b1_low_branch_matlab_config(data_root, {'Bi1117_20w'});
cfg.version = 'B1-low-MATLAB-Bi1117-20261117';
cfg.preparation.train_frames = 300;
if ~isnan(h)
    assert(h >= 0, 'b1_low_branch_bi1117_config:H', 'h must be nonnegative');
    cfg.session.h = h;
end
end
