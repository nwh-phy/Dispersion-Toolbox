function cfg = b1_low_branch_matlab_config(data_root, ids)
% MATLAB-only B1 analysis, using the existing toolbox rather than a port.
% ids: datasets of bisb_datasets (default the three 20260120 sessions). The
% session fields keep the order of the stored prepared.mat signatures.
arguments
 data_root (1,:) char {mustBeFolder}
 ids (1,:) cell = {'590','n0','20w'}
end
cfg.version='B1-low-MATLAB-20261024';cfg.seed=20261024;
cfg.data_root=data_root;
cfg.session=struct('id',{},'folder',{},'file',{},'dq',{},'zero',{},'h',{},'sha256',{});
for k=1:numel(ids)
 d=bisb_datasets(ids{k});
 assert(isfinite(d.h),'b1_low_branch_matlab_config:NoSlot', ...
  'Dataset %s has no slot acceptance h; enter it in bisb_datasets first',d.id);
 cfg.session(k)=struct('id',d.short,'folder',d.folder,'file',d.file,'dq',d.dq,'zero',d.zero, ...
  'h',d.h,'sha256',d.sha256);
end
cfg.preparation=struct('native_max',.006,'q_max',.020,'q_min',.0005, ...
 'high_channels',[3 5],'reference_q',.001,'alignment_window',[-100 100], ...
 'raw_zlp_search',[-350 100],'energy_range',[-300 1800],'do_q_drift',true, ...
 'q_centroid_range',.002,'q_centroid_energy',16,'q_median_frames',11,'train_frames',250);
base=b1_fullq_config();
cfg.fit=struct('fit_window',[-300 1800],'signal_window',base.window, ...
 'peak_model','lorentz','n_zlp',base.n_zlp,'aux_windows',base.aux{1,2}, ...
 'n_starts',4,'prefactor_floor',base.prefactor_floor,'beam_kV',30, ...
 'peak_E0_bounds',[300 1800],'peak_width_max',3000,'max_iterations',250, ...
 'n_sub',53,'n_probe',25,'seed',cfg.seed);
cfg.bootstrap=struct('targets',[.001 .005],'repeats',32,'block_frames',10,'n_starts',1);
cfg.review=struct('weak_area_fraction',.05,'merge_gap_meV',100);
end
