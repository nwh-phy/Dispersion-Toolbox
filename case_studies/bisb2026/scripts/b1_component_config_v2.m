function cfg = b1_component_config_v2()
% Project-only policy; registry remains the source of session dq and paths.
cfg=struct(); cfg.version='2.0'; cfg.sessions=thesis_sessions();
cfg.q_range_Ainv=[-0.015 0.015]; cfg.q_skip_Ainv=0.0005;
cfg.native_N=[1 3 5]; cfg.matched_20w_N=[6 10];
cfg.energy_windows_meV=[300 1800;300 2000;300 2100];
cfg.models={'lorentz_symmetric','lorentz'}; cfg.n_starts=24; cfg.seed=20260912;
cfg.baseline_mode='power_law'; cfg.noise_model='unweighted_LS_unknown_correlated_detector_noise';
cfg.align_zlp=false; cfg.runFits=false; cfg.jump_repair=false;
cfg.invalid_count_threshold=double(intmax('uint32'));
cfg.detector_mask_policy='any frame uint32 saturation sentinel or nonfinite excludes entire native q channel';
cfg.weak_peak_deletion_fraction=0; cfg.trend_constraints=false;
cfg.processing_level='L1'; cfg.response='unknown_no_deconvolution';
cfg.representative_q_Ainv=[-0.0025 0.0075 0.0125];
cfg.representative_policy='fixed low/mid/high signed-q targets before fitting; nearest valid full bin';
end
