function cfg = b1_p1_config()
% P1 of the low-q review (2026-10-09, scope approved by the user after the
% GPT-Pro review): background constrained by a Stage A prefit (ZLP + aux
% peaks with their detailed-balance gain side, on the gain side, the core and
% the loss side up to 150 meV); the 150-250 meV band held out of every
% data-dependent term (prefit, joint data term, overshoot, noise estimate)
% and used only to validate the extrapolated background; loss-peak bounds
% separated from the data windows; narrow aux peaks; instrument kernel fixed
% to the Stage A ZLP. Everything else as the previous baseline (stage2_2d_h5:
% 2D kinematic factor, h = 0.0005, prefactor held below 250 meV).
cfg.version = 'P1-2026-10-09';
cfg.fit_window = [-Inf 1800];
cfg.signal_window = [250 1800];          % B1 statistics (data term: everything outside the holdout)
cfg.holdout_windows = [150 250];         % validation band, no data-dependent term
cfg.peak_model = 'lorentz'; cfg.n_peaks = 2;
cfg.peak_E0_bounds = [150 2500];         % parameter bounds, independent of the data windows
cfg.peak_width_max = 3000;
cfg.aux_windows = [30 60; 80 150];       % ~43 meV and the 90-140 meV features
cfg.aux_width_max = [40; 60];
cfg.aux_model = 'lorentz';
cfg.n_zlp = 2; cfg.asym_sigma = 0.15; cfg.gain_weight = 3; cfg.core_weight = 0.05;
cfg.overshoot = true; cfg.overshoot_sigma = 3;
cfg.prefit_aux = true; cfg.prefit_loss_max = 150;
cfg.kernel_mode = 'prefit'; cfg.kernel_halfwidth = 300;
cfg.kinematic_form = '2d'; cfg.h_perp = 0.0005; cfg.beam_kV = 30;
cfg.prefactor_floor = 250;               % unchanged in P1 (its removal belongs to P2)
cfg.prefactor_on_aux = false;
E0 = [400 650; 600 1000; 800 1600];      % deterministic starts (meV): E0 pairs x width pairs
W = [150 250; 400 650; 900 1500];
[i, j] = ndgrid(1:size(E0, 1), 1:size(W, 1));
cfg.start_E0 = E0(i(:), :); cfg.start_width = W(j(:), :);
cfg.add_baseline_start = true;           % the baseline E0/width as one more start
cfg.n_random_starts = 8;
cfg.seed = 20261009;
cfg.n_boot = 20; cfg.boot_random_starts = 2;
cfg.score_windows = struct('name', {'gain_-180_-30', 'loss_30_150', 'holdout_150_250', 'b1_250_300', 'b1_250_1800', 'b1_300_1800'}, ...
    'range', {[-180 -30], [30 150], [150 250], [250 300], [250 1800], [300 1800]});
cfg.baseline_stage = 'stage2_2d_h5';
cfg.retardation_flag_ratio = 0.5;        % flag q_gamma(E)/|q| above this as "to check"
end
