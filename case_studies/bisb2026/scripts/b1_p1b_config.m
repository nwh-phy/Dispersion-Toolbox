function cfg = b1_p1b_config()
% P1b of the low-q review (2026-10-10; plan: docs/plans/lowq_background_20261010/
% PLAN.md, from GPT Pro round 2). Same data as P1 (590, |q| = 0.0005, 0.001,
% 0.002, 0.005, both signs). Changes against P1:
%  - one response for elastic and inelastic signal (zlp_mode 'prefit_shape'):
%    Stage B ZLP = a*R(E - c), inelastic terms convolved with the same R at the
%    same c; R fixed from Stage A;
%  - Stage A back to the gain side + elastic core, aux peaks with their full
%    +/- model, overshoot restricted to the Stage A data; two kernel versions:
%    KA = core-constrained (gain side >= -60 meV, loss side <= 20 meV),
%    KB = gain-side-constrained full kernel (whole gain side, loss side up to
%    one FWHM, as the original Stage A);
%  - noise: common-beam reference (whole frame, holdout excluded) and the
%    first-difference variance coefficient (slow drifts removed), holdout
%    masked in the running median;
%  - 250 meV prefactor floor kept for the main fits (f250); floor 0 only as a
%    fixed-parameter forward comparison, a quadrature check and short fits (f0).
% Unchanged: 2D kinematic factor, h = 0.0005, holdout 150-250 meV, B1 bounds.
cfg.version = 'P1b-2026-10-10';
cfg.fit_window = [-1e9 1800];            % finite lower bound (JSON-safe), = everything below 1800 meV
cfg.signal_window = [250 1800];
cfg.holdout_windows = [150 250];
cfg.peak_model = 'lorentz'; cfg.n_peaks = 2;
cfg.peak_E0_bounds = [150 2500]; cfg.peak_width_max = 3000;
cfg.aux_windows = [35 60; 70 150];       % ~45 meV feature; 80-140 meV features (weakly constrained)
cfg.aux_width_max = [30; 60];
cfg.aux_model = 'lorentz';
cfg.n_zlp = 2; cfg.asym_sigma = 0.15; cfg.gain_weight = 3; cfg.core_weight = 0.05;
cfg.overshoot = true; cfg.overshoot_sigma = 3;
cfg.zlp_mode = 'prefit_shape'; cfg.kernel_iterations = 5; cfg.kernel_halfwidth = 300;
cfg.prefit_aux = true; cfg.prefit_overshoot_local = true;
cfg.kernels = struct('name', {'KA', 'KB'}, 'prefit_gain_min', {-60, -Inf}, 'prefit_loss_max', {20, NaN}, ...
    'note', {'core-constrained: gain side >= -60 meV, loss side <= 20 meV', ...
             'gain-side-constrained full kernel: whole gain side, loss side <= one FWHM'});
cfg.kinematic_form = '2d'; cfg.h_perp = 0.0005; cfg.beam_kV = 30;
cfg.floors = [250 0];                    % main fits at 250 (as P1); 0 = short fits after the checks
cfg.prefactor_on_aux = false;
cfg.quadrature_lo = struct('n_sub', 9, 'n_probe', 13, 'probe_extent', 3, 'q_round', 1e-6);
cfg.quadrature_hi = struct('n_sub', 27, 'n_probe', 41, 'probe_extent', 4.5, 'q_round', 1e-7);
cfg.quadrature_tol_sigma = 0.1;          % acceptance: |change| < 0.1 sigma in fit and holdout regions
cfg.noise = struct('reference', 'whole frame, valid channels, all energies except the holdout', ...
    'estimator', 'first differences of common-beam corrected frames', 'median_points', 21);
E0 = [400 650; 600 1000; 800 1600]; W = [150 250; 400 650; 900 1500];
[i, j] = ndgrid(1:size(E0, 1), 1:size(W, 1));
cfg.start_E0 = E0(i(:), :); cfg.start_width = W(j(:), :);
cfg.add_p1_start = true;
cfg.n_random_starts = 4; cfg.n_random_starts_short = 2;
cfg.seed = 20261010;
cfg.score_windows = struct('name', {'gain_-180_-30', 'loss_30_150', 'holdout_150_250', 'b1_250_300', 'b1_250_1800', 'b1_300_1800'}, ...
    'range', {[-180 -30], [30 150], [150 250], [250 300], [250 1800], [300 1800]});
cfg.share_energies = [198 300 392];
cfg.p1_dir = 'p1/590_PL2_10w_12f9932f';
cfg.baseline_stage = 'stage2_2d_h5';
cfg.perturbation = struct('q', {0.0005, -0.002}, 'scale', 1.3, 'offset_counts', 50);
end
