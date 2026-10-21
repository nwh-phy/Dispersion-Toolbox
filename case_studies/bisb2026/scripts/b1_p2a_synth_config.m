function cfg = b1_p2a_synth_config()
% Configuration of P2a step 2.4 (injection and recovery with the 590 block operators).
cfg.version = 'P2a-step2.4d-2026-10-15';   % 20 noise seeds   % + weighting test of the noisy fits
cfg.fit_window = [250 1800];                 % B1 data term
cfg.holdout = [150 250];                     % extrapolation check, not fitted
cfg.gen = struct('n_a', 8, 'n_s', 17, 'n_gh', 61, 'Qf', 2);   % finer generator discretisation
cfg.variants = {'slice', 'integrated', 'exit'};
cfg.lambda_curv = 10;                        % weak curvature penalty on log E0^2 (second differences over the nodes)
cfg.bounds_E0 = [300 1300; 900 2500];        % meV, low / high branch
cfg.bounds_G = [50 3000];                    % meV
cfg.start_perturb = struct('E0', 1.05, 'G', 1.2);
cfg.noise_seeds = 1:20;
cfg.noise_class = 3;
cfg.noise_variants = {'slice', 'exit'};
cfg.noise_weights = {'data', 'smoothed', 'oracle'};   % sigma^2 = counts / 21-point moving mean of the counts / expected counts
cfg.smooth_points = 21;
cfg.check = struct('start_perturb', [1.0 1.0], 'FunctionTolerance', 1e-15, 'StepTolerance', 1e-12, 'MaxIterations', 500, ...
    'weights', 'oracle');            % convergence check of the noisy fits: start at the truth, tight tolerances
cfg.classes = struct( ...
    'name', {'fixed E0, widths and weights vary with Q', 'known dispersion', 'reference L (P1b-like)'}, ...
    'E0', {[760 760 760 760 760; 1620 1620 1620 1620 1620], [500 650 720 760 780; 1300 1450 1550 1620 1650], []}, ...
    'G', {[300 450 600; 600 1200 1800], [400 500 550; 900 1300 1600], []}, ...
    'a', {[0.004 0.0015 0.0006 0.08 0.04 0.025], [], []});
cfg.report_q = [0.0005 0.001 0.002 0.003 0.005];
end
