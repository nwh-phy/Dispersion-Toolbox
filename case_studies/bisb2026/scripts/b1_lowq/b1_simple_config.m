function cfg = b1_simple_config()
% Configuration of the simplified pipeline check (docs/plans/simple_pipeline_20261017/PLAN.md):
% SI-style processing of He et al. (frame alignment, sum, per-q single-spectrum fit) plus the
% single-channel kinematic factor, compared with P2b. Settings shared with P2b are read from
% b1_p2b_config so that the comparison only tests the simplifications.
p2b = b1_p2b_config();
cfg.version = 'simple-S0-2026-10-17';
cfg.p2b_dir = fullfile('p2b', 'fit_590_691c36ad');
cfg.window = p2b.variants(1).window;        % 250-1800 meV, as P2b main
cfg.tail = p2b.tail;                        % (E/250)^-4 per spectrum, projected out
cfg.noise_smooth = p2b.noise_smooth;        % sigma^2 = 21-point moving mean of the counts
cfg.bounds_E0 = p2b.bounds_E0; cfg.bounds_G = p2b.bounds_G;
cfg.starts = p2b.starts; cfg.nodes = p2b.nodes;   % starts: P2b node values, interpolated to each channel
cfg.kT = 25.852;                            % meV, thermal factor before the convolution (as P2b)
cfg.res_fwhm = 23.2;                        % meV, angle-integrated zero-loss line (P2a step 1), Gaussian here
cfg.res_half = 72;                          % meV, kernel support
cfg.dq = 0.0005;
% S0: channels |q| >= 0.0015 of the P2b truth; quantities judged against the P2b statistical error
cfg.s0_qmin_test = 0.0015;
cfg.judged = {'E0_b1', 'apex_b1', 'apex_b2'};      % E0 of the high branch is reported, not judged (PLAN S3)
% simplified angular response, first row = the plan: single channel, Gaussian probe = narrow
% component of the two-component profile (P2a step 1: FWHM 8.0e-4 -> sigma 3.40e-4), no broad
% component. Other rows change one thing each, to attribute a bias if one appears.
cfg.s0_variants = struct( ...
    'name', {'plan', 'y as P2a (h 0.0005, sigma_y 4.5e-4)', 'residue prop. to Q inside the bin', 'no probe spread'}, ...
    'sigma', {3.40e-4, 3.40e-4, 3.40e-4, 0}, ...
    'h_perp', {0, 0.0005, 0, 0}, ...
    'sigma_y', {NaN, 4.5e-4, NaN, NaN}, ...
    'qweight', {false, false, true, false});
cfg.fig_q = [0.0015 0.0025 0.0035 0.005];   % channels shown in the spectra figure (both signs)
end
