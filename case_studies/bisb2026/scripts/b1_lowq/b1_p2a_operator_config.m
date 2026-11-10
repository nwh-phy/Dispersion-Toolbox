function cfg = b1_p2a_operator_config()
% Configuration of P2a step 2.2-2.5 (operator, checks, fixed-L forward diagnostics).
cfg.version = 'P2a-step2.2c-2026-10-15';   % y nodes 31 (main), delta-chi2 before/after residue refit
cfg.exit_q = [-0.005 -0.003 -0.0025 -0.002 -0.0015 -0.001 -0.0005 0.0005 0.001 0.0015 0.002 0.0025 0.003 0.005];
cfg.h = 0.0005;                    % half acceptance across the dispersion direction (A^-1), assumption
cfg.sigma_y = 4.5e-4;              % Gaussian y spread of the incident beam (A^-1), main assumption (review round 4)
cfg.n_a = 4;                       % Gauss-Legendre nodes per incident slice
cfg.n_sx = 9; cfg.n_sy = 9;        % exit-box samples
cfg.n_gh = 31;                     % Gauss-Hermite nodes in y (Gaussian beam)
cfg.Qm_spec = [0 1e-5 2e-4; 2e-4 2.5e-5 2e-3; 2e-3 1e-4 1e-2; 1e-2 5e-4 6e-2];   % radial grid segments [from step to]
cfg.tail_edges = [0.00675 0.00875 0.01075 0.01475 0.02 0.03 0.04];               % incident tail slices beyond the measured range
cfg.omega = [-400 4 2400];         % loss grid (meV)
cfg.eps = [-72 4 72];              % kernel offset grid (meV, absolute energy of the slice lines)
cfg.support = 60;                  % kernel support around each slice line (meV)
cfg.ref_E0 = [815 792 756 729 713; 1520 1580 1660 1670 1594];   % reference L: E0 at the E0 nodes (P1b KB, meV)
cfg.ref_G = [360 510 580; 600 1400 1775];                         % reference L: widths at the width nodes (meV)
cfg.ref_fit = [300 1800];          % window of the residue fit that sets the reference amplitudes
cfg.regions = [250 1800; 150 250; 40 150];
cfg.old_sigma = 4.12e-4;           % old W: centred Gaussian probe (stage 1 sigma_probe), single slice
cfg.conv = struct('n_a', 8, 'n_s', 17, 'n_gh', 61, 'Qm_factor', 2, 'fine_step', 1, 'support', [40 100], 'tail_max', 0.08);
cfg.maps_q = [-0.002 -0.0005 0.0005 0.002];
cfg.maps_E = [300 600 1000 1500];
end
