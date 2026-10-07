function cfg = b1_fullq_config()
% Shared settings of the B1 full-q task. Values marked "stage 0" were fixed
% from the six-target check (see the run report).
cfg.q_grid = 0.0025:0.0015:0.0595;   % |q| bin centres (A^-1), N=3, non-overlapping
cfg.reference_abs_q = 0.001;          % ZLP alignment reference, as v7
cfg.zlp_window = [-100 100];
cfg.window = [300 1800];              % B1 signal window (meV)
cfg.beam_kV = 30;
cfg.models = {'lorentz', 'lorentz_symmetric'};   % DL (main), symmetric Lorentzian (cross-check)
cfg.n_zlp = 2;                        % stage 0
cfg.aux = {'aux2f', [30 80; 80 300]; 'aux1', [30 300]};   % first row = main
cfg.kinematic = {'3d', '2d'};         % first = main (stage 0: 2D degenerates at |q| = 0.0025)
cfg.h_main = 0.002;                   % perpendicular half-acceptance (A^-1), stage 0
cfg.h_variants = [0.0005 0.008];
cfg.prefactor_floor = 250;            % meV, stage 0 (no low-energy amplification of B1 tails)
cfg.prefactor_on_aux = false;
cfg.h_variant_qmax = 0.0065;          % h variants only where they matter
cfg.n_starts = 12;                    % main configuration
cfg.n_starts_other = 8;               % all other configurations (plus warm starts)
cfg.seed = 20261007;
cfg.n_boot = 30;                      % frame bootstrap replicas, main configuration
cfg.boot_starts = 3;
cfg.n_synth = 20;                     % stage 3 replicas per selected bin
end
