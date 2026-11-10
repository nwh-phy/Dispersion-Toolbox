function cfg = b1_adaptq_config()
% B1 adaptive-binning task: native channels (N=1) where the counts allow it
% (low q), wider non-overlapping bins only at high q to reduce noise.
% Fit settings are those of b1_fullq (stage 0 decisions) unless noted.
f = b1_fullq_config();
% Bin layout per side, as [lo, hi, N]: native channels with lo < |q| <= hi,
% grouped by N; the channel spacing comes from the data set (the worker sets
% cfg.dq from the channel file). Counts below are for dq = 0.0005.
cfg.layout = [0      0.0060 1;    % every native channel with 0 < |q| <= 0.006 (12)
              0.0060 0.0140 3;    % 5 bins, centres 0.0070 ... 0.0130
              0.0140 0.0240 5;    % 4 bins, centres 0.0155 ... 0.0230
              0.0240 0.0310 7];   % 2 bins, centres 0.0260, 0.0295
cfg.window = f.window; cfg.zlp_window = f.zlp_window; cfg.beam_kV = f.beam_kV;
cfg.n_zlp = f.n_zlp; cfg.aux = f.aux(1,:); cfg.prefactor_floor = f.prefactor_floor;
cfg.prefactor_on_aux = f.prefactor_on_aux; cfg.h_main = f.h_main; cfg.h_variants = f.h_variants;
cfg.h_variant_qmax = 0.0040;      % h matters most where q < h
% Configurations: DL n=2 / n=1, 3D (main) and 2D kinematic forms.
cfg.configs = {'lorentz',2,'3d'; 'lorentz',2,'2d'; 'lorentz',1,'3d'; 'lorentz',1,'2d'};
cfg.n_starts = 12; cfg.n_starts_other = 8;
cfg.seed = 20261009;
cfg.n_boot = 10; cfg.boot_starts = 3;
end
