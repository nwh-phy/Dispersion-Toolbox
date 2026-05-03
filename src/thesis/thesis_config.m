function cfg = thesis_config()
%THESIS_CONFIG  Conservative reproducible q-EELS pipeline configuration.
%   The thesis pipeline is intentionally separate from exploratory GUI defaults
%   and legacy scripts.  All Chapter 3 numerical results should be generated
%   from this configuration or from explicitly named sensitivity variants.

root_dir = fileparts(fileparts(fileparts(mfilename('fullpath'))));

cfg = struct();
cfg.pipeline = struct( ...
    'name', 'thesis_qeels_pipeline', ...
    'version', '0.1.0', ...
    'description', ['Conservative BiSb q-EELS thesis pipeline: ZLP peak ', ...
                    'normalization, Wiener2D denoising, quasi-elastic ', ...
                    'background subtraction, pre-subtracted Drude-Lorentz ', ...
                    'fitting, deterministic branch assignment.']);

cfg.project_root = root_dir;
cfg.analysis = struct( ...
    'q_range_Ainv', [-0.15 0.15], ...
    'E_range_meV', [300 3836], ...
    'q_skip_Ainv', 0.005);

cfg.preprocess = struct();
cfg.preprocess.do_despike = false;
cfg.preprocess.do_normalize = true;
cfg.preprocess.norm_method = 'ZLP Peak';
cfg.preprocess.norm_min = -50;
cfg.preprocess.norm_max = 50;
cfg.preprocess.do_denoise = true;
cfg.preprocess.denoise_method = 'Wiener2D';
cfg.preprocess.denoise_sigma = 0;
cfg.preprocess.sg_order = 3;
cfg.preprocess.sg_framelen = 11;
cfg.preprocess.do_bg_sub = true;
cfg.preprocess.bg_method = 'Power';
cfg.preprocess.bg_win_lo = [50 300];
cfg.preprocess.bg_win_hi = [];
cfg.preprocess.bg_iterative = false;
cfg.preprocess.do_deconv = false;
cfg.preprocess.deconv_iter = 5;

cfg.fit = struct();
cfg.fit.E_min = cfg.analysis.E_range_meV(1);
cfg.fit.E_max = cfg.analysis.E_range_meV(2);
cfg.fit.q_start = cfg.analysis.q_range_Ainv(1);
cfg.fit.q_end = cfg.analysis.q_range_Ainv(2);
cfg.fit.prominence = 0.05;  % relative prominence, multiplied by per-spectrum max
cfg.fit.smooth_width = 25;
cfg.fit.max_peaks = 2;
cfg.fit.peak_model = 'lorentz';
cfg.fit.pre_subtracted = true;
cfg.fit.guesses = [];
cfg.fit.seed_idx = [];
cfg.fit.max_shift = 350;
cfg.fit.R2_threshold = 0.3;

cfg.branch = struct();
cfg.branch.q_skip_Ainv = cfg.analysis.q_skip_Ainv;
cfg.branch.max_gamma_ratio = 2.0;
cfg.branch.min_R2 = cfg.fit.R2_threshold;
cfg.branch.score_column = 12;  % prefer raw peak height when available
cfg.branch.low = struct( ...
    'name', 'low', ...
    'energy_window_meV', [500 2100], ...
    'max_delta_meV', 350, ...
    'max_q_gap_Ainv', 0.015);
cfg.branch.high = struct( ...
    'name', 'high', ...
    'energy_window_meV', [2800 3800], ...
    'max_delta_meV', 350, ...
    'max_q_gap_Ainv', 0.015);

cfg.symmetry = struct();
cfg.symmetry.q_tolerance_Ainv = 1e-9;
cfg.symmetry.min_weight = 1e-9;
cfg.symmetry.weight_mode = 'r2_damping_height';

cfg.models = struct();
cfg.models.low = 'quasi2d_plasmon';
cfg.models.high = 'optical_constant';
cfg.models.epsilon_s = 1;

cfg.output = struct();
cfg.output.root = fullfile(root_dir, 'paper_results', 'thesis_qeels_pipeline');
cfg.output.baseline_dir = fullfile(cfg.output.root, 'baseline');
cfg.output.sensitivity_dir = fullfile(cfg.output.root, 'sensitivity');

cfg.notes = struct();
cfg.notes.sample_label = 'encapsulated two-dimensional BiSb alloy sample';
cfg.notes.naming_caution = ['Some data folders retain the historical "Bi" label; ', ...
    'the thesis interpretation should use BiSb alloy wording when supported by sample context.'];
cfg.notes.physics_caution = ['rho0 and E_flat are effective fit parameters under ', ...
    'the selected empirical model and dielectric-background assumption, not direct ', ...
    'intrinsic BiSb constants.'];
end
