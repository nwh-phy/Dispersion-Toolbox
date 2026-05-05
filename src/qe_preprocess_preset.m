function opts = qe_preprocess_preset(name, overrides)
%QE_PREPROCESS_PRESET  Named preprocessing option bundles.
%   OPTS = QE_PREPROCESS_PRESET(NAME) returns a struct accepted by
%   qe_preprocess.  Named presets keep GUI, scripts, and thesis pipelines from
%   drifting when a data contract has been validated.

arguments
    name {mustBeTextScalar} = "bisb_lowq_stable_v1"
    overrides struct = struct()
end

preset_name = lower(strtrim(char(name)));
switch preset_name
    case 'bisb_lowq_stable_v1'
        opts = local_base_bisb_opts();
        opts.bg_method = 'Auto';
        opts.bg_candidate_methods = {'Power', 'Exp2', 'ExpPoly3', 'Pearson', 'PearsonVII'};
        opts.bg_win_hi = [2800 3400];
        opts.bg_auto_group_qmax = 0.01;
        opts.bg_auto_group_ranges = zeros(0, 2);

    case 'thesis_baseline_v1'
        opts = local_base_bisb_opts();
        opts.bg_method = 'Power';
        opts.bg_candidate_methods = {};
        opts.bg_win_hi = [];

    otherwise
        error('qe_preprocess_preset:UnknownPreset', ...
            'Unknown preprocessing preset "%s".', name);
end

opts = local_apply_overrides(opts, overrides);
end


function opts = local_base_bisb_opts()
opts = struct();
opts.do_despike = false;
opts.do_normalize = true;
opts.norm_method = 'ZLP Peak';
opts.norm_min = -50;
opts.norm_max = 50;
opts.do_denoise = true;
opts.denoise_method = 'Wiener2D';
opts.denoise_sigma = 0;
opts.sg_order = 3;
opts.sg_framelen = 11;
opts.do_bg_sub = true;
opts.bg_method = 'Power';
opts.bg_candidate_methods = {};
opts.bg_win_lo = [50 300];
opts.bg_win_hi = [];
opts.bg_iterative = false;
opts.bg_auto_group_ranges = zeros(0, 2);
opts.do_deconv = false;
opts.deconv_iter = 5;
end


function opts = local_apply_overrides(opts, overrides)
fields = fieldnames(overrides);
for i = 1:numel(fields)
    opts.(fields{i}) = overrides.(fields{i});
end
end
