function variants = thesis_sensitivity_configs(cfg)
%THESIS_SENSITIVITY_CONFIGS Small audit matrix around the thesis baseline.
%   These variants are intentionally few.  They test whether branch existence,
%   energy separation, and qualitative dispersion hierarchy survive the most
%   important analysis choices: deconvolution, background model, and prominence.

if nargin < 1 || isempty(cfg)
    cfg = thesis_config();
end

variants = repmat(struct('name', '', 'description', '', 'cfg', cfg), 1, 5);

variants(1).name = 'baseline';
variants(1).description = 'Baseline thesis configuration.';
variants(1).cfg = cfg;

variants(2).name = 'deconv_iter5';
variants(2).description = 'Same as baseline but enables Lucy-Richardson deconvolution with 5 iterations.';
variants(2).cfg = cfg;
variants(2).cfg.preprocess.do_deconv = true;
variants(2).cfg.preprocess.deconv_iter = 5;

variants(3).name = 'bg_auto';
variants(3).description = 'Auto background model selection with the same primary window.';
variants(3).cfg = cfg;
variants(3).cfg.preprocess.bg_method = 'Auto';
variants(3).cfg.preprocess.bg_candidate_methods = {'Power', 'Exp2', 'ExpPoly3', 'Pearson', 'PearsonVII'};

variants(4).name = 'prominence_strict';
variants(4).description = 'Stricter relative prominence threshold for peak detection.';
variants(4).cfg = cfg;
variants(4).cfg.fit.prominence = 0.15;

variants(5).name = 'bg_pearsonvii';
variants(5).description = 'PearsonVII background model in the same low-energy fit window.';
variants(5).cfg = cfg;
variants(5).cfg.preprocess.bg_method = 'PearsonVII';
end
