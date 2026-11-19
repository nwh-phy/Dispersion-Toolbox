function cfg = b1_forward_20w_config()
%B1_FORWARD_20W_CONFIG  S1–S3 inputs and conventions for the two-layer run.
cfg.version = 'P1-2026-11-19';
cfg.paths = b1_forward_config('no_PL2_20w_2film');
cfg.bad_column = 461;
cfg.nodes = b1_p2b_config().nodes;
cfg.nodes.E0 = sort([cfg.nodes.E0 0.0015 0.0025]);
cfg.variants = {'main','y: Gaussian','kernel: exit', ...
    'broad component off for B1','low-Q extension: gapless'};
cfg.h_scan = [0.0015 0.0025 0.005];
cfg.q_min = 0.001;
cfg.layer_qmin = [0.002 0.005];
cfg.baseline_fits = {'fit_590_3fb982f2','fit_590_691c36ad'};
end
