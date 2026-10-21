function out = b1_p2b_hscan_run(h, run_root)
% One point of the slot-acceptance scan of P2b: the unchanged fit code b1_p2b_fit with
% h (A^-1) and the extended node grid of b1_p2b_slot_config('hscan', h).
% Run one MATLAB per h; compare chi2 in the window and in the 150-250 meV holdout.
root = b1_simple_root();
if nargin < 2, run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
[cfg, ocfg] = b1_p2b_slot_config('hscan', h);
out = b1_p2b_fit(run_root, cfg, ocfg);
fprintf('h scan point h = %.4f written to %s\n', h, out);
end
