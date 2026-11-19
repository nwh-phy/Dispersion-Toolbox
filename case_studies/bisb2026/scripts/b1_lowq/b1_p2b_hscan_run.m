function out = b1_p2b_hscan_run(h, run_root, cfg, ocfg)
% One point of the slot-acceptance scan of P2b: the unchanged fit code b1_p2b_fit with
% h (A^-1) and the extended node grid of b1_p2b_slot_config('hscan', h).
% Run one MATLAB per h; compare chi2 in the window and in the 150-250 meV holdout.
if nargin < 2, run_root = []; end
if nargin < 3, cfg = []; end
if nargin < 4, ocfg = []; end
[cfg, ocfg] = b1_p2b_slot_config('hscan', h, cfg, ocfg);
out = b1_p2b_fit(run_root, cfg, ocfg);
fprintf('h scan point h = %.4f written to %s\n', h, out);
end
