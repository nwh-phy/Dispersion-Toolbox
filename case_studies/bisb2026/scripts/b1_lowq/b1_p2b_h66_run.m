function out = b1_p2b_h66_run(run_root,cfg,ocfg)
% P2b rerun with the slot acceptance of the 0.125 mm slot (h = 0.0066 A^-1 at 10w instead
% of the assumed 0.0005): the unchanged fit code b1_p2b_fit with the configurations
% b1_p2b_slot_config('h66') (h, n_sy, version tag, three variants).
if nargin < 1, run_root = []; end
if nargin < 2, cfg = []; end
if nargin < 3, ocfg = []; end
[cfg, ocfg] = b1_p2b_slot_config('h66', [], cfg, ocfg);
out = b1_p2b_fit(run_root, cfg, ocfg);
end
