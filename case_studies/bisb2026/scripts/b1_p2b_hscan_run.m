function out = b1_p2b_hscan_run(h, run_root)
% One point of the slot-acceptance scan of P2b: the unchanged main-checkout b1_p2b_fit
% with h (A^-1) and the extended node grid of ./hscan (configs read h from B1_P2B_H).
% Run one MATLAB per h; compare chi2 in the window and in the 150-250 meV holdout.
here = fileparts(mfilename('fullpath')); root = b1_simple_root(); addpath(fullfile(here, 'hscan'));
assert(contains(which('b1_p2a_operator_config'), [filesep 'hscan' filesep]), 'hscan override not first on the path');
assert(contains(which('b1_p2b_config'), [filesep 'hscan' filesep]), 'hscan override not first on the path');
if nargin < 2, run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
setenv('B1_P2B_H', sprintf('%.6g', h));
out = b1_p2b_fit(run_root);
fprintf('h scan point h = %.4f written to %s\n', h, out);
end
