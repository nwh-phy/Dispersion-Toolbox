function out = b1_p2b_h66_run(run_root)
% P2b rerun with the slot acceptance of the 0.125 mm slot (h = 0.0066 A^-1 at 10w instead
% of the assumed 0.0005). The fit code is the unchanged main-checkout b1_p2b_fit; only
% b1_p2a_operator_config (h, n_sy) and b1_p2b_config (version tag, three variants) are
% overridden by the copies in ./h66, which this function puts first on the path.
here = fileparts(mfilename('fullpath')); root = b1_simple_root(); addpath(fullfile(here, 'h66'));
assert(contains(which('b1_p2a_operator_config'), [filesep 'h66' filesep]), 'h66 override not first on the path');
assert(contains(which('b1_p2b_config'), [filesep 'h66' filesep]), 'h66 override not first on the path');
if nargin < 1, run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
out = b1_p2b_fit(run_root);
end
