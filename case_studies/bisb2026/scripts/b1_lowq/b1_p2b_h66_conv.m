function T = b1_p2b_h66_conv(run_root, cfg, ocfg0, grid_s)
% Forward check for the h = 0.0066 rerun of P2b: the P2b main operator (slice kernels,
% broad component, y mixture) evaluated at the P2b main L (nodes and residues fixed) for
%   h 0.0005 / n_sy 9 (P2b as run), h 0.0066 / n_sy 9, 53, 105.
% Reports per exit channel, in the B1 window 250-1800 meV and summed over the six blocks,
%   conv:   max |mu(n_sy) - mu(105)| / sqrt(mu(105))                (quadrature)
%   effect: max |shape(h 0.0005) - shape(h 0.0066, 105)| / sqrt(mu)  (shape = counts scaled to
%           the same window sum, i.e. what a residue refit cannot absorb)
% Output: <run_root>/p2b/h66_conv/.
if nargin < 1, run_root = []; end
if nargin < 2 || isempty(cfg), cfg = b1_p2b_config(); end
paths = b1_forward_paths(cfg, run_root);
if nargin < 3 || isempty(ocfg0)
    [~, ocfg0] = b1_p2b_slot_config('h66', [], cfg);
end
if nargin < 4, grid_s = [0.0005 9; ocfg0.h 9; ocfg0.h ocfg0.n_sy; ocfg0.h 2*ocfg0.n_sy-1]; end
t0 = tic;
out = fullfile(paths.p2b_dir, 'h66_conv');
if isfield(cfg,'output_dir'), out = cfg.output_dir; end
S = load(fullfile(paths.fit_dir, 'p2b_results.mat'), 'res', 'cfg'); F = S.res(1).fit;
cfg.nodes = S.cfg.nodes; cfg.exit_q = S.cfg.exit_q;
ind = paths.input_dir; rd = paths.response_dir;
L = load(fullfile(ind, 'p2a_cube.mat'), 'cube'); cube = L.cube; q = cube.q(:).'; dq = cube.dq;
RR = load(fullfile(rd, 'response.mat'), 'R', 'J6');
ie = arrayfun(@(v) find(abs(q - v) < 1e-9), cfg.exit_q); ni = numel(ie);
intr = qe_lowq_joint_forward('intrinsic', 'Qn_E0', cfg.nodes.E0, 'Qn_G', cfg.nodes.G, 'Qn_A', cfg.nodes.A);
th = b1_p2b_decode(F.u, cfg.nodes);
mu = cell(size(grid_s, 1), 1);
for s = 1:size(grid_s, 1)
    oc = ocfg0; oc.h = grid_s(s, 1); oc.n_sy = grid_s(s, 2); oc.exit_q = cfg.exit_q;
    st = struct('blocking', 'b6', 'variant', 'slice', 'ytype', 'mixture', 'wide', 'on', 'ext', 'const', 'n_a', oc.n_a, ...
        'n_sx', oc.n_sx, 'n_sy', oc.n_sy, 'n_gh', oc.n_gh, 'Qf', 1, 'step', oc.omega(2), 'support', oc.support, ...
        'tail_max', oc.tail_edges(end), 'exits', 1:ni, 'blocks', []);
    [ops, K, sw, ex] = b1_p2a_build_ops(st, oc, cube, RR.R, RR.J6, q, dq, ie);
    fo = struct('omega', oc.omega(1):st.step:oc.omega(3), 'eps', -72:st.step:72, 'variant', 'slice', ...
        'slice_weight', {sw}, 'exit_slice', ex);
    cache = qe_lowq_joint_forward('build', ops, K, fo); clear ops K
    cols = qe_lowq_joint_forward('eval', cache, th, intr);
    mu{s} = squeeze(sum(sum(cols .* reshape(F.a, 1, [], 1, 1), 2), 4));      % n_E x n_exit, blocks summed
    Eo = cache.E_out(:); clear cache cols
    fprintf('h %.4f n_sy %3d built and evaluated (%.0f s)\n', grid_s(s, 1), grid_s(s, 2), toc(t0));
end
window=S.cfg.variants(1).window; m = Eo >= window(1) & Eo <= window(2); ref = mu{end}(m, :); sref = sqrt(max(ref, 1));
shape = @(x) x(m, :) ./ sum(x(m, :), 1) .* sum(ref, 1);
T = table(cfg.exit_q(:), 'VariableNames', {'q'});
T.conv_nsy9 = max(abs(mu{2}(m, :) - ref) ./ sref, [], 1).';
T.conv_nsy53 = max(abs(mu{3}(m, :) - ref) ./ sref, [], 1).';
T.effect_h = max(abs(shape(mu{1}) - ref) ./ sref, [], 1).';
T.acceptance_ratio = (sum(mu{4}(m, :), 1) ./ sum(mu{1}(m, :), 1)).';       % window counts, h 0.0066 / 0.0005 at fixed residues
if ~isfolder(out), mkdir(out); end
writetable(T, fullfile(out, 'h66_conv.csv')); disp(T);
fprintf('h66 forward check written to %s (%.0f s)\n', out, toc(t0));
end
