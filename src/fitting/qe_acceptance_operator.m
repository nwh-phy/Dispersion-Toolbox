function op = qe_acceptance_operator(q_exit, ax, wx, rlab, ay, wy, Qm, options)
%QE_ACCEPTANCE_OPERATOR  Momentum-transfer weights of exit channels for sliced incident beams.
%
%   op = qe_acceptance_operator(q_exit, ax, wx, rlab, ay, wy, Qm, dq=0.0005, h=0.0005)
%
%   For exit channel i (box s_x in [q_i - dq/2, q_i + dq/2], s_y in [-h, h]) and
%   incident slice r (x nodes ax with weights wx and slice labels rlab; y nodes
%   ay with weights wy), the integral
%       G_ir(omega) = int_{box i} d^2s  sum_{n in r} wx_n sum_k wy_nk C(|s - a_nk|, omega)
%   is written as G_ir(omega) = sum_m W(i, r, m) C(Q_m, omega) on the radial grid Qm
%   (linear interpolation in Q between the nodes, so a C that is linear in Q
%   between nodes is integrated exactly). The exit box is sampled with n_sx x n_sy
%   midpoints and carries its area dq * 2h.
%
%   wy is a column (the same y distribution for every x node: separable beam) or a
%   matrix with one row per x node (a non-separable beam, e.g. a narrow and a broad
%   component with their own y widths). Weights are not normalised here: wx carries
%   the slice weights (e.g. elastic counts of the slice), so sum_m W(i, r, m) is
%   box area x (slice weight) x (y weight sum).
%
%   op.W       n_exit x n_slice x numel(Qm) weights
%   op.Qm      radial grid (A^-1), op.out_of_range  weight beyond Qm(end) (added to the last node)
%
%   See also: qe_lowq_joint_forward, qe_kinematic_prefactor

arguments
    q_exit (:,1) double
    ax (:,1) double
    wx (:,1) double
    rlab (:,1) double {mustBeInteger, mustBePositive}
    ay (:,1) double
    wy double
    Qm (:,1) double
    options.dq (1,1) double {mustBePositive} = 0.0005
    options.h (1,1) double {mustBeNonnegative} = 0.0005
    options.n_sx (1,1) double {mustBePositive, mustBeInteger} = 9
    options.n_sy (1,1) double {mustBePositive, mustBeInteger} = 9
end
assert(numel(ax) == numel(wx) && numel(ax) == numel(rlab), 'qe_acceptance_operator:Size', 'ax, wx, rlab differ in size');
assert(Qm(1) == 0 && all(diff(Qm) > 0), 'qe_acceptance_operator:Grid', 'Qm must start at 0 and increase');
if isvector(wy), wy = repmat(wy(:).', numel(ax), 1); end
assert(isequal(size(wy), [numel(ax) numel(ay)]), 'qe_acceptance_operator:Wy', 'wy must be numel(ay) or numel(ax) x numel(ay)');

nr = max(rlab); M = numel(Qm); ni = numel(q_exit);
bx = ((1:options.n_sx) - 0.5) / options.n_sx - 0.5; by = ((1:options.n_sy) - 0.5) / options.n_sy - 0.5;
sy = 2 * options.h * by; if options.h == 0, sy = 0; end
area = options.dq * max(2 * options.h, eps) / (options.n_sx * numel(sy));
if options.h == 0, area = options.dq / options.n_sx; end
W = zeros(ni, nr, M); oor = 0;
w_xy = wx .* wy;                                     % n_x x n_y
for i = 1:ni
    sx = q_exit(i) + options.dq * bx;                % 1 x n_sx
    dx = sx(:) - ax(:).';                            % n_sx x n_x
    dy = sy(:) - ay(:).';                            % n_sy x n_y
    Q = sqrt(reshape(dx, [numel(sx) 1 numel(ax) 1]) .^ 2 + reshape(dy, [1 numel(sy) 1 numel(ay)]) .^ 2);
    w = repmat(reshape(w_xy, [1 1 numel(ax) numel(ay)]), numel(sx), numel(sy)) * area;
    r = repmat(reshape(rlab, [1 1 numel(ax) 1]), numel(sx), numel(sy), 1, numel(ay));
    Q = Q(:); w = w(:); r = r(:);
    hi = Q >= Qm(end); oor = oor + sum(w(hi)); Q(hi) = Qm(end) * (1 - 1e-12);
    j = discretize(Q, Qm); t = (Q - Qm(j)) ./ (Qm(j + 1) - Qm(j));
    W(i, :, :) = reshape(accumarray([r j; r j + 1], [w .* (1 - t); w .* t], [nr M]), [1 nr M]);
end
op = struct('W', W, 'Qm', Qm, 'q_exit', q_exit, 'dq', options.dq, 'h', options.h, 'n_sx', options.n_sx, ...
    'n_sy', options.n_sy, 'out_of_range', oor);
end
