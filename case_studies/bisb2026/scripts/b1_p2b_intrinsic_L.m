function [L, Lj] = b1_p2b_intrinsic_L(fit_dir, Q, E, variant)
%B1_P2B_INTRINSIC_L  Intrinsic loss function of a stored P2b fit.
%   [L, Lj] = b1_p2b_intrinsic_L(fit_dir, Q, E) evaluates the fitted two-branch
%   Drude-Lorentz loss function L(Q,E) of the main variant in
%   fit_dir/p2b_results.mat, without the kinematic factor, the thermal factor and
%   the measurement operator. Q (A^-1) and E (meV, uniform step, E > 0) are vectors;
%   L is numel(Q) x numel(E), Lj is numel(Q) x numel(E) x 2 (low, high branch).
%   Same nodes, extension and residues as the fit: qe_lowq_joint_forward('C', ...)
%   divided by the kinematic and thermal factors that 'C' multiplies in.
%   variant selects res(variant) (default 1 = main).
arguments
    fit_dir char
    Q (:, 1) double {mustBePositive}
    E (1, :) double {mustBePositive}
    variant (1, 1) double {mustBePositive, mustBeInteger} = 1
end
S = load(fullfile(fit_dir, 'p2b_results.mat'), 'cfg', 'res'); cfg = S.cfg; R = S.res(variant); F = R.fit;
V = cfg.variants(strcmp({cfg.variants.name}, R.name));
intr = qe_lowq_joint_forward('intrinsic', 'Qn_E0', cfg.nodes.E0, 'Qn_G', cfg.nodes.G, 'Qn_A', cfg.nodes.A, 'ext', V.ext);
theta = struct('E0sq', exp(reshape(F.u(1:10), 5, 2).'), 'G', exp(reshape(F.u(11:16), 3, 2).'));
C = qe_lowq_joint_forward('C', Q, E, theta, intr);                     % numel(Q) x numel(E) x n_col
qE = E * intr.qE_per_meV; K = Q ./ (Q .^ 2 + qE .^ 2) .^ 2;          % as in 'C' (E > step / 4)
g = E ./ (-expm1(-E / intr.kT));
undo = E ./ (K .* g);                                                   % C = K .* (L / E) .* g per residue column
nk = numel(intr.Qn_A); Lj = zeros(numel(Q), numel(E), intr.n_branch);
for j = 1:intr.n_branch
    for k = 1:nk, Lj(:, :, j) = Lj(:, :, j) + F.a((j - 1) * nk + k) * C(:, :, (j - 1) * nk + k); end
    Lj(:, :, j) = Lj(:, :, j) .* undo;
end
L = sum(Lj, 3);
end
