function varargout = qe_lowq_joint_forward(mode, varargin)
%QE_LOWQ_JOINT_FORWARD  Forward counts of the shared low-q response through sliced acceptance operators.
%
%   cache = qe_lowq_joint_forward('build', ops, R, options)
%       ops{b}.W   n_exit x n_slice x M radial weights (qe_acceptance_operator), same Qm
%       R{b}       n_slice x n_eps energy kernels of the slices (sampled at the eps grid,
%                  multiplied by the energy step, so each row sums to about 1)
%       options:   omega (loss grid, uniform), eps (kernel offset grid, same step),
%                  variant 'slice' | 'integrated' | 'exit', exit_slice (slice index
%                  of every exit channel, for 'exit'), slice_weight{b} (for 'integrated')
%   cols = qe_lowq_joint_forward('eval', cache, theta, intr)
%       forward of every linear column (branch j, residue basis k) on the output grid
%       cache.E_out: n_E x n_col x n_exit x n_block. theta.E0sq (2 x 5, meV^2),
%       theta.G (2 x 3, meV) at the nodes of intr (see 'intrinsic').
%   C = qe_lowq_joint_forward('C', Qm, omega, theta, intr)   M x n_omega x n_col
%   intr = qe_lowq_joint_forward('intrinsic', name, value, ...)  defaults below.
%
%   Model (P2a, review rounds 3-4): two Drude-Lorentz branches
%       L_j(Q, w) = A_j(Q) w G_j / ((w^2 - E0_j^2)^2 + w^2 G_j^2),
%       A_j(Q) = (Q / Q_ref) sum_k a_jk B_k(Q)   (a_jk >= 0, linear; 2D optical limit),
%   log E0_j^2 and log G_j linear in log Q between the nodes, constant outside;
%   below the first node E0^2 is constant ('const') or proportional to Q ('gapless').
%   B_k are hats in log Q on the residue nodes, flat outside. The thermal factor acts
%   before the convolution: S(w) = (L/w) * w / (1 - exp(-w/kT)) (expm1; limit kT at w = 0).
%   C(Q, w) = K(Q, |w|) S(Q, w), K = Q / (Q^2 + qE^2)^2, qE = |w| / (hbar v). At w = 0 the
%   loss energy in K and L/w is taken as step/4 (cell average); this affects only |E| < step.
%   Counts: mu_ib(E) = sum_r [R_br * (sum_m W_irm C(Q_m, .))](E) (variant 'slice'); the
%   kernels carry the energy step, W carries the slice weights.
%
%   See also: qe_acceptance_operator, qe_elastic_line_model, qe_kinematic_prefactor

switch mode
    case 'intrinsic'
        varargout{1} = local_intrinsic(varargin{:});
    case 'C'
        varargout{1} = local_C(varargin{:});
    case 'build'
        varargout{1} = local_build(varargin{:});
    case 'eval'
        varargout{1} = local_eval(varargin{:});
    otherwise
        error('qe_lowq_joint_forward:Mode', 'Unknown mode %s', mode);
end
end

function intr = local_intrinsic(varargin)
gam = 1 + 30e3 / 510998.95; v = sqrt(1 - 1 / gam ^ 2) * 299792458;
intr = struct('Qn_E0', [0.00025 0.00075 0.0015 0.003 0.005], 'Qn_G', [0.00025 0.0015 0.005], ...
    'Qn_A', [0.00025 0.0015 0.005], 'Q_ref', 0.001, 'ext', 'const', 'kT', 25.852, ...
    'qE_per_meV', 1e-3 / 6.582119569e-16 / v * 1e-10, 'n_branch', 2);
for k = 1:2:numel(varargin), intr.(varargin{k}) = varargin{k + 1}; end
end

function C = local_C(Qm, omega, theta, intr)
Q = Qm(:); w = omega(:).'; dw = abs(w(2) - w(1)); aw = max(abs(w), dw / 4);
qE = aw * intr.qE_per_meV; K = Q ./ (Q .^ 2 + qE .^ 2) .^ 2; K(Q == 0, :) = 0;
g = w ./ (-expm1(-w / intr.kT)); g(w == 0) = intr.kT;
nk = numel(intr.Qn_A); C = zeros(numel(Q), numel(w), intr.n_branch * nk);
lQ = log(max(Q, realmin));
B = local_hats(lQ, log(intr.Qn_A)) .* (Q / intr.Q_ref);
for j = 1:intr.n_branch
    lE = local_interp(lQ, log(intr.Qn_E0), log(theta.E0sq(j, :)));
    if strcmp(intr.ext, 'gapless')
        lo = Q < intr.Qn_E0(1); lE(lo) = log(theta.E0sq(j, 1)) + lQ(lo) - log(intr.Qn_E0(1));
    end
    E0sq = exp(lE); G = exp(local_interp(lQ, log(intr.Qn_G), log(theta.G(j, :))));
    Lw = G ./ ((aw .^ 2 - E0sq) .^ 2 + aw .^ 2 .* G .^ 2);      % L / w at |w|
    core = K .* Lw .* g;
    for k = 1:nk, C(:, :, (j - 1) * nk + k) = core .* B(:, k); end
end
end

function y = local_interp(x, xn, yn)
% piecewise linear in x through the nodes, constant outside
y = interp1(xn(:), yn(:), min(max(x, xn(1)), xn(end)), 'linear');
end

function B = local_hats(x, xn)
nk = numel(xn); B = zeros(numel(x), nk); xc = min(max(x, xn(1)), xn(end));
for k = 1:nk
    e = zeros(1, nk); e(k) = 1; B(:, k) = interp1(xn(:), e(:), xc, 'linear');
end
end

function cache = local_build(ops, R, options)
nb = numel(ops); [ni, ~, M] = size(ops{1}.W);
w = options.omega(:).'; e = options.eps(:).'; step = w(2) - w(1);
assert(abs(e(2) - e(1) - step) < 1e-9, 'qe_lowq_joint_forward:Grid', 'omega and eps need the same step');
nout = numel(w) + numel(e) - 1; N = 2 ^ nextpow2(nout); Nf = N / 2 + 1;
% double precision: C spans many decades between omega -> 0 and the B1 region, and the
% FFT spreads rounding errors of the largest values over all energies
H = complex(zeros(ni * nb, M, Nf));
for b = 1:nb
    Rb = R{b}; Rh = fft(Rb, N, 2); Rh = Rh(:, 1:Nf);                 % n_slice x Nf
    for i = 1:ni
        Wi = reshape(ops{b}.W(i, :, :), size(ops{b}.W, 2), M);          % n_slice x M
        switch options.variant
            case 'slice'
                Hi = Wi.' * Rh;                                          % M x Nf
            case 'integrated'
                pw = options.slice_weight{b}(:); Ri = (pw.' * Rh) / sum(pw);
                Hi = sum(Wi, 1).' * Ri;
            case 'exit'
                Hi = sum(Wi, 1).' * Rh(options.exit_slice(i), :);
            otherwise
                error('qe_lowq_joint_forward:Variant', 'Unknown variant %s', options.variant);
        end
        H((b - 1) * ni + i, :, :) = reshape(Hi, [1 M Nf]);
    end
end
cache = struct('H', H, 'Qm', ops{1}.Qm, 'omega', w, 'eps', e, 'N', N, 'Nf', Nf, 'n_exit', ni, 'n_block', nb, ...
    'E_out', w(1) + e(1) + step * (0:nout - 1), 'variant', options.variant);
end

function cols = local_eval(cache, theta, intr)
C = local_C(cache.Qm, cache.omega, theta, intr);                        % M x n_w x n_col
ncol = size(C, 3); M = numel(cache.Qm); N = cache.N; Nf = cache.Nf;
Ch = fft(C, N, 2); Ch = permute(Ch(:, 1:Nf, :), [1 3 2]);                % M x n_col x Nf
X = pagemtimes(cache.H, Ch);                                            % (n_exit n_block) x n_col x Nf
X = cat(3, X, conj(X(:, :, Nf - 1:-1:2)));
y = real(ifft(X, N, 3)); nout = numel(cache.E_out);
cols = permute(reshape(y(:, :, 1:nout), [cache.n_exit cache.n_block ncol nout]), [4 3 1 2]);
end
