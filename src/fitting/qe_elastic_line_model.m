function varargout = qe_elastic_line_model(mode, varargin)
%QE_ELASTIC_LINE_MODEL  Elastic line of one angular slice: two asymmetric Pearson VII
%   components sharing one centre (the ZLP form of qe_zlp_joint_fit).
%
%   y = qe_elastic_line_model('eval', E, th)        counts at the energies E (meV)
%   [th, info] = qe_elastic_line_model('fit', E, Y, window, opts)
%       fits the counts Y (column, at E) inside [c0 + window(1), c0 + window(2)]
%       (c0 = log-parabola centre of the maximum within |E| < 40 meV), Poisson
%       weights; th = struct with c, A(2), sL(2), sR(2), mL(2), mR(2). Constraints:
%       m >= 2, 2 <= s <= 40 (60) meV for the first (second) component, sR/sL in [1/3, 3].
%   k = qe_elastic_line_model('kernel', x, th, support)
%       unit-area kernel at the offsets x (meV, absolute energy - 0), truncated to
%       |E - c| <= support and normalised on a 0.25 meV grid; returns zero outside.
%   a = qe_elastic_line_model('area', th, lim)      integral over c + [lim(1), lim(2)]
%
%   Component k: A_k (1 + x^2 / (m s^2))^(-m), s = sL, m = mL for x < 0 and sR, mR
%   otherwise, x = E - c (the parametrisation of local_pearson in qe_zlp_joint_fit).
%
%   See also: qe_zlp_joint_fit, b1_p2a_response

switch mode
    case 'eval'
        varargout{1} = local_eval(varargin{1}, varargin{2});
    case 'fit'
        [varargout{1}, varargout{2}] = local_fit(varargin{:});
    case 'kernel'
        x = varargin{1}; th = varargin{2}; S = varargin{3};
        xf = (th.c - S:0.25:th.c + S).'; nrm = trapz(xf, local_eval(xf, th));
        k = local_eval(x, th) / nrm; k(abs(x - th.c) > S) = 0;
        varargout{1} = k;
    case 'area'
        th = varargin{1}; lim = varargin{2};
        xf = (th.c + lim(1):0.25:th.c + lim(2)).';
        varargout{1} = trapz(xf, local_eval(xf, th));
    otherwise
        error('qe_elastic_line_model:Mode', 'Unknown mode %s', mode);
end
end

function y = local_eval(E, th)
x = E - th.c; y = zeros(size(x));
for k = 1:numel(th.A)
    y = y + th.A(k) * local_pearson(x, th.sL(k), th.sR(k), th.mL(k), th.mR(k));
end
end

function y = local_pearson(x, sL, sR, mL, mR)
y = zeros(size(x)); l = x < 0;
y(l) = (1 + (x(l) / sL) .^ 2 / mL) .^ (-mL);
y(~l) = (1 + (x(~l) / sR) .^ 2 / mR) .^ (-mR);
end

function [th, info] = local_fit(E, Y, window, opts)
E = E(:); Y = Y(:); dE = E(2) - E(1);
[~, i0] = max(Y .* (abs(E) < 40)); ly = log(max(Y(i0-1:i0+1), 1));
c0 = E(i0) + (ly(3) - ly(1)) / (2 * (2 * ly(2) - ly(1) - ly(3))) * dE;
m = E >= c0 + window(1) & E <= c0 + window(2); x = E(m); y = Y(m); w = 1 ./ sqrt(max(y, 1));
% light tails (m >= 2) and bounded left/right asymmetry: outside the window the
% model only extrapolates, and the instrument tails are weak (the gain side of the
% central channel is 2.5e-5 of the peak at -60..-150 meV in 590)
u0 = [c0 log([10 10 4 4 20 20 2.5 2.5])];
lb = [c0 - 8 log([2 2 2 2 5 5 2 2])]; ub = [c0 + 8 log([40 40 50 50 60 60 50 50])];
Ai = zeros(4, 9); Ai(1, [2 3]) = [-1 1]; Ai(2, [2 3]) = [1 -1]; Ai(3, [6 7]) = [-1 1]; Ai(4, [6 7]) = [1 -1];
bA = log(3) * ones(4, 1);          % |log(sR/sL)| <= log 3 for both components
f = @(u) local_res(u, x, y, w);
u = lsqnonlin(f, u0, lb, ub, Ai, bA, [], [], [], opts);
[r, A] = f(u); th = local_decode(u, A);
info = struct('c0', c0, 'chi2dof', sum(r .^ 2) / max(numel(y) - 11, 1), 'n', numel(y), 'window', window, ...
    'at_bound', any(abs(u - lb) < 1e-6 | abs(u - ub) < 1e-6));
end

function [r, A] = local_res(u, x, y, w)
th = local_decode(u, [1 1]); B = zeros(numel(x), 2);
for k = 1:2, B(:, k) = local_pearson(x - th.c, th.sL(k), th.sR(k), th.mL(k), th.mR(k)); end
A = lsqnonneg(B .* w, y .* w); r = (y - B * A) .* w;
end

function th = local_decode(u, A)
v = exp(u(2:end));
th = struct('c', u(1), 'A', A(:).', 'sL', v([1 5]), 'sR', v([2 6]), 'mL', v([3 7]), 'mR', v([4 8]));
end
