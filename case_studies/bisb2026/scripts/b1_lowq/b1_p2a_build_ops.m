function varargout = b1_p2a_build_ops(varargin)
% P2a step 2: block operators of 590 from the slice response (shared by
% b1_p2a_operator and b1_p2a_synth).
%   [ops, K, sw, ex] = b1_p2a_build_ops(st, cfg, cube, R, J6, q, dq, ie)
%       st: setting struct (blocking, ytype, wide, n_a, n_sx, n_sy, n_gh, Qf, step,
%       support, tail_max, exits, blocks); returns per block the operator
%       (qe_acceptance_operator), the slice kernels (rows, times the energy step),
%       the slice weights and the exit-slice index of every exit channel.
%   Qm = b1_p2a_build_ops('Qm', spec, factor)    radial grid
%   [x, w] = b1_p2a_build_ops('gh', n)           Gauss-Hermite nodes
if ischar(varargin{1})
    switch varargin{1}
        case 'Qm', varargout{1} = local_Qm(varargin{2}, varargin{3});
        case 'gh', [varargout{1}, varargout{2}] = local_gh(varargin{2});
    end
    return
end
[varargout{1:nargout}] = local_build(varargin{:});
end

function Qm = local_Qm(spec, f)
Qm = [];
for k = 1:size(spec, 1), Qm = [Qm, spec(k, 1):spec(k, 2) / f:spec(k, 3)]; end %#ok<AGROW>
Qm = unique(round(Qm(:), 12));
end

function [x, w] = local_gl(n)
b = (1:n - 1) ./ sqrt(4 * (1:n - 1) .^ 2 - 1); [V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D)); w = 2 * V(1, i) .^ 2; w = w(:);
end

function [x, w] = local_gh(n)
b = sqrt((1:n - 1) / 2); [V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D)); w = sqrt(pi) * V(1, i) .^ 2; w = w(:);
end

function [ops, K, sw, ex] = local_build(st, cfg, cube, R, J6, q, dq, ie)
% incident nodes per block, y nodes, operators, slice kernels
Rb = R.(st.blocking); nbk = size(Rb.area, 1); blocks = st.blocks; if isempty(blocks), blocks = 1:nbk; end
if strcmp(st.blocking, 'b6')
    pos = J6.pos; logA = J6.logA;
else
    pos = cube.b12.pos_narrow_mode + mean(J6.pos) - mean(cube.b6.pos_narrow_mode); logA = [];
end
step = st.step; eps = -72:step:72; Qm = local_Qm(cfg.Qm_spec, st.Qf);
te = cfg.tail_edges(cfg.tail_edges <= st.tail_max + 1e-12); if te(end) < st.tail_max - 1e-12, te = [te st.tail_max]; end
[gx, gw] = local_gl(st.n_a); nk = numel(q); nexit = numel(st.exits);
ops = cell(1, numel(blocks)); K = ops; sw = ops; ex = zeros(nexit, 1);
for bb = 1:numel(blocks)
    b = blocks(bb);
    if isempty(logA)
        p0 = [0 pos(b) / dq J6.shape];
        s0 = sum(qe_angular_profile_model('channel', q(:), dq, p0)); p = [log(sum(Rb.area(b, :)) / s0) pos(b) / dq J6.shape];
    else
        p = [logA(b) pos(b) / dq J6.shape];
    end
    % slices: measured channels, then tails on both sides
    lo = [q - dq / 2, -fliplr(te(2:end)), te(1:end - 1)]; hi = [q + dq / 2, -fliplr(te(1:end - 1)), te(2:end)];
    kidx = [1:nk, ones(1, numel(te) - 1), nk * ones(1, numel(te) - 1)];
    ax = []; wx = []; rl = [];
    for s = 1:numel(lo)
        a = (lo(s) + hi(s)) / 2 + (hi(s) - lo(s)) / 2 * gx; cw = (hi(s) - lo(s)) / 2 * gw;
        if strcmp(st.wide, 'off')
            c = qe_angular_profile_model('components', a, dq, p); d = c(:, 1);
        else
            d = qe_angular_profile_model('density', a, dq, p);
        end
        wgt = cw .* d / dq;                                               % counts
        if s <= nk && strcmp(st.wide, 'on'), wgt = wgt / sum(wgt) * Rb.area(b, s); end
        ax = [ax; a]; wx = [wx; wgt]; rl = [rl; s * ones(numel(a), 1)]; %#ok<AGROW>
    end
    [ay, wy] = local_y(st, cfg, J6, dq, ax, p);
    op = qe_acceptance_operator(cfg.exit_q(st.exits).', ax, wx, rl, ay, wy, Qm, dq=dq, h=cfg.h, n_sx=st.n_sx, n_sy=st.n_sy);
    ops{bb} = op;
    Kb = zeros(numel(lo), numel(eps));
    for s = 1:numel(lo), Kb(s, :) = qe_elastic_line_model('kernel', eps, Rb.th{b, kidx(s)}, st.support) * step; end
    K{bb} = Kb; sw{bb} = accumarray(rl, wx, [numel(lo) 1]);
end
for i = 1:nexit, ex(i) = ie(st.exits(i)); end
end

function [ay, wy] = local_y(st, cfg, J6, dq, ax, p)
% y nodes and weights; weights include 1/T (slot transmission) so that the x
% weights (slot-filtered elastic counts) are turned into all incident electrons
h = cfg.h;
switch st.ytype
    case 'gauss'
        [x, w] = local_gh(st.n_gh); s = cfg.sigma_y; ay = sqrt(2) * s * x; T = erf(h / (sqrt(2) * s));
        wy = w / sqrt(pi) / T;
    case {'narrow_pearson', 'mixture'}
        ay = unique([-0.03:5e-4:-0.003, -0.003:5e-5:0.003, 0.003:5e-4:0.03]).'; dy = gradient(ay);
        gn = mean(exp(J6.shape(3:4))) * dq; mn = exp(J6.shape(5));
        pn = (1 + (ay / gn) .^ 2) .^ (-mn); pn = pn / sum(pn .* dy); Tn = sum(pn(abs(ay) <= h) .* dy(abs(ay) <= h));
        if strcmp(st.ytype, 'narrow_pearson')
            wy = pn .* dy / Tn;
        else
            gh = mean(exp(J6.shape(6:7))) * dq; mh = exp(J6.shape(8));
            ph = (1 + (ay / gh) .^ 2) .^ (-mh); ph = ph / sum(ph .* dy); Th = sum(ph(abs(ay) <= h) .* dy(abs(ay) <= h));
            c = qe_angular_profile_model('components', ax, dq, p); f = c ./ sum(c, 2);   % slot-filtered shares at each x node
            wy = (f(:, 1) / Tn) .* (pn .* dy).' + (f(:, 2) / Th) .* (ph .* dy).';
        end
end
end

