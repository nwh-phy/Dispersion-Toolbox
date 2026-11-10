function out = b1_simple_s0(run_root)
% S0 of the simplified pipeline check: how far down in q is the simple model safe?
% Truth: the P2b main model counts (block operators, broad component, shared L(Q,E),
% fitted tails) of the 20 native channels, summed over the six blocks (no noise).
% Each channel with |q| >= 0.0015 is fitted with the simple model:
%   counts(E) = [K_c(|E|) * sum_j a_j G_j / ((E^2 - E0_j^2)^2 + E^2 G_j^2) * g(E)] (x) R + b (E/250)^-4
% K_c = single-channel 2D kinematic factor with a Gaussian probe (qe_kinematic_prefactor),
% g = thermal factor, R = Gaussian of the angle-integrated zero-loss width; window, tail,
% variance and starts as P2b. Bias = simple fit - P2b curve at the nominal |q|, judged
% against the P2b statistical error; the expected per-channel error of the simple fit
% (Jacobian with Poisson variance at the truth) and chi2/point of the noiseless fit
% (shape mismatch) are reported too. Variants change one element of the simple model.
% Output: <run_root>/simple/s0_<config hash>/.
root = b1_simple_root();
if nargin < 1, run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
cfg = b1_simple_config(); tag = b1_sha256(cfg); t0 = tic;
out = fullfile(run_root, 'simple', sprintf('s0_%s', tag(1:8))); if ~isfolder(out), mkdir(out); end
S = load(fullfile(run_root, cfg.p2b_dir, 'p2b_results.mat'), 'cfg', 'res', 'mainEo', 'mainMu');
Eo = S.mainEo(:); qx = S.cfg.exit_q(:).'; C = S.res(1).curves;
sel = find(abs(qx) >= cfg.s0_qmin_test - 1e-9); q = qx(sel);
Y = reshape(sum(S.mainMu(:, sel, :), 3), numel(Eo), []);
S2 = max(movmean(Y, cfg.noise_smooth, 1), 1);
mw = Eo >= cfg.window(1) & Eo <= cfg.window(2); Rk = local_kernel(cfg, Eo(2) - Eo(1));
opts = optimoptions('lsqnonlin', 'Display', 'off', 'FunctionTolerance', 1e-12, 'StepTolerance', 1e-10, ...
    'OptimalityTolerance', 1e-10, 'MaxIterations', 400, 'MaxFunctionEvaluations', 4000);
names = {'E0_b1', 'G_b1', 'apex_b1', 'E0_b2', 'G_b2', 'apex_b2'};
rows = {}; nv = numel(cfg.s0_variants); fits = cell(nv, numel(sel));
for v = 1:nv
    V = cfg.s0_variants(v);
    for k = 1:numel(sel)
        F = local_fit(Eo, Y(:, k), S2(:, k), mw, local_K(q(k), V, cfg), Rk, q(k), cfg, opts); fits{v, k} = F;
        T = local_truth(C, abs(q(k)));
        for n = 1:numel(names)
            nm = names{n}; jd = any(strcmp(nm, cfg.judged)); b = F.(nm) - T.(nm);
            rows(end + 1, :) = {V.name, q(k), nm, F.(nm), F.(['s' nm]), T.(nm), b, T.(['s' nm]), b / T.(['s' nm]), ...
                jd, abs(b) < T.(['s' nm]), F.chi2n, F.start}; %#ok<AGROW>
        end
    end
    fprintf('%-40s done (%.0f s)\n', V.name, toc(t0));
end
R = cell2table(rows, 'VariableNames', {'variant', 'q', 'quantity', 'simple', 'simple_err_expected', 'p2b', 'bias', ...
    'p2b_stat_err', 'bias_over_p2b_err', 'judged', 'within_p2b_err', 'chi2_per_point_noiseless', 'best_start'});
writetable(R, fullfile(out, 's0_bias.csv'));

% q_min per variant: smallest |q| from which every judged quantity on both signs is within the P2b error
lev = unique(round(abs(q(:)), 6)); srow = {};
for v = 1:nv
    nm = cfg.s0_variants(v).name; ok = false(size(lev));
    for l = 1:numel(lev)
        x = R(strcmp(R.variant, nm) & abs(abs(R.q) - lev(l)) < 1e-9 & R.judged, :); ok(l) = all(x.within_p2b_err);
    end
    i0 = find(arrayfun(@(l) all(ok(l:end)), 1:numel(lev)), 1); qmin = NaN; if ~isempty(i0), qmin = lev(i0); end
    for l = 1:numel(lev)
        x = R(strcmp(R.variant, nm) & abs(abs(R.q) - lev(l)) < 1e-9, :);
        g = @(f) max(abs(x.bias(strcmp(x.quantity, f))));
        srow(end + 1, :) = {nm, lev(l), ok(l), qmin, g('E0_b1'), g('apex_b1'), g('apex_b2'), g('E0_b2'), max(x.chi2_per_point_noiseless)}; %#ok<AGROW>
    end
end
Ssum = cell2table(srow, 'VariableNames', {'variant', 'abs_q', 'judged_all_within', 'q_min', 'max_abs_bias_E0_b1', ...
    'max_abs_bias_apex_b1', 'max_abs_bias_apex_b2', 'max_abs_bias_E0_b2', 'max_chi2_per_point_noiseless'});
writetable(Ssum, fullfile(out, 's0_summary.csv'));

local_fig_bias(out, R, cfg, lev);
local_fig_spectra(out, Eo, Y, S2, mw, fits(1, :), q, cfg);
save(fullfile(out, 's0_results.mat'), 'cfg', 'R', 'Ssum', 'fits', 'q', 'Eo', 'Y', 'S2');
meta = struct('config', cfg, 'config_sha256', tag, 'code_sha256', b1_sha256({'file', [mfilename('fullpath') '.m']}), ...
    'truth', fullfile(run_root, cfg.p2b_dir, 'p2b_results.mat'), 'runtime_s', toc(t0), ...
    'created', char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
fid = fopen(fullfile(out, 'config.json'), 'w'); fprintf(fid, '%s', jsonencode(meta, 'PrettyPrint', true)); fclose(fid);
fprintf('S0 written to %s (%.0f s)\n', out, toc(t0));
end


% =============================== model ===============================
function Rk = local_kernel(cfg, step)
e = -cfg.res_half:step:cfg.res_half; s = cfg.res_fwhm / (2 * sqrt(2 * log(2)));
Rk = exp(-e(:) .^ 2 / (2 * s ^ 2)); Rk = Rk / sum(Rk);
end

function Kf = local_K(qc, V, cfg)
[Kf, info] = qe_kinematic_prefactor(qc, dq = cfg.dq, form = '2d', sigma_probe = V.sigma, ...
    h_perp = V.h_perp, sigma_probe_y = V.sigma_y);
if V.qweight
    % residue proportional to Q inside the bin (2D optical limit, as the P2b intrinsic model)
    qn = info.q_nodes(:).'; wn = info.q_weights(:); qE = info.qE_per_meV; qm = qn * wn;
    Kf = @(E) reshape(((qn .^ 2 ./ (qn .^ 2 + (qE * abs(E(:))) .^ 2) .^ 2) * wn) / qm, size(E));
end
end

function F = local_fit(Eo, y, s2, mw, Kf, Rk, qc, cfg, opts)
step = Eo(2) - Eo(1); aw = max(abs(Eo), step / 4); Kv = Kf(aw);
g = Eo ./ (-expm1(-Eo / cfg.kT)); g(Eo == 0) = cfg.kT;
T = (Eo(mw) / cfg.tail.E_ref) .^ (-cfg.tail.p); w = 1 ./ sqrt(s2(mw));
D = struct('aw', aw, 'Kg', Kv .* g, 'Rk', Rk, 'mw', mw, 'yw', y(mw) .* w, 'w', w, 'tw', T .* w);
lb = log([cfg.bounds_E0(:, 1).' cfg.bounds_G(1) cfg.bounds_G(1)]); ub = log([cfg.bounds_E0(:, 2).' cfg.bounds_G(2) cfg.bounds_G(2)]);
lq = log(abs(qc)); best = [];
for s = 1:numel(cfg.starts)
    st = cfg.starts(s); e0 = zeros(1, 2); g0 = e0;
    for j = 1:2
        e0(j) = interp1(log(cfg.nodes.E0), st.E0(j, :), min(max(lq, log(cfg.nodes.E0(1))), log(cfg.nodes.E0(end))));
        g0(j) = interp1(log(cfg.nodes.G), st.G(j, :), min(max(lq, log(cfg.nodes.G(1))), log(cfg.nodes.G(end))));
    end
    u0 = min(max(log([e0 g0]), lb + 1e-6), ub - 1e-6);
    [u, rn, ~, flag, ~, ~, J] = lsqnonlin(@(u) local_res(u, D), u0, lb, ub, opts);
    if isempty(best) || rn < best.rn, best = struct('u', u, 'rn', rn, 'J', J, 'flag', flag, 'start', st.name); end
end
u = best.u; [~, a, b, mu] = local_res(u, D);
Cov = inv(full(best.J.' * best.J));                 % expected covariance: Poisson variance at the truth
v0 = local_derived(u); gr = zeros(numel(v0), 4); h = 1e-6;
for k = 1:4, up = u; up(k) = up(k) + h; gr(:, k) = (local_derived(up) - v0) / h; end
sv = sqrt(max(sum((gr * Cov) .* gr, 2), 0));
nm = {'E0_b1', 'G_b1', 'apex_b1', 'E0_b2', 'G_b2', 'apex_b2'}; F = struct();
for k = 1:6, F.(nm{k}) = v0(k); F.(['s' nm{k}]) = sv(k); end
F.a = a; F.b = b; F.mu = mu; F.chi2n = best.rn / nnz(mw); F.flag = best.flag; F.start = best.start; F.u = u;
end

function [r, a, b, mu] = local_res(u, D)
E0 = exp(u(1:2)); G = exp(u(3:4)); A = zeros(nnz(D.mw), 2); full_cols = zeros(numel(D.aw), 2);
for j = 1:2
    core = D.Kg .* G(j) ./ ((D.aw .^ 2 - E0(j) ^ 2) .^ 2 + D.aw .^ 2 * G(j) ^ 2);
    c = conv(core, D.Rk, 'same'); full_cols(:, j) = c; A(:, j) = c(D.mw);
end
Aw = A .* D.w; tt = D.tw.' * D.tw;
ywp = D.yw - D.tw * (D.tw.' * D.yw) / tt; Awp = Aw - D.tw * (D.tw.' * Aw) / tt;
sc = max(abs(Awp), [], 1); sc(sc == 0) = 1; a = lsqnonneg(Awp ./ sc, ywp) ./ sc(:);
r = Awp * a - ywp;
if nargout > 2
    b = (D.tw.' * (D.yw - Aw * a)) / tt; mu = full_cols * a; mu(D.mw) = mu(D.mw) + b * D.tw ./ D.w;
end
end

function v = local_derived(u)
E0 = exp(u(1:2)); G = exp(u(3:4)); x = 2 * E0 .^ 2 - G .^ 2; apex = sqrt((x + sqrt(x .^ 2 + 12 * E0 .^ 4)) / 6);
v = [E0(1); G(1); apex(1); E0(2); G(2); apex(2)];
end

function T = local_truth(C, aq)
lQ = log(C.Q); x = log(aq); T = struct();
for j = 1:2
    for f = {'E0', 'G', 'apex'}
        T.(sprintf('%s_b%d', f{1}, j)) = interp1(lQ, C.(f{1})(j, :), x);
        T.(sprintf('s%s_b%d', f{1}, j)) = interp1(lQ, C.(['s' f{1}])(j, :), x);
    end
end
end


% =============================== figures ===============================
function local_fig_bias(out, R, cfg, lev)
qs = {'E0_b1', 'apex_b1', 'apex_b2', 'E0_b2'};
ttl = {'low branch E0 (judged)', 'low branch apex (judged)', 'high branch apex (judged)', 'high branch E0 (reported only)'};
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1600 1000]); tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
cl = lines(numel(cfg.s0_variants));
for t = 1:4
    ax = nexttile(tl); hold(ax, 'on');
    x = R(strcmp(R.variant, cfg.s0_variants(1).name) & strcmp(R.quantity, qs{t}), :);
    band = arrayfun(@(l) mean(x.p2b_stat_err(abs(abs(x.q) - l) < 1e-9)), lev);
    fill(ax, [lev; flipud(lev)] * 1e3, [band; -flipud(band)], [0.6 0.6 0.6], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'DisplayName', '+/- P2b stat. error');
    for v = 1:numel(cfg.s0_variants)
        x = R(strcmp(R.variant, cfg.s0_variants(v).name) & strcmp(R.quantity, qs{t}), :);
        p = x.q > 0; m = x.q < 0; o = (v - 2.5) * 0.03;
        if v == 1
            errorbar(ax, abs(x.q(p)) * 1e3 + o, x.bias(p), x.simple_err_expected(p), 'o', 'Color', cl(v, :), 'LineWidth', 1.3, 'DisplayName', [cfg.s0_variants(v).name ', +q (bar: expected error)']);
            errorbar(ax, abs(x.q(m)) * 1e3 + o, x.bias(m), x.simple_err_expected(m), 's', 'Color', cl(v, :), 'LineWidth', 1.3, 'DisplayName', [cfg.s0_variants(v).name ', -q']);
        else
            mb = arrayfun(@(l) mean(x.bias(abs(abs(x.q) - l) < 1e-9)), lev);
            plot(ax, lev * 1e3 + o, mb, '.-', 'Color', cl(v, :), 'DisplayName', [cfg.s0_variants(v).name ' (mean of +/-q)']);
        end
    end
    yline(ax, 0, 'k-', 'HandleVisibility', 'off'); grid(ax, 'on'); xlim(ax, [1.2 5.3]);
    xlabel(ax, '|q| (10^{-3} A^{-1})'); ylabel(ax, 'simple - P2b (meV)'); title(ax, ttl{t});
    if t == 1, legend(ax, 'FontSize', 7, 'Location', 'best'); end
end
exportgraphics(f, fullfile(out, 'fig_s0_bias.png'), 'Resolution', 110); close(f);
end

function local_fig_spectra(out, Eo, Y, S2, mw, fits, q, cfg)
show = find(ismember(round(abs(q), 6), round(cfg.fig_q, 6)));
f = figure('Visible', 'off', 'Theme', 'light', 'Position', [0 0 1800 800]); tl = tiledlayout(f, 2, numel(show) / 2, 'TileSpacing', 'compact');
for k = show(:).'
    ax = nexttile(tl); F = fits{k}; m = mw;
    yyaxis(ax, 'left'); plot(ax, Eo(m), Y(m, k), 'k-', Eo(m), F.mu(m), 'r--', 'LineWidth', 1.2); ylabel(ax, 'counts (black: truth, red: simple fit)');
    ax.YAxis(1).Color = 'k';
    yyaxis(ax, 'right'); plot(ax, Eo(m), (Y(m, k) - F.mu(m)) ./ sqrt(S2(m, k)), '-', 'Color', [0 0.45 0.74]); ylabel(ax, '(truth - simple) / sqrt(counts)');
    ax.YAxis(2).Color = [0 0.45 0.74];
    grid(ax, 'on'); xlim(ax, cfg.window);
    title(ax, sprintf('q = %+.4f: chi2/pt %.2f, E0 %.0f / %.0f meV', q(k), F.chi2n, F.E0_b1, F.E0_b2), 'FontSize', 9);
end
exportgraphics(f, fullfile(out, 'fig_s0_spectra.png'), 'Resolution', 100); close(f);
end
