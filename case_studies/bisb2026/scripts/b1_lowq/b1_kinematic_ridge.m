function T = b1_kinematic_ridge(run_root, out_dir)
% Does the sqrt(q)-like ridge of the column-normalised raw q-E map (undergraduate thesis,
% defence slide 11) come from the kinematic factor? Measured: 590 native channels, raw
% counts, each column normalised to its 200-1700 meV area. Synthetic: one loss function
% for every q (no dispersion at all) - the two DL branches of the P2b h-scan fit with the
% real slot (h = 0.0066) taken at Q = 0.005 - times the single-channel 2D kinematic factor
% of each channel (Gaussian probe 3.4e-4, h_perp = 0.0066) and the thermal factor, then
% normalised the same way. Ridge = maximum of each column after a 100 meV moving mean,
% 300-1700 meV, the same rule for both. Writes kinematic_ridge_590.png/.pdf to out_dir.
root = b1_simple_root();
if nargin < 1 || isempty(run_root), run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
if nargin < 2, out_dir = fullfile(run_root, 'stage5', 'figures'); end
qmax = 0.015; win = [200 1700]; rwin = [300 1700]; Q0 = 0.005; kT = 25.852;
% measured map
cf = fullfile(run_root, 'stage1', 'channels_590_PL2_10w.mat'); S = load(cf, 'C'); C = S.C; mf = matfile(cf);
sel = find(abs(C.q) <= qmax + 1e-9 & abs(C.q) > 1e-9); q = C.q(sel); E = C.E(:);
Y = zeros(numel(E), numel(sel)); for k = 1:numel(sel), Y(:, k) = sum(double(mf.frames(:, :, sel(k))), 2); end
% fixed loss function: P2b h-scan main fit at h = 0.0066, evaluated at Q0
P = load(fullfile(run_root, 'p2b', 'fit_590_79ee195c', 'p2b_results.mat'), 'cfg', 'res'); F = P.res(1).fit; nd = P.cfg.nodes;
E0sq = exp(reshape(F.u(1:10), 5, 2).'); G = exp(reshape(F.u(11:16), 3, 2).');
lq = log(Q0); nk = numel(nd.A); B = zeros(1, nk);
for k = 1:nk, e = zeros(1, nk); e(k) = 1; B(k) = interp1(log(nd.A), e, lq); end
L = zeros(numel(E), 1); par = zeros(2, 3); Ep = max(abs(E), 1);
for j = 1:2
    e0 = sqrt(exp(interp1(log(nd.E0), log(E0sq(j, :)), lq))); g = exp(interp1(log(nd.G), log(G(j, :)), lq));
    A = (Q0 / 0.001) * (B * F.a((j - 1) * nk + (1:nk))); par(j, :) = [e0 g A];
    L = L + A * g ./ ((Ep .^ 2 - e0 ^ 2) .^ 2 + Ep .^ 2 * g ^ 2);          % L / E, as C in qe_lowq_joint_forward
end
th = Ep ./ (-expm1(-Ep / kT)); Ys = zeros(size(Y));
for k = 1:numel(sel)
    Kf = qe_kinematic_prefactor(q(k), dq = C.dq, form = '2d', sigma_probe = 3.4e-4, h_perp = 0.0066);
    Ys(:, k) = Kf(Ep) .* L .* th;
end
m = E >= win(1) & E <= win(2); nrm = @(Z) Z ./ trapz(E(m), Z(m, :), 1);
Yn = nrm(Y); Ysn = nrm(Ys);
ridge = @(Z) local_ridge(E, Z, rwin);
rm = ridge(Yn); rs = ridge(Ysn);
% the fixed L has the same peak at every q: its maximum (same smoothing rule) and the DL apex of each branch
Lw = L .* Ep; Lpk = local_ridge(E, Lw, rwin);
x = 2 * par(:, 1) .^ 2 - par(:, 2) .^ 2; apx = sqrt((x + sqrt(x .^ 2 + 12 * par(:, 1) .^ 4)) / 6);
T = table(q(:), rm(:), rs(:), 'VariableNames', {'q', 'ridge_measured_meV', 'ridge_flat_model_meV'});
writetable(T, fullfile(out_dir, 'kinematic_ridge_590.csv'));
% figure
blue = [42 120 214] / 255; orange = [235 104 52] / 255;
f = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', 'Position', [0 0 1650 520]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
ttl = {'(a) measured 590, column-normalised raw counts', '(b) non-dispersing L \times kinematic factor'};
Z = {Yn, Ysn}; R = {rm, rs}; cl = {blue, orange};
for p = 1:2
    ax = nexttile(tl); imagesc(ax, q * 1e3, E(m), Z{p}(m, :)); set(ax, 'YDir', 'normal'); colormap(ax, 'sky');
    clim(ax, prctile(Z{p}(m, :), [2 99.5], 'all')); hold(ax, 'on');
    plot(ax, q * 1e3, R{p}, 'o', 'MarkerFaceColor', cl{p}, 'MarkerEdgeColor', 'w', 'MarkerSize', 5, 'LineStyle', 'none');
    xlabel(ax, 'q (10^{-3} Å^{-1})'); ylabel(ax, 'energy loss (meV)'); title(ax, ttl{p}, 'FontWeight', 'normal');
    ax.FontSize = 11; ylim(ax, win);
end
ax = nexttile(tl); hold(ax, 'on');
aq = unique(round(abs(q), 6)); am = arrayfun(@(a) mean(rm(abs(abs(q) - a) < 1e-9)), aq); as = arrayfun(@(a) mean(rs(abs(abs(q) - a) < 1e-9)), aq);
plot(ax, aq * 1e3, am, 'o', 'Color', blue, 'MarkerFaceColor', blue, 'MarkerSize', 6, 'DisplayName', 'measured ridge');
plot(ax, aq * 1e3, as, 's', 'Color', orange, 'MarkerFaceColor', orange, 'MarkerSize', 6, 'DisplayName', 'ridge of the non-dispersing model');
yline(ax, Lpk, '-', 'Color', [0.32 0.32 0.31], 'LineWidth', 1.2, 'DisplayName', sprintf('peak of the model L itself (every q): %.2f eV', Lpk / 1e3));
yline(ax, apx(1), ':', 'Color', [0.32 0.32 0.31], 'LineWidth', 1.2, 'DisplayName', sprintf('branch apexes of the model L: %.2f, %.2f eV', apx(1) / 1e3, apx(2) / 1e3));
yline(ax, apx(2), ':', 'Color', [0.32 0.32 0.31], 'LineWidth', 1.2, 'HandleVisibility', 'off');
grid(ax, 'on'); ax.GridColor = [0.85 0.85 0.85]; ax.GridAlpha = 1; box(ax, 'off'); ax.FontSize = 11;
xlim(ax, [0 qmax * 1e3]); ylim(ax, win); xlabel(ax, '|q| (10^{-3} Å^{-1})'); ylabel(ax, 'ridge energy (meV)');
title(ax, '(c) ridge: measured vs non-dispersing model', 'FontWeight', 'normal'); legend(ax, 'Location', 'southeast', 'Box', 'off', 'FontSize', 10);
exportgraphics(f, fullfile(out_dir, 'kinematic_ridge_590.png'), 'Resolution', 200);
exportgraphics(f, fullfile(out_dir, 'kinematic_ridge_590.pdf'), 'ContentType', 'vector'); close(f);
fprintf('fixed L at Q = %.3f: low E0 %.0f G %.0f apex %.0f, high E0 %.0f G %.0f apex %.0f; peak of L %.0f meV\n', Q0, par(1, 1), par(1, 2), apx(1), par(2, 1), par(2, 2), apx(2), Lpk);
fprintf('kinematic ridge figure written to %s\n', out_dir); disp(T(abs(T.q) <= 0.004 | mod(round(abs(T.q) * 1e4), 20) == 0, :));
end

function r = local_ridge(E, Z, w)
s = movmean(Z, 25, 1); m = E >= w(1) & E <= w(2); s(~m, :) = -Inf; [~, i] = max(s, [], 1); r = E(i).';
end
