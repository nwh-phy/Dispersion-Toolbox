function T = b1_p2b_lossmap_fig(fit_dir, q_single, fig_cfg)
%B1_P2B_LOSSMAP_FIG  Before and after the kinematic factor, with the P2b forward model.
%   T = b1_p2b_lossmap_fig(fit_dir) uses the main variant stored in
%   fit_dir/p2b_results.mat, with all stored time blocks summed.
%   fig_p2b_lossmap: (a) measured B1 counts, (b) forward-model counts, (c) intrinsic
%     loss function L(Q,E) of the same fit at each channel's |q| (kinematic factor,
%     thermal factor and measurement operator removed, b1_p2b_intrinsic_L), (d) shape
%     residual (data - s*model)/sigma in 40 meV windows, s = one energy-independent
%     scale per channel (data/model in the window, column scale_data_over_model of T;
%     the fit shares the residues between +q and -q, and in 590 the -q side is 2-5 %
%     above the model and the +q side 1-4 % below). B1 counts = counts minus the fitted per-spectrum
%     power-law tail. (a)-(c): 100 meV moving mean, each column normalised to its
%     maximum in 250-1800 meV; dots in (a), (b) = column maxima (what peak tracking
%     sees); lines in (c) = DL apex (solid) and E0 (dashed) of the two branches; the
%     shaded columns |q| <= 0.0015 have the high branch at the window edge.
%   fig_p2b_single_spectrum: one channel before (counts, model) and after (L at the
%     channel's |q|, its two branches, the point-Q kinematic factor). q_single defaults
%     to the channel whose counts maximum lies furthest from the maximum of L, among
%     |q| >= 0.002 with the high-branch apex inside the window: below 0.002 the low
%     branch depends on the slot acceptance h and the high branch runs to the window
%     edge (p2b_hscan_20261017.md), so L there is not a reliable "after".
%   T: per channel, the maxima of counts, model and L and the branch apexes and E0
%     (meV), also written to p2b_lossmap_peaks.csv. Figures go to fit_dir.
arguments
    fit_dir char = ''
    q_single double = []
    fig_cfg struct = struct()
end
paths=b1_forward_paths(fig_cfg);
if isempty(fit_dir), fit_dir=paths.fit_dir; end
out=fit_dir;
if isfield(fig_cfg,'output_dir'), out=fig_cfg.output_dir; end
S = load(fullfile(fit_dir, 'p2b_results.mat'), 'cfg', 'res', 'mainEo', 'mainY', 'mainS2', 'mainMu');
cfg = S.cfg; F = S.res(1).fit; Cv = S.res(1).curves; E = S.mainEo(:); q = cfg.exit_q(:).';
[~, ni, nb] = size(S.mainY);
tl = (E / cfg.tail.E_ref) .^ (-cfg.tail.p); tl(E < cfg.holdout(1)) = 0;
tail = F.tailOn * tl .* reshape(F.b, 1, ni, nb);
Y = sum(S.mainY - tail, 3); M = sum(S.mainMu - tail, 3); V = sum(S.mainS2, 3);   % nE x ni
window=cfg.variants(1).window; w = E >= window(1) & E <= window(2); Ew = E(w);
fit_paths=b1_forward_paths(cfg); d=bisb_datasets(fit_paths.session);
if isempty(fieldnames(fig_cfg)) && ~isfield(cfg,'paths'), fig_cfg.edge_q=0.0015; end
nsm = round(100 / median(diff(E)));
Ys = movmean(Y, nsm, 1); Ms = movmean(M, nsm, 1);
[L, Lj] = b1_p2b_intrinsic_L(fit_dir, abs(q(:)), Ew.');
L = L.'; Lj = permute(Lj, [2 1 3]);                                   % nEw x ni (x 2 branches)
[~, iy] = max(Ys(w, :), [], 1); [~, im] = max(Ms(w, :), [], 1); [~, il] = max(L, [], 1);
cl = @(f, j) interp1(Cv.Q, Cv.(f)(j, :), min(max(abs(q), Cv.Q(1)), Cv.Q(end))).';
T = table(q(:), Ew(iy), Ew(im), Ew(il), cl('apex', 1), cl('apex', 2), cl('E0', 1), cl('E0', 2), ...
    'VariableNames', {'q', 'counts_max_meV', 'model_max_meV', 'L_max_meV', 'apex_low_meV', 'apex_high_meV', ...
    'E0_low_meV', 'E0_high_meV'});
T.shift_meV = T.L_max_meV - T.counts_max_meV;
T.scale_data_over_model = (sum(Y(w, :), 1) ./ sum(M(w, :), 1)).';
if ~isfolder(out), mkdir(out); end
writetable(T, fullfile(out, 'p2b_lossmap_peaks.csv'));
if isempty(q_single)
    c = find(abs(q(:)) >= 0.002 - 1e-9 & T.apex_high_meV < Ew(end) - 4);
    if isempty(c), [~, imax]=max(abs(q)); c=imax; end
    [~, k] = max(abs(T.shift_meV(c)) + 1e-6 * (q(c).' > 0)); q_single = q(c(k));
end
local_maps(out, q, d.dq, Ew, Ys(w, :), Ms(w, :), L, Y(w, :) - T.scale_data_over_model.' .* M(w, :), V(w, :), T, Cv, fig_cfg);
local_single(out, q, Ew, Y(w, :), Ys(w, :), M(w, :), L, Lj, T, q_single, nb);
disp(T);
end


function local_maps(fit_dir, q, dq, Ew, Ys, Ms, L, R, V, T, Cv, cfg)
grid = round(min(q)/dq):round(max(q)/dq); qg=grid*dq; [~, col] = ismember(round(q/dq),grid);
put = @(Z) local_grid(Z, col, numel(qg));
nrm = @(Z) Z ./ max(Z, [], 1);
Rw = movsum(R, 10, 1) ./ sqrt(movsum(V, 10, 1));                       % 40 meV windows
maps = {put(nrm(Ys)), put(nrm(Ms)), put(nrm(L)), put(Rw)};
mx = {T.counts_max_meV, T.model_max_meV, [], []};
ttl = {'(a) measured counts, B1', '(b) P2b forward model, counts', ...
    '(c) intrinsic L(Q,E): apex (solid), E_0 (dashed)', '(d) shape residual / \sigma, 40 meV'};
f = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 42 11.5], 'Color', 'w');
if isprop(f, 'Theme'), f.Theme = 'light'; end
tlo = tiledlayout(f, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:4
    ax = nexttile(tlo); Z = maps{p};
    imagesc(ax, qg * 1e3, Ew / 1e3, Z, 'AlphaData', ~isnan(Z)); axis(ax, 'xy'); hold(ax, 'on');
    if p < 4
        colormap(ax, turbo); clim(ax, [0.3 1]);
        if ~isempty(mx{p}), plot(ax, q * 1e3, mx{p} / 1e3, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 4); end
    else
        colormap(ax, local_diverging()); clim(ax, [-4 4]);
    end
    if p == 3
        keep = Cv.Q <= max(abs(q))+dq/2;
        for j = 1:2
            for s = [-1 1]
                plot(ax, s * Cv.Q(keep) * 1e3, Cv.apex(j, keep) / 1e3, '-', 'Color', 'w', 'LineWidth', 1.5);
                plot(ax, s * Cv.Q(keep) * 1e3, Cv.E0(j, keep) / 1e3, '--', 'Color', 'w', 'LineWidth', 1.0);
            end
        end
        if isfield(cfg,'edge_q') && cfg.edge_q>0
            edge=(cfg.edge_q+dq/2)*1e3; el=[Ew(1) Ew(end)]/1e3;
            patch(ax,[-edge edge edge -edge],el([1 1 2 2]), 'w','FaceAlpha',0.6,'EdgeColor','none');
            text(ax,0,mean(el),'high branch at window edge','Rotation',90,'HorizontalAlignment','center', ...
                'VerticalAlignment','middle','FontSize',8,'Color',[0.2 0.2 0.2]);
        end
    end
    set(ax, 'Color', 'w', 'FontName', 'Helvetica', 'FontSize', 9, 'TickDir', 'out', 'Layer', 'top'); box(ax, 'on');
    xlim(ax, ([qg(1) qg(end)]+[-1 1]*dq/2)*1e3); ylim(ax,[Ew(1) Ew(end)]/1e3); xlabel(ax, 'q (10^{-3} Å^{-1})');
    if p == 1, ylabel(ax, 'energy loss (eV)'); end
    title(ax, ttl{p}, 'FontWeight', 'normal', 'FontSize', 10);
    cb = colorbar(ax); cb.FontSize = 8;
    if p == 3, cb.Label.String = 'normalised in each column'; elseif p == 4, cb.Label.String = '\sigma'; end
end
exportgraphics(f, fullfile(fit_dir, 'fig_p2b_lossmap.png'), 'Resolution', 200, 'BackgroundColor', 'white');
exportgraphics(f, fullfile(fit_dir, 'fig_p2b_lossmap.pdf'), 'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
close(f);
end


function local_single(fit_dir, q, Ew, Y, Ys, M, L, Lj, T, q1, nb)
i = find(abs(q - q1) < 1e-9, 1); e = Ew / 1e3; Q = abs(q1);
intr = qe_lowq_joint_forward('intrinsic');
K = Q ./ (Q ^ 2 + (Ew * intr.qE_per_meV) .^ 2) .^ 2;                  % point-Q 2D kinematic factor
f = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 30 13], 'Color', 'w');
if isprop(f, 'Theme'), f.Theme = 'light'; end
tlo = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tlo); hold(ax, 'on');
plot(ax, e, Y(:, i), '.', 'Color', [0.72 0.72 0.72], 'MarkerSize', 5, 'DisplayName', 'counts, 4 meV');
plot(ax, e, Ys(:, i), 'k-', 'LineWidth', 1.6, 'DisplayName', 'counts, 100 meV mean');
plot(ax, e, M(:, i), 'r-', 'LineWidth', 1.2, 'DisplayName', 'P2b forward model');
xline(ax, T.counts_max_meV(i) / 1e3, 'k:', sprintf('maximum %.2f eV', T.counts_max_meV(i) / 1e3), 'LineWidth', 1.2, ...
    'LabelVerticalAlignment', 'top', 'FontSize', 9, 'HandleVisibility', 'off');
local_axes(ax, sprintf('(a) before: measured spectrum, q = %+.4f Å^{-1}', q1), sprintf('counts (B1, %d blocks, per %g meV)',nb,median(diff(Ew))), 'south');
ax = nexttile(tlo); hold(ax, 'on'); s = max(L(:, i)); cb = [0.10 0.35 0.85]; co = [0.85 0.35 0.10];
plot(ax, e, Ys(:, i) / max(Ys(:, i)), '-', 'Color', [0.65 0.65 0.65], 'LineWidth', 1.4, 'DisplayName', 'measured (a), normalised');
plot(ax, e, K / K(1), ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.3, 'DisplayName', 'kinematic factor K(Q,E) / K(Q, 0.25 eV)');
plot(ax, e, Lj(:, i, 1) / s, '--', 'Color', cb, 'LineWidth', 1.3, 'DisplayName', 'L, low branch');
plot(ax, e, Lj(:, i, 2) / s, '--', 'Color', co, 'LineWidth', 1.3, 'DisplayName', 'L, high branch');
plot(ax, e, L(:, i) / s, 'k-', 'LineWidth', 1.8, 'DisplayName', 'L(Q,E), total');
xline(ax, T.counts_max_meV(i) / 1e3, ':', sprintf('before: %.2f eV', T.counts_max_meV(i) / 1e3), 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom', 'FontSize', 9, 'HandleVisibility', 'off');
xline(ax, T.apex_low_meV(i) / 1e3, ':', sprintf('low apex %.2f eV', T.apex_low_meV(i) / 1e3), 'Color', cb, ...
    'LineWidth', 1.2, 'LabelVerticalAlignment', 'top', 'FontSize', 9, 'HandleVisibility', 'off');
xline(ax, T.apex_high_meV(i) / 1e3, ':', sprintf('high apex %.2f eV', T.apex_high_meV(i) / 1e3), 'Color', co, ...
    'LineWidth', 1.2, 'LabelVerticalAlignment', 'top', 'FontSize', 9, 'HandleVisibility', 'off');
local_axes(ax, sprintf('(b) after: intrinsic loss function at Q = %.4f Å^{-1}', Q), 'normalised', 'southoutside');
ylim(ax, [0 1.05]);
tag = strrep(sprintf('%+.4f', q1), '.', 'p');
exportgraphics(f, fullfile(fit_dir, ['fig_p2b_single_spectrum_q' tag '.png']), 'Resolution', 200, 'BackgroundColor', 'white');
exportgraphics(f, fullfile(fit_dir, ['fig_p2b_single_spectrum_q' tag '.pdf']), 'ContentType', 'vector', 'BackgroundColor', 'white');
close(f);
end


function local_axes(ax, ttl, ylab, loc)
set(ax, 'FontName', 'Helvetica', 'FontSize', 9, 'TickDir', 'out', 'Box', 'on'); grid(ax, 'on');
xlim(ax, [0.25 1.8]); xlabel(ax, 'energy loss (eV)'); ylabel(ax, ylab);
title(ax, ttl, 'FontWeight', 'normal', 'FontSize', 10);
lg = legend(ax, 'Location', loc, 'FontSize', 8, 'Box', 'off');
if strcmp(loc, 'southoutside'), lg.NumColumns = 2; end
end


function G = local_grid(Z, col, n)
G = nan(size(Z, 1), n); G(:, col) = Z;
end


function c = local_diverging()
t = linspace(0, 1, 129).';
c = [interp1([0 0.5 1], [0.13 1 0.70], t), interp1([0 0.5 1], [0.40 1 0.15], t), interp1([0 0.5 1], [0.75 1 0.12], t)];
end
