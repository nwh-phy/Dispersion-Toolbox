function info = qe_plot_heatmap(E, q, M, out_base, opts)
%QE_PLOT_HEATMAP  Column-normalised q-E heatmap, exported as PNG (300 dpi) and PDF.
%   info = qe_plot_heatmap(E, q, M, out_base) crops M (E x q counts) to
%   q_window x e_window, divides every column by its trapz area in norm_window
%   (default: the displayed e_window), and draws it with turbo and colour limits
%   at the percentiles pct of the displayed values (or fixed clim), no overlays.
%   log = true shows log10 of the normalised map (use a norm_window outside the
%   ZLP tail, e.g. [300 1700], when e_window starts below ~200 meV).
%   This is the display rule of the B1 heatmaps (b1_qe_heatmaps_three_dataset_260521,
%   2026-11-07 Bi figures). info: q, E, map, clim.
%
%   See also: qe_loss_map, build_comparison_qe
arguments
    E (:,1) double
    q (1,:) double
    M {mustBeNumeric}
    out_base (1,:) char
    opts.q_window (1,2) double = [-0.015 0.015]
    opts.e_window (1,2) double = [200 1700]
    opts.norm_window (1,2) double = [NaN NaN]
    opts.pct (1,2) double = [2 98]
    opts.clim (1,2) double = [NaN NaN]
    opts.log (1,1) logical = false
    opts.title (1,:) char = ''
    opts.visible (1,1) logical = false
end
assert(isequal(size(M), [numel(E) numel(q)]), 'qe_plot_heatmap:Size', 'M must be numel(E) x numel(q)');
nw = opts.norm_window; if any(isnan(nw)), nw = opts.e_window; end
dq = median(diff(q));
qm = q >= opts.q_window(1) - 1e-12 & q <= opts.q_window(2) + 1e-12;
em = E >= opts.e_window(1) & E <= opts.e_window(2); nm = E >= nw(1) & E <= nw(2);
H = double(M(:, qm));
a = trapz(E(nm), H(nm, :), 1); a(~(isfinite(a) & a > 0)) = 1;
H = H(em, :) ./ a; if opts.log, H = log10(max(H, eps)); end
qs = q(qm); es = E(em);
cl = opts.clim;
if any(isnan(cl))
    v = sort(H(isfinite(H)));
    cl = [v(max(1, round(opts.pct(1)/100*numel(v)))) v(max(1, round(opts.pct(2)/100*numel(v))))];
end
vis = 'off'; if opts.visible, vis = 'on'; end
fig = figure('Visible', vis, 'Color', 'w', 'Units', 'pixels', 'Position', [100 100 900 760]);
if isprop(fig, 'Theme'), fig.Theme = 'light'; end   % white export also in macOS dark mode
ax = axes(fig); imagesc(ax, qs, es, H);
set(ax, 'YDir', 'normal', 'FontSize', 16, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'on', 'Layer', 'top');
colormap(ax, turbo); clim(ax, cl); cb = colorbar(ax); cb.FontSize = 14;
if opts.log
    cb.Label.String = sprintf('log_{10}(intensity / area in %g–%g meV)', nw);
else
    cb.Label.String = sprintf('Intensity / area in %g–%g meV (meV^{-1})', nw);
end
xlabel(ax, 'q (Å^{-1})'); ylabel(ax, 'Energy loss (meV)');
xlim(ax, [min(qs) max(qs)] + [-0.5 0.5]*dq); ylim(ax, [min(es) max(es)]);
if ~isempty(opts.title), title(ax, opts.title, 'FontWeight', 'normal', 'Interpreter', 'none'); end
folder = fileparts(out_base); if ~isempty(folder) && ~isfolder(folder), mkdir(folder); end
exportgraphics(fig, [out_base '.png'], 'Resolution', 300);
exportgraphics(fig, [out_base '.pdf'], 'ContentType', 'vector');
close(fig);
info = struct('q', qs, 'E', es, 'map', H, 'clim', cl, 'norm_window', nw);
end
