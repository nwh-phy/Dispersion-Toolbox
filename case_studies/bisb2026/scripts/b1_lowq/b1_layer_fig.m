function b1_layer_fig(run_root, out_dir)
% Figure of the 1L / 2L comparison (reports/layer_1L_2L_20261017.md): DL apex of the low
% and high B1 branch against |q| for 590 (1L), n0 (1L, 0.85 um away) and 20w (2L), from the
% 2d_h5 adaptive-binning tables in <run_root>/stage5. Each point is the mean of the +q and
% -q bins at that |q|; the bar is half their difference (the +/-q asymmetry, as in the
% error budget of He et al., SI Fig. S7). Colour encodes the layer number, the marker fill
% the 1L region. Writes layer_1L_2L.png / .pdf to out_dir (default <run_root>/stage5/figures).
root = b1_simple_root();
if nargin < 1 || isempty(run_root), run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
if nargin < 2, out_dir = fullfile(run_root, 'stage5', 'figures'); end
if ~isfolder(out_dir), mkdir(out_dir); end
blue = [42 120 214] / 255; orange = [235 104 52] / 255;    % categorical slots 1-2, validated (light)
S = struct('file', {'590_PL2_10w', 'n0_PL2_10w_repeat', 'no_PL2_20w_2film'}, ...
    'label', {'1L (590)', '1L (n0, 0.85 \mum away)', '2L (20w)'}, 'color', {blue, blue, orange}, ...
    'marker', {'o', 'o', 's'}, 'fill', {blue, 'w', orange});
qr = [0.002 0.0125]; key = {'apex_low', 'apex_high'}; ttl = {'(a) low branch', '(b) high branch'};
f = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', 'Position', [0 0 1200 470]);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for p = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for s = 1:numel(S)
        T = readtable(fullfile(run_root, 'stage5', sprintf('e0_stage2_2d_h5_%s.csv', S(s).file)));
        T = T(strcmp(T.status, 'converged') & abs(T.q) >= qr(1) - 1e-9 & abs(T.q) <= qr(2) + 1e-9, :);
        aq = unique(round(abs(T.q), 5)); y = nan(size(aq)); e = y;
        for k = 1:numel(aq)
            v = T.(key{p})(abs(round(abs(T.q), 5) - aq(k)) < 1e-9);
            y(k) = mean(v); if numel(v) == 2, e(k) = abs(diff(v)) / 2; end
        end
        errorbar(ax, aq * 1e3, y / 1e3, e / 1e3, S(s).marker, 'Color', S(s).color, 'MarkerFaceColor', S(s).fill, ...
            'MarkerSize', 7, 'LineWidth', 1.2, 'CapSize', 0, 'DisplayName', S(s).label);
        plot(ax, aq * 1e3, y / 1e3, '-', 'Color', [S(s).color 0.35], 'LineWidth', 1, 'HandleVisibility', 'off');
    end
    ax.GridColor = [0.85 0.85 0.85]; ax.GridAlpha = 1; grid(ax, 'on'); box(ax, 'off');
    ax.FontSize = 12; ax.XColor = [0.32 0.32 0.31]; ax.YColor = [0.32 0.32 0.31];
    xlim(ax, [1.5 13]); xlabel(ax, '|q| (10^{-3} Å^{-1})', 'Color', 'k'); ylabel(ax, 'DL apex (eV)', 'Color', 'k');
    title(ax, ttl{p}, 'FontWeight', 'normal', 'Color', 'k');
    if p == 1, legend(ax, 'Location', 'northwest', 'Box', 'off', 'FontSize', 11); end
end
exportgraphics(f, fullfile(out_dir, 'layer_1L_2L.png'), 'Resolution', 200);
exportgraphics(f, fullfile(out_dir, 'layer_1L_2L.pdf'), 'ContentType', 'vector');
close(f); fprintf('layer figure written to %s\n', out_dir);
end
