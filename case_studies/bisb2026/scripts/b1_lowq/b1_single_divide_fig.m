function out = b1_single_divide_fig(fit_dir, q1, h_perp, sigma_probe)
%B1_SINGLE_DIVIDE_FIG  One channel divided by the kinematic factor, no fitting.
%   b1_single_divide_fig(fit_dir, q1) takes the measured counts of channel q1 from
%   fit_dir/p2b_results.mat (590, the 6 time blocks summed; nothing subtracted) and
%   divides them by the 2D kinematic factor of that channel, averaged over the channel
%   width, the slot acceptance +/- h_perp and a Gaussian probe sigma_probe
%   (qe_kinematic_prefactor, defaults h_perp = 0.0066, sigma_probe = 3.4e-4 A^-1 as in
%   b1_kinematic_ridge). No model and no background enter. Both curves get a 100 meV
%   moving mean and are normalised to their maximum in 250-1800 meV.
%   Writes fig_single_divide_q<q>.png/.pdf to fit_dir.
arguments
    fit_dir char
    q1 (1, 1) double = -0.004
    h_perp (1, 1) double = 0.0066
    sigma_probe (1, 1) double = 3.4e-4
end
S = load(fullfile(fit_dir, 'p2b_results.mat'), 'cfg', 'mainEo', 'mainY');
E = S.mainEo(:); i = find(abs(S.cfg.exit_q - q1) < 1e-9, 1);
w = E >= 250 & E <= 1800; Ew = E(w); e = Ew / 1e3;
y = sum(S.mainY(w, i, :), 3);
Kf = qe_kinematic_prefactor(q1, dq = 0.0005, form = '2d', h_perp = h_perp, sigma_probe = sigma_probe, beam_kV = 30);
K = Kf(Ew); z = y ./ K;
n = round(100 / median(diff(Ew))); ys = movmean(y, n); zs = movmean(z, n);
[~, iy] = max(ys); [~, iz] = max(zs);
f = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 17 11], 'Color', 'w');
if isprop(f, 'Theme'), f.Theme = 'light'; end
ax = axes(f); hold(ax, 'on');
plot(ax, e, y / ys(iy), '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 4, 'HandleVisibility', 'off');
plot(ax, e, z / zs(iz), '.', 'Color', [0.95 0.75 0.65], 'MarkerSize', 4, 'HandleVisibility', 'off');
plot(ax, e, ys / ys(iy), 'k-', 'LineWidth', 1.8, 'DisplayName', 'measured counts');
plot(ax, e, zs / zs(iz), '-', 'Color', [0.85 0.30 0.10], 'LineWidth', 1.8, 'DisplayName', 'counts \div kinematic factor');
plot(ax, e, K / K(1), ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.3, 'DisplayName', 'kinematic factor K(E) / K(0.25 eV)');
xline(ax, e(iy), 'k:', sprintf('%.2f eV', e(iy)), 'LineWidth', 1.2, 'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
xline(ax, e(iz), ':', sprintf('%.2f eV', e(iz)), 'Color', [0.85 0.30 0.10], 'LineWidth', 1.2, ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
set(ax, 'FontName', 'Helvetica', 'FontSize', 10, 'TickDir', 'out', 'Box', 'on'); grid(ax, 'on');
xlim(ax, [0.25 1.8]); ylim(ax, [0 1.1]); xlabel(ax, 'energy loss (eV)'); ylabel(ax, 'normalised');
title(ax, sprintf('590, q = %+.4f Å^{-1}: counts before and after dividing by K (no fit)', q1), 'FontWeight', 'normal', 'FontSize', 10);
legend(ax, 'Location', 'southwest', 'FontSize', 9, 'Box', 'off');
tag = strrep(sprintf('%+.4f', q1), '.', 'p'); out = fullfile(fit_dir, ['fig_single_divide_q' tag]);
exportgraphics(f, [out '.png'], 'Resolution', 200, 'BackgroundColor', 'white');
exportgraphics(f, [out '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
close(f);
fprintf('maximum: counts %.0f meV, counts/K %.0f meV; K(1.8 eV)/K(0.25 eV) = %.3f\n', Ew(iy), Ew(iz), K(end) / K(1));
end
