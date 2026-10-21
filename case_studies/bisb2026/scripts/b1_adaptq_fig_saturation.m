function out = b1_adaptq_fig_saturation(run_root, session)
% Measured E0 of both B1 branches (2D, h = 0.0005 fits, stage5/e0_stage2_2d_h5_
% <session>.csv) against the quasi-2D plasmon E(q) = Einf*sqrt(rho*q/(1+rho*q))
% for several interband screening lengths rho (illustration, not a fit), with
% the light cone and q < q_E (30 kV) at Einf shaded. Points with |q| < 0.0025
% (q < q_E) are drawn open. Writes stage5/figures_stage2_2d_h5/
% dispersion_saturation_<session>.png.
if nargin < 2, session = '590_PL2_10w'; end
T = readtable(fullfile(run_root, 'stage5', sprintf('e0_stage2_2d_h5_%s.csv', session)));
q = T.q; aq = abs(q); ok = strcmp(T.status, 'converged') & aq <= 0.0231;
hc = 1973.27;                       % eV*A
beta = 0.3284;                      % v/c at 30 kV
Einf = 1.60;                        % eV, upper-branch plateau
qL = Einf / hc; qE = Einf / (hc * beta);
f = figure('Visible', 'off', 'Units', 'centimeters', 'Position', [2 2 17 11], 'Color', 'w');
if isprop(f, 'Theme'), f.Theme = 'light'; end
ax = axes(f); hold(ax, 'on'); box(ax, 'on');
patch(ax, [0 qL qL 0], [0 0 2.1 2.1], [0.80 0.80 0.80], 'EdgeColor', 'none', 'HandleVisibility', 'off');
patch(ax, [qL qE qE qL], [0 0 2.1 2.1], [0.92 0.92 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
text(ax, qL / 2, 0.06, 'light cone', 'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7.5, 'Color', [0.25 0.25 0.25]);
text(ax, (qL + qE) / 2, 0.06, 'q < q_E (30 kV)', 'Rotation', 90, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'FontSize', 7.5, 'Color', [0.25 0.25 0.25]);
qq = linspace(1e-6, 0.024, 600);
rhos = [100 300 1000 3000]; dy = [0 -0.015 0 0.045]; cm = [0.75 0.75 0.75; 0.55 0.55 0.55; 0.35 0.35 0.35; 0.10 0.10 0.10];
for i = 1:numel(rhos)
    E = Einf * sqrt(rhos(i) * qq ./ (1 + rhos(i) * qq));
    plot(ax, qq, E, '-', 'Color', cm(i, :), 'LineWidth', 1.1, 'HandleVisibility', 'off');
    text(ax, 0.0241, E(end) + dy(i), sprintf(' \\rho_0 = %d Å', rhos(i)), 'FontSize', 8, 'Color', cm(i, :), 'VerticalAlignment', 'middle');
end
rel = aq >= 0.0025;
cU = [0.80 0.25 0.10]; cL = [0.10 0.35 0.75];
m = ok & rel;
errorbar(ax, aq(m), T.E0_high(m) / 1e3, T.E0_high_sd(m) / 1e3, 'o', 'Color', cU, 'MarkerFaceColor', cU, 'MarkerSize', 4, 'CapSize', 2, 'DisplayName', 'upper plasmon E_0');
errorbar(ax, aq(m), T.E0_low(m) / 1e3, T.E0_low_sd(m) / 1e3, 'o', 'Color', cL, 'MarkerFaceColor', cL, 'MarkerSize', 4, 'CapSize', 2, 'DisplayName', 'lower plasmon E_0');
u = ok & ~rel;
plot(ax, aq(u), T.E0_high(u) / 1e3, 'o', 'Color', cU, 'MarkerSize', 4, 'DisplayName', 'E_0 at q < q_E (fit unreliable)');
plot(ax, aq(u), T.E0_low(u) / 1e3, 'o', 'Color', cL, 'MarkerSize', 4, 'HandleVisibility', 'off');
a = ok & rel & T.apex_measured_LF < 1790;
plot(ax, aq(a), T.apex_measured_LF(a) / 1e3, 'x', 'Color', [0.2 0.2 0.2], 'MarkerSize', 5, 'DisplayName', 'peak of measured loss function');
xlim(ax, [0 0.0292]); ylim(ax, [0 2.0]);
set(ax, 'XTick', 0:0.005:0.025, 'FontName', 'Helvetica', 'FontSize', 9, 'TickDir', 'in', 'LineWidth', 0.6);
xlabel(ax, '|\itq\rm| (Å^{-1})'); ylabel(ax, 'Energy (eV)');
legend(ax, 'Location', 'southeast', 'FontSize', 8, 'Box', 'off');
title(ax, 'Grey lines: E(q) = E_\infty [\rho_0 q / (1 + \rho_0 q)]^{1/2},  E_\infty = 1.6 eV', 'FontSize', 9, 'FontWeight', 'normal');
outd = fullfile(run_root, 'stage5', 'figures_stage2_2d_h5'); if ~isfolder(outd), mkdir(outd); end
out = fullfile(outd, sprintf('dispersion_saturation_%s.png', session));
exportgraphics(f, out, 'Resolution', 250, 'BackgroundColor', 'white');
close(f); fprintf('qL = %.5f, qE = %.5f; saved %s\n', qL, qE, out);
end
