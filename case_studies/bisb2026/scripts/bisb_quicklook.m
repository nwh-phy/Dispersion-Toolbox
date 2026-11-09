function out = bisb_quicklook(id, opts)
%BISB_QUICKLOOK  Standard first figures of a registered dataset in one call.
%   out = bisb_quicklook(id) writes to paper_results/quicklook/<id>_<time>/:
%     heatmap_q015   column-normalised q-E map, |q| <= 0.015, 200-1700 meV
%     heatmap_q005   the same for |q| <= 0.005
%     stack          (+q,-q) mean spectra at |q| = 0, 0.001, ..., 0.015 (44 meV
%                    moving mean, each divided by its mean in 400-1700 meV, offset 1)
%     q0_spectrum    q = 0 channel: full spectrum (log) and 150-1800 meV with the
%                    power-law tail of every bg window subtracted
%     summary.txt    alignment, ZLP centre and FWHM at q = 0, peak - ZLP from the
%                    raw counts and after each tail subtraction
%   A new acquisition needs only an entry in bisb_datasets. Raw counts, no
%   kinematic factor. The ~130 meV feature is the phonon of amorphous carbon
%   built up during acquisition: keep bg windows above it (default 150-250 and
%   200-300 meV); a tail fitted below 150 meV runs above the data.
%
%   See also: bisb_frame_map, qe_plot_heatmap, qe_pm_spectra, qe_peak_apex, qe_powerlaw_tail
arguments
    id (1,:) char
    opts.bg_windows (:,2) double = [150 250; 200 300]
    opts.out_dir (1,:) char = ''
end
S = bisb_frame_map(id); d = S.dataset; E = S.E(:); q = S.q;
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
out_dir = opts.out_dir;
if isempty(out_dir)
    out_dir = fullfile(root, 'paper_results', 'quicklook', ...
        sprintf('%s_%s', d.id, char(datetime('now', 'Format', 'yyyyMMdd''T''HHmmss'))));
end
if ~isfolder(out_dir), mkdir(out_dir); end
ttl = sprintf('%s: %s, %g μm defocus, %d × %g s', d.id, d.sample, d.defocus_um, d.n_frames, d.frame_s);

% heatmaps
h1 = qe_plot_heatmap(E, q, S.M, fullfile(out_dir, 'heatmap_q015'), q_window=[-0.015 0.015]);
h2 = qe_plot_heatmap(E, q, S.M, fullfile(out_dir, 'heatmap_q005'), q_window=[-0.005 0.005]);

% stacked spectra
qa = 0:0.001:0.015; [P, ok] = qe_pm_spectra(S.M, q, qa); qa = qa(ok); P = P(:, ok);
for i = 1:numel(qa)
    y = movmean(P(:, i), 11); P(:, i) = y / mean(y(E >= 400 & E <= 1700));
end
st = qe_prepare_stacked_spectra(E, P, qa, struct('q_start', 0, 'q_end', 0.015, 'q_step', 0.001, 'offset', 1.0));
fig = local_figure([50 50 800 1000]); ax = axes(fig); hold(ax, 'on');
w = E >= 150 & E <= 1800; c = turbo(numel(st.q_values) + 2); iend = find(E <= 1800, 1, 'last');
for i = 1:numel(st.q_values)
    plot(ax, E(w), st.shifted_traces(w, i), 'Color', c(i + 1, :), 'LineWidth', 1.3);
    text(ax, 1810, st.offsets(i) + st.traces(iend, i), sprintf('%.3f', st.q_values(i)), 'FontSize', 10, 'Color', c(i + 1, :) * 0.8);
end
xlim(ax, [150 1950]); ylim(ax, [0, st.offsets(end) + 2.5]); box(ax, 'on'); set(ax, 'FontSize', 13, 'TickDir', 'out');
xlabel(ax, 'Energy loss (meV)'); ylabel(ax, 'counts / mean(400–1700 meV) + offset');
title(ax, [ttl ', (+q,−q) mean'], 'FontWeight', 'normal', 'FontSize', 11, 'Interpreter', 'none');
text(ax, 1810, st.offsets(end) + 2.2, '|q| (Å^{-1})', 'FontSize', 10);
local_export(fig, fullfile(out_dir, 'stack'));

% q = 0 spectrum, tail subtraction and peak - ZLP
y = qe_pm_spectra(S.M, q, 0); [zc, zfw] = qe_zlp_center(y, E);
raw_apex = qe_peak_apex(E, y, [250 1000]) - zc;
nb = size(opts.bg_windows, 1); bg_apex = nan(1, nb); bg_r = nan(1, nb); bg_neg = nan(1, nb); sub = zeros(numel(E), nb);
for k = 1:nb
    [bg, bg_r(k)] = qe_powerlaw_tail(E, y, opts.bg_windows(k, :)); sub(:, k) = y - bg;
    bg_apex(k) = qe_peak_apex(E, sub(:, k), [250 1000]) - zc;
    m = E > opts.bg_windows(k, 2) & E <= 1000; bg_neg(k) = min(movmean(sub(m, k), 11)) / max(movmean(sub(m, k), 11));
end
fig = local_figure([50 50 1400 560]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl); m = E >= -200 & E <= 1800;
semilogy(ax, E(m) - zc, y(m), 'Color', [0.80 0.15 0.10], 'LineWidth', 1.1); grid(ax, 'on');
xlim(ax, [-200 1800]); xlabel(ax, 'Energy loss (meV)'); ylabel(ax, 'counts');
title(ax, [ttl ', q = 0'], 'FontWeight', 'normal', 'FontSize', 11, 'Interpreter', 'none');
ax = nexttile(tl); hold(ax, 'on'); m = E >= 150 & E <= 1800; cc = lines(nb);
for k = 1:nb
    plot(ax, E(m) - zc, movmean(sub(m, k), 11), 'Color', cc(k, :), 'LineWidth', 1.5, 'DisplayName', ...
        sprintf('tail fitted %g–%g meV (r = %.2f), peak %.0f meV', opts.bg_windows(k, :), bg_r(k), bg_apex(k)));
    if isfinite(bg_apex(k)), xline(ax, bg_apex(k), ':', 'Color', cc(k, :), 'LineWidth', 1.2, 'HandleVisibility', 'off'); end
end
yline(ax, 0, 'k-', 'HandleVisibility', 'off'); xlim(ax, [150 1800]); grid(ax, 'on'); box(ax, 'on');
xlabel(ax, 'Energy loss (meV)'); ylabel(ax, 'counts − power-law tail (44 meV moving mean)');
legend(ax, 'Location', 'northeast'); title(ax, 'q = 0, ZLP tail subtracted', 'FontWeight', 'normal');
local_export(fig, fullfile(out_dir, 'q0_spectrum'));

% summary
fid = fopen(fullfile(out_dir, 'summary.txt'), 'w');
fprintf(fid, '%s\n%s\nraw file: %s\nframe map: %s\n', ttl, d.dq_status, d.path, S.file);
fprintf(fid, 'dq %.5g 1/A, q = 0 at native channel %d, h %g 1/A\n', d.dq, d.zero, d.h);
fprintf(fid, 'alignment: energy shifts %.0f..%+.0f meV, beam drift %.2f..%+.2f channels, sigma_probe %.3g 1/A\n', ...
    min(S.offsets_meV), max(S.offsets_meV), min(S.drift_Ainv)/d.dq, max(S.drift_Ainv)/d.dq, S.sigma_probe);
fprintf(fid, 'q = 0: ZLP centre %.1f meV, FWHM %.1f meV\n', zc, zfw);
fprintf(fid, 'q = 0 peak - ZLP, raw counts: %.0f meV\n', raw_apex);
for k = 1:nb
    fprintf(fid, 'q = 0 peak - ZLP, tail %g-%g meV (r = %.2f) subtracted: %.0f meV; lowest residual above the window %.2f of the peak\n', ...
        opts.bg_windows(k, :), bg_r(k), bg_apex(k), bg_neg(k));
end
fprintf(fid, 'heatmap colour limits: q015 [%.3g %.3g], q005 [%.3g %.3g]\n', h1.clim, h2.clim);
fclose(fid);
out = struct('folder', out_dir, 'zlp_center_meV', zc, 'zlp_fwhm_meV', zfw, 'raw_peak_minus_zlp_meV', raw_apex, ...
    'bg_windows', opts.bg_windows, 'bg_exponent', bg_r, 'bg_peak_minus_zlp_meV', bg_apex, 'bg_lowest_residual', bg_neg);
fprintf('bisb_quicklook %s -> %s\n', d.id, out_dir);
end

function fig = local_figure(position)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', position);
if isprop(fig, 'Theme'), fig.Theme = 'light'; end   % white export also in macOS dark mode
end

function local_export(fig, base)
exportgraphics(fig, [base '.png'], 'Resolution', 200);
exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
close(fig);
end
