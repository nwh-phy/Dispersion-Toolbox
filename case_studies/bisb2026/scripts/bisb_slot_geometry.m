function out = bisb_slot_geometry(result_dir)
%BISB_SLOT_GEOMETRY slot-v4 detector-geometry corrections of a saved slot run.
%   bisb_slot_geometry(result_dir) reads analysis.mat written by
%   bisb_slot_analyze (it needs the per-channel scan-row sums Fr) and corrects
%   three detector effects that v3 left in the data:
%   1. The elastic zero depends on the detector channel (curved elastic line).
%      It is measured per channel on the aligned sum, with the diagonal streak
%      through the direct beam masked, interpolated across the channels the
%      streak crosses and smoothed.
%   2. The q position of the loss signal moves with energy. The q apex of the
%      direct beam and of the +/-G pair, after step 1, is fitted with
%      c = c0 + k + (gamma + beta*k)*E: a shift plus a magnification change.
%   3. The channel axis is resampled on integer offsets k from the direct beam
%      and binned in odd-width bins centred on k = 0, so +q and -q are compared
%      at equal |q|. Time reversal gives S(q,E) = S(-q,E) for single
%      scattering; the paired +/-q difference is reported as a systematic error.
%   The acceptance test is that the background-subtracted spectra at the Bragg
%   spots, which are double-scattering replicas of the central spectrum, line up
%   with the central spectrum in energy.
%   Writes figures 11-17, CSV tables, geometry.mat and v4_summary.txt into
%   result_dir. q in 1/A uses cfg.geometry.dq_Ainv, a provisional calibration.
arguments
    result_dir (1,:) char {mustBeFolder}
end
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
a = load(fullfile(result_dir, 'analysis.mat'), 'cfg', 'E', 'q', 'M', 'Fr');
g = a.cfg.geometry; E = a.E(:); q = double(a.q(:)).'; M = double(a.M); nq = numel(q);

% Elastic anchors: direct beam, strongest maximum on each side, outer left spot.
el = sum(M(abs(E) <= 40, :), 1);
[~, i] = max(el); c0 = local_apex3(q, el, i);
iL = local_argmax(el, q < c0 - 20); iR = local_argmax(el, q > c0 + 20);
iO = local_argmax(el, q < q(iL) - 15);
anchors = [local_apex3(q, el, iL) c0 local_apex3(q, el, iR) local_apex3(q, el, iO)];
names = {'left_G', 'direct', 'right_G', 'outer_left'};

% 1. Per-channel elastic zero. The streak is fitted on the gain side (k < 0),
% where it is separated from the loss signal.
ks = g.streak_fit_k; er = nan(size(ks));
for t = 1:numel(ks)
    j = find(q == round(c0 + ks(t)));
    y = movmean(sum(M(:, j + (-2:2)), 2), 5);
    w = find(E >= g.streak_search_slope(2) * ks(t) & E <= g.streak_search_slope(1) * ks(t));
    [~, p] = max(y(w)); er(t) = E(w(p));
end
sp = polyfit(ks, er, 1); k = q - c0; es = polyval(sp, k);
zraw = nan(1, nq); fw = nan(1, nq); ok = false(1, nq);
for j = 1:nq
    allow = abs(E) <= g.zlp_window_meV;
    if abs(k(j)) > g.streak_core_k, allow = allow & abs(E - es(j)) > g.streak_mask_meV; end
    w = find(allow); [~, i0] = max(M(w, j)); i0 = w(i0); jj = i0 + (-2:2);
    if i0 < 3 || i0 > numel(E) - 2 || ~all(allow(jj)), continue, end
    b = [ones(5, 1) E(jj) E(jj).^2] \ log(max(M(jj, j), 1));
    if b(3) >= 0, continue, end
    zraw(j) = -b(2) / (2 * b(3)); ok(j) = true; fw(j) = local_fwhm(E, M(:, j), i0);
end
ok = ok & ~(abs(k) > g.streak_core_k & abs(es) < g.streak_exclusion_meV);
z = zraw; z(~ok) = interp1(q(ok), zraw(ok), q(~ok), 'linear');
z = movmedian(fillmissing(z, 'nearest'), g.zlp_smooth);

% 2. q apex against energy after the energy correction; fit shift + magnification.
% Off-axis anchors are re-measured on the corrected elastic line: their elastic
% peaks sit up to 40 meV from the v3 zero and a fixed window would clip them.
Mr = bisb_slot_regrid(E, q, M, z, c0, k);
elr = sum(Mr(abs(E) <= 40, :), 1, 'omitnan');
anchors([1 3 4]) = [local_apex3(q, elr, iL) local_apex3(q, elr, iR) local_apex3(q, elr, iO)];
bands = g.drift_bands; Eb = mean(bands, 2); A = nan(numel(Eb) + 1, 3);
for b = 0:numel(Eb)
    if b == 0, m = abs(E) <= 40; else, m = E >= bands(b, 1) & E <= bands(b, 2); end
    I = sum(Mr(m, :), 1, 'omitnan');
    for s = 1:3, A(b + 1, s) = local_apex5(q, I, anchors(s)); end
end
D = A(2:end, :) - A(1, :); K = repmat(anchors(1:3) - c0, numel(Eb), 1); EE = repmat(Eb / 1000, 1, 3);
X = [EE(:) EE(:) .* K(:)]; P = X \ D(:); gamma = P(1); beta = P(2); drift_rms = rms(D(:) - X * P);

% 3. Resample on k, bin symmetrically about the direct beam.
kg = ceil(g.channel_range(1) - c0):floor(g.channel_range(2) - c0);
Y = bisb_slot_regrid(E, q, a.Fr, z, c0, kg, gamma, beta);
okE = all(isfinite(reshape(Y, numel(E), [])), 2); E4 = E(okE); Y = Y(okE, :, :);
Y0 = bisb_slot_regrid(E, q, a.Fr, zeros(1, nq), c0, kg); Y0 = Y0(okE, :, :);   % v3 geometry on the same grid
nr = size(Y, 3); w = g.bin_width; h = (w - 1) / 2;
assert(mod(w, 2) == 1, 'bisb_slot_geometry:BinWidth', 'Odd bin width keeps a bin centred on q = 0');
kc = (ceil((kg(1) + h) / w):floor((kg(end) - h) / w)) * w; nb = numel(kc);
S4 = zeros(numel(E4), nb, nr);
for b = 1:nb, S4(:, b, :) = sum(Y(:, abs(kg - kc(b)) <= h, :), 2); end
B4 = sum(S4, 3); qA = kc * g.dq_Ainv;

% Paired +/-q test with a shared scan-row bootstrap. The direct-beam profile
% itself is skewed (more intensity at +k), so the test is whether I(+q)/I(-q)
% depends on energy, which would move peaks; an energy-independent ratio
% leaves the spectral shape, and every peak energy, unchanged.
rng(a.cfg.seed); pk = kc(kc > 0 & ismember(-kc, kc)); np = numel(pk); nboot = a.cfg.bootstrap;
v = nan(np, 9); lo = nan(np, 9); hi = nan(np, 9); sym = zeros(numel(E4), np + 1);
sym(:, 1) = B4(:, kc == 0); syst = zeros(np + 1, 1);
for p = 1:np
    L = squeeze(S4(:, kc == -pk(p), :)); R = squeeze(S4(:, kc == pk(p), :));
    v(p, :) = local_pair(E4, sum(L, 2), sum(R, 2), g); bt = nan(nboot, 9);
    for t = 1:nboot
        ix = randi(nr, 1, nr); bt(t, :) = local_pair(E4, sum(L(:, ix), 2), sum(R(:, ix), 2), g);
    end
    lo(p, :) = prctile(bt, 2.5, 1); hi(p, :) = prctile(bt, 97.5, 1);
    sym(:, p + 1) = (sum(L, 2) + sum(R, 2)) / 2;
    % Apex difference that the measured ratio slope alone produces on the symmetrized spectrum.
    tilt = exp(v(p, 9) / 100 * E4 / 2000);
    syst(p + 1) = abs(qe_peak_apex(E4, sym(:, p + 1) .* tilt, [400 800], smooth = 25, half_width = 80) - ...
        qe_peak_apex(E4, sym(:, p + 1) ./ tilt, [400 800], smooth = 25, half_width = 80)) / 2;
end
pass = abs(v(:, 9)) <= g.sym_slope_tol;
Tsym = table(pk(:), pk(:) * g.dq_Ainv, v(:, 8), lo(:, 8), hi(:, 8), v(:, 9), lo(:, 9), hi(:, 9), ...
    v(:, 1), v(:, 2), v(:, 3), lo(:, 3), hi(:, 3), v(:, 6), v(:, 7), lo(:, 7), hi(:, 7), syst(2:end), pass, ...
    'VariableNames', {'k_channel', 'abs_q_Ainv', 'intensity_plus_over_minus', 'intensity_lo', 'intensity_hi', ...
    'ratio_slope_pct_per_eV', 'slope_lo', 'slope_hi', 'low_apex_minus_q', 'low_apex_plus_q', ...
    'apex_plus_minus_minus', 'apex_diff_lo', 'apex_diff_hi', 'band_ratio_plus_over_minus', ...
    'shape_shift_meV', 'shift_lo', 'shift_hi', 'systematic_meV', 'pass'});
writetable(Tsym, fullfile(result_dir, 'symmetry_test.csv'));
symq = [0 pk] * g.dq_Ainv; sapex = nan(np + 1, 2);
for p = 1:np + 1
    sapex(p, 1) = qe_peak_apex(E4, sym(:, p), [400 800], smooth = 25, half_width = 80);
    sapex(p, 2) = qe_peak_apex(E4, sym(:, p), [800 1600], smooth = 25, half_width = 80);
end
Tsp = table([0 pk].', symq.', sapex(:, 1), sapex(:, 2), syst, [true; pass], ...
    'VariableNames', {'k_channel', 'abs_q_Ainv', 'low_apex_meV', 'high_apex_meV', 'systematic_meV', 'pass'});
writetable(Tsp, fullfile(result_dir, 'symmetrized_apex.csv'));

% Acceptance test: Bragg replicas against the central spectrum, v3 and v4 geometry.
% Every window is sampled at the same offsets -3..3 from its own centre.
Tr = table(); spots = anchors([1 3 4]); sname = names([1 3 4]); fig17 = cell(2, numel(spots));
for s = 1:numel(spots)
    for ver = 1:2
        ks = spots(s) - c0; kw = [(-3:3), ks + (-3:3), ks - 10 + (-3:3), ks + 10 + (-3:3)];
        if ver == 1, Z = bisb_slot_regrid(E, q, M, zeros(1, nq), c0, kw);
        else, Z = bisb_slot_regrid(E, q, M, z, c0, kw, gamma, beta); end
        Z = reshape(Z(okE, :), numel(E4), 7, 4); cen = sum(Z(:, :, 1), 2); spot = sum(Z(:, :, 2), 2);
        inside = c0 + ks + [-10 10] - 3 >= q(1) & c0 + ks + [-10 10] + 3 <= q(end);
        nbk = squeeze(sum(Z(:, :, 3:4), 2)); nbk = nbk(:, inside);
        rep = spot - mean(nbk, 2);
        [dlt, amp] = local_shift(E4, cen, rep, g.replica_window, g.shift_grid);
        elr = sum(rep(abs(E4) <= 40)) / sum(cen(abs(E4) <= 40));
        Tr = [Tr; table(sname(s), spots(s), ver + 2, dlt, amp, elr, amp / elr, size(nbk, 2), ...
            qe_peak_apex(E4, spot, [400 800], smooth = 25, half_width = 80), ...
            'VariableNames', {'spot', 'channel', 'version', 'shift_vs_centre_meV', 'amplitude', ...
            'elastic_ratio', 'inelastic_over_elastic_ratio', 'background_windows', 'spot_raw_low_apex_meV'})]; %#ok<AGROW>
        fig17{ver, s} = struct('cen', cen, 'rep', rep, 'amp', amp, 'dlt', dlt);
    end
end
writetable(Tr, fullfile(result_dir, 'replica_check.csv'));

% Calibrated maps: full v4 map, energy-only Savitzky-Golay displays, symmetrized map.
lab = sprintf('q (Å^{-1}); provisional dq = %.5f', g.dq_Ainv);
base = qe_plot_heatmap(E4, qA, B4, fullfile(result_dir, '13_v4_dispersion'), q_window = [min(qA) max(qA)], ...
    e_window = [200 1500], norm_window = [400 1500], q_label = lab, title = 'slot-v4: per-channel zero, q distortion, bins centred on q = 0');
qe_plot_heatmap(E4, qA, B4, fullfile(result_dir, '13_v4_dispersion_wide'), q_window = [min(qA) max(qA)], ...
    e_window = [200 3000], norm_window = [400 3000], q_label = lab, title = 'slot-v4, 0.2-3 eV');
for fr = g.sg_frames
    pp = qe_preprocess(struct('intensity', B4, 'energy_meV', E4), struct('do_denoise', true, ...
        'denoise_method', 'SavGol', 'sg_order', 3, 'sg_framelen', fr));
    qe_plot_heatmap(E4, qA, pp.intensity, fullfile(result_dir, sprintf('14_v4_sg%d', fr)), q_window = [min(qA) max(qA)], ...
        e_window = [200 1500], norm_window = [400 1500], clim = base.clim, q_label = lab, ...
        title = sprintf('slot-v4; energy-only Savitzky-Golay %d channels (display)', fr));
end
qe_plot_heatmap(E4, symq, sym, fullfile(result_dir, '16_v4_symmetrized'), q_window = [0 max(symq)], ...
    e_window = [200 1500], norm_window = [400 1500], q_label = sprintf('|q| (Å^{-1}); provisional dq = %.5f', g.dq_Ainv), ...
    title = '(I(+q) + I(-q)) / 2, slot-v4');

% Diagnostic figures.
f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 900]); tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
ax = nexttile(tl, [1 2]); m = E >= -400 & E <= 500;
imagesc(ax, q, E(m), log10(max(M(m, :), 1))); axis(ax, 'xy'); hold(ax, 'on'); colorbar(ax);
plot(ax, q(ok), zraw(ok), 'w.', q, z, 'k-', q, es, 'm--', 'MarkerSize', 6, 'LineWidth', 1);
xlabel(ax, 'Detector channel'); ylabel(ax, 'Energy after v3 alignment (meV)');
title(ax, 'log_{10} aligned sum; white: measured elastic zero; black: zero used; magenta: streak fit', 'FontWeight', 'normal');
ax = nexttile(tl); plot(ax, q, z, 'k-', q(ok), zraw(ok), '.', q(~ok), zraw(~ok), 'x'); grid(ax, 'on');
xline(ax, anchors, ':'); xlabel(ax, 'Detector channel'); ylabel(ax, 'Elastic zero (meV)'); legend(ax, 'used', 'measured', 'masked', 'Location', 'best');
ax = nexttile(tl); semilogy(ax, q, el); grid(ax, 'on'); xline(ax, anchors, ':'); xlabel(ax, 'Detector channel'); ylabel(ax, 'Counts within \pm40 meV');
exportgraphics(f, fullfile(result_dir, '11_elastic_zero.png'), 'Resolution', 180); close(f);

f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 900 600]); ax = axes(f); hold(ax, 'on');
cols = lines(3); ee = linspace(0, max(Eb), 50) / 1000;
for s = 1:3
    plot(ax, Eb / 1000, D(:, s), 'o', 'Color', cols(s, :), 'MarkerFaceColor', cols(s, :));
    plot(ax, ee, (gamma + beta * (anchors(s) - c0)) * ee, '-', 'Color', cols(s, :));
end
grid(ax, 'on'); box(ax, 'on'); xlabel(ax, 'Energy loss (eV)'); ylabel(ax, 'q apex minus elastic apex (channels)');
legend(ax, {'-G data', '-G model', 'direct data', 'direct model', '+G data', '+G model'}, 'Location', 'northwest');
title(ax, sprintf('c = c_0 + k + (%.3f %+.5f k) E;  fixed point k = %.0f;  rms %.3f ch', gamma, beta, -gamma / beta, drift_rms), 'FontWeight', 'normal');
exportgraphics(f, fullfile(result_dir, '12_q_distortion.png'), 'Resolution', 180); close(f);

% Same pair metrics in the v3 geometry, for comparison (point estimates).
B0 = zeros(numel(E4), nb); for b = 1:nb, B0(:, b) = sum(sum(Y0(:, abs(kg - kc(b)) <= h, :), 2), 3); end
v0 = nan(np, 9); for p = 1:np, v0(p, :) = local_pair(E4, B0(:, kc == -pk(p)), B0(:, kc == pk(p)), g); end
Tsym.ratio_slope_v3_geometry = v0(:, 9); writetable(Tsym, fullfile(result_dir, 'symmetry_test.csv'));
f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 900]); tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact');
ax = nexttile(tl); hold(ax, 'on'); plot(ax, pk * g.dq_Ainv, v0(:, 9), 'x', 'Color', [.6 .6 .6]);
errorbar(ax, pk * g.dq_Ainv, v(:, 9), v(:, 9) - lo(:, 9), hi(:, 9) - v(:, 9), 'o'); yline(ax, 0, ':');
yline(ax, g.sym_slope_tol * [-1 1], '--'); grid(ax, 'on'); box(ax, 'on'); legend(ax, 'v3 geometry', 'v4', 'Location', 'best');
xlabel(ax, '|q| (Å^{-1})'); ylabel(ax, 'd ln(I(+q)/I(-q)) / dE  (% per eV)');
title(ax, sprintf('Energy dependence of the +/-q ratio, %g-%g meV', g.shape_window), 'FontWeight', 'normal');
ax = nexttile(tl); plot(ax, pk * g.dq_Ainv, v(:, 1), 'o', pk * g.dq_Ainv, v(:, 2), 's', symq, sapex(:, 1), 'k.-'); grid(ax, 'on');
xlabel(ax, '|q| (Å^{-1})'); ylabel(ax, 'Low-band apex, 400-800 meV search (meV)'); legend(ax, '-q', '+q', 'symmetrized', 'Location', 'best');
ax = nexttile(tl); errorbar(ax, pk * g.dq_Ainv, v(:, 8), v(:, 8) - lo(:, 8), hi(:, 8) - v(:, 8), 'o'); yline(ax, 1, ':');
grid(ax, 'on'); xlabel(ax, '|q| (Å^{-1})'); ylabel(ax, sprintf('I(+q) / I(-q), %g-%g meV', g.shape_window));
title(ax, 'Energy-independent part: the direct-beam profile is itself skewed', 'FontWeight', 'normal');
ax = nexttile(tl); Lm = fliplr(B4(:, ismember(kc, -pk))); Rm = B4(:, ismember(kc, pk)); m = E4 >= 200 & E4 <= 1600;
lr = log(Rm(m, :) ./ Lm(m, :)); ms = E4(m) >= g.shape_window(1); lr = lr - mean(lr(ms, :), 1);
imagesc(ax, pk * g.dq_Ainv, E4(m), lr); axis(ax, 'xy'); colormap(ax, redblue_local()); clim(ax, [-0.2 0.2]); colorbar(ax);
xlabel(ax, '|q| (Å^{-1})'); ylabel(ax, 'Energy loss (meV)');
title(ax, sprintf('ln(I(+q)/I(-q)) minus its mean over %g-%g meV', g.shape_window), 'FontWeight', 'normal');
exportgraphics(f, fullfile(result_dir, '15_symmetry_test.png'), 'Resolution', 180); close(f);

f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1500 800]); tl = tiledlayout(f, 2, numel(spots), 'TileSpacing', 'compact');
for ver = 1:2
    for s = 1:numel(spots)
        ax = nexttile(tl); r = fig17{ver, s}; mm = E4 >= 200 & E4 <= 2400;
        plot(ax, E4(mm) / 1000, movmean(r.rep(mm), 9), 'LineWidth', 1.1); hold(ax, 'on');
        plot(ax, E4(mm) / 1000, r.amp * movmean(r.cen(mm), 9), 'k--', 'LineWidth', 1);
        grid(ax, 'on'); xlabel(ax, 'Energy loss (eV)'); ylabel(ax, 'Counts');
        title(ax, sprintf('v%d %s: shift %+.1f meV', ver + 2, strrep(sname{s}, '_', ' '), r.dlt), 'FontWeight', 'normal');
    end
end
legend(ax, 'spot minus neighbours', 'scaled centre', 'Location', 'best');
exportgraphics(f, fullfile(result_dir, '17_replica_check.png'), 'Resolution', 180); close(f);

% MoS2 cross-check of the momentum scale, if the outer left spot is MoS2(100).
G_mos2 = 4 * pi / (sqrt(3) * g.mos2_a_A); dq_mos2 = G_mos2 / (c0 - anchors(4));
Tz = table(q.', k.', zraw.', z.', ok.', fw.', el.', 'VariableNames', ...
    {'channel', 'k_channel', 'zero_measured_meV', 'zero_used_meV', 'measured_valid', 'zlp_fwhm_meV', 'elastic_counts'});
writetable(Tz, fullfile(result_dir, 'elastic_zero.csv'));
Td = table(Eb, D(:, 1), D(:, 2), D(:, 3), 'VariableNames', {'band_centre_meV', 'shift_left_G', 'shift_direct', 'shift_right_G'});
writetable(Td, fullfile(result_dir, 'q_distortion.csv'));
geometry = struct('c0', c0, 'anchors', anchors, 'anchor_names', {names}, 'streak_fit', sp, 'zero', z, ...
    'gamma_ch_per_eV', gamma, 'beta_per_eV', beta, 'fixed_point_k', -gamma / beta, 'drift_rms_ch', drift_rms, ...
    'dq_Ainv', g.dq_Ainv, 'dq_source', g.dq_source, 'dq_mos2_Ainv', dq_mos2);
save(fullfile(result_dir, 'geometry.mat'), 'geometry', 'E4', 'kg', 'kc', 'qA', 'S4', 'B4', 'sym', 'symq', ...
    'Tsym', 'Tsp', 'Tr', 'Tz', 'Td', '-v7.3');

fid = fopen(fullfile(result_dir, 'v4_summary.txt'), 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'slot-v4 geometry corrections (%s)\n', a.cfg.version);
fprintf(fid, 'Figures 01-10 and their tables in this folder keep the v3 geometry (same alignment); 11-17 are slot-v4.\n');
fprintf(fid, 'Elastic anchors (channels): -G %.3f, direct %.3f, +G %.3f, outer left %.3f\n', anchors([1 2 3 4]));
fprintf(fid, 'Streak: E = %.3f meV/ch * k %+.1f meV (gain-side fit); masked within %g meV\n', sp, g.streak_mask_meV);
fprintf(fid, 'Elastic zero used: %.1f to %.1f meV over channels %d-%d; %d/%d channels measured directly\n', ...
    min(z(q >= g.channel_range(1) & q <= g.channel_range(2))), max(z(q >= g.channel_range(1) & q <= g.channel_range(2))), g.channel_range, nnz(ok), nq);
fprintf(fid, 'q distortion: gamma = %.4f ch/eV, beta = %.6f /eV, fixed point k = %.1f (channel %.1f); rms residual %.3f ch over %d points\n', ...
    gamma, beta, -gamma / beta, c0 - gamma / beta, drift_rms, numel(D));
fprintf(fid, 'Grid k = %d..%d; %d-channel bins centred on k = 0; energy support %.0f..%.0f meV\n', kg([1 end]), w, E4([1 end]));
fprintf(fid, 'Momentum: dq = %.10f 1/A (%s); MoS2(100) cross-check with a = %.3f A: dq = %.10f 1/A (%+.2f%%)\n', ...
    g.dq_Ainv, g.dq_source, g.mos2_a_A, dq_mos2, 100 * (dq_mos2 / g.dq_Ainv - 1));
fprintf(fid, '+/-q pairs with |d ln(I+/I-)/dE| <= %g %%/eV over %g-%g meV: %d/%d (v3 geometry: %d/%d); median |slope| %.1f %%/eV (v3 %.1f)\n', ...
    g.sym_slope_tol, g.shape_window, nnz(pass), np, nnz(abs(v0(:, 9)) <= g.sym_slope_tol), np, median(abs(v(:, 9))), median(abs(v0(:, 9))));
fprintf(fid, 'Systematic apex error from the slope (half +/-q difference on the symmetrized spectrum): median %.1f meV, max %.1f meV\n', ...
    median(syst(2:end), 'omitnan'), max(syst(2:end)));
elk = sum(sum(Y(abs(E4) <= 40, :, :), 3), 1);
fprintf(fid, 'Elastic I(+k)/I(-k) at k = 3, 5, 7 channels: %.2f, %.2f, %.2f (skewed direct-beam profile)\n', ...
    elk(kg == 3) / elk(kg == -3), elk(kg == 5) / elk(kg == -5), elk(kg == 7) / elk(kg == -7));
fprintf(fid, 'Replica shifts vs centre (meV), v3 -> v4:\n');
for s = 1:numel(spots)
    r3 = Tr(strcmp(Tr.spot, sname{s}) & Tr.version == 3, :); r4 = Tr(strcmp(Tr.spot, sname{s}) & Tr.version == 4, :);
    fprintf(fid, '  %s: %+.1f -> %+.1f; inelastic/elastic ratio %.3f -> %.3f\n', sname{s}, r3.shift_vs_centre_meV, ...
        r4.shift_vs_centre_meV, r3.inelastic_over_elastic_ratio, r4.inelastic_over_elastic_ratio);
end
fprintf(fid, 'Bootstrap: %d draws over %d scan rows, shared between +q and -q. Not independent acquisitions.\n', nboot, nr);
disp(Tsym); disp(Tr); fprintf('%s', fileread(fullfile(result_dir, 'v4_summary.txt')));
out = geometry;
end

function i = local_argmax(y, mask)
w = find(mask); [~, p] = max(y(w)); i = w(p);
end

function c = local_apex3(q, y, i)
ly = log(max(y(i-1:i+1), 1)); c = q(i) + (ly(3) - ly(1)) / (2 * (2 * ly(2) - ly(1) - ly(3))) * (q(i+1) - q(i));
end

function c = local_apex5(q, I, c_guess)
jj = find(abs(q - c_guess) <= 8); [~, i0] = max(I(jj)); i0 = jj(i0); kk = i0 + (-2:2);
b = [ones(5, 1) q(kk).' q(kk).'.^2] \ log(I(kk)).'; c = -b(2) / (2 * b(3));
end

function f = local_fwhm(E, y, i0)
h = y(i0) / 2; l = find(y(1:i0) < h, 1, 'last'); r = i0 - 1 + find(y(i0:end) < h, 1, 'first'); f = NaN;
if ~isempty(l) && ~isempty(r)
    f = interp1(y([r-1 r]), E([r-1 r]), h) - interp1(y([l l+1]), E([l l+1]), h);
end
end

function v = local_pair(E, l, r, g)
% [apex -q, apex +q, apex difference, band ratio -q, band ratio +q, band ratio +/-,
%  shape shift (NaN at the grid edge), intensity +/-, slope of ln(I+/I-) in % per eV]
al = qe_peak_apex(E, l, [400 800], smooth = 25, half_width = 80);
ar = qe_peak_apex(E, r, [400 800], smooth = 25, half_width = 80);
lo = E >= 400 & E <= 800; hi = E >= 1200 & E <= 1600; m = E >= g.shape_window(1) & E <= g.shape_window(2);
rl = trapz(E(lo), l(lo)) / trapz(E(hi), l(hi)); rr = trapz(E(lo), r(lo)) / trapz(E(hi), r(hi));
d = local_shift(E, l, r, g.shape_window, g.shift_grid);
if abs(d) >= max(abs(g.shift_grid)) - 1, d = NaN; end
% ln(I+/I-) against energy on 100 meV blocks, so noise in single channels does not dominate.
ed = g.shape_window(1):100:g.shape_window(2); x = nan(numel(ed) - 1, 1); y = x;
for b = 1:numel(ed) - 1
    mb = E >= ed(b) & E < ed(b + 1); x(b) = mean(E(mb)) / 1000; y(b) = log(sum(r(mb)) / sum(l(mb)));
end
ok = isfinite(y); p = polyfit(x(ok), y(ok), 1);
v = [al ar ar - al rl rr rr / rl d trapz(E(m), r(m)) / trapz(E(m), l(m)) 100 * p(1)];
end

function [d, A] = local_shift(E, ref, y, window, grid)
% Energy shift d and amplitude A minimising sum (y(E) - A ref(E - d))^2 in window;
% d > 0 means the features of y lie at higher energy than those of ref.
m = E >= window(1) & E <= window(2); res = inf(size(grid)); amp = nan(size(grid));
for t = 1:numel(grid)
    x = interp1(E, ref, E(m) - grid(t), 'linear', NaN);
    if any(~isfinite(x)), continue, end
    amp(t) = (x' * y(m)) / (x' * x); res(t) = sum((y(m) - amp(t) * x).^2);
end
[~, t] = min(res); d = grid(t); A = amp(t);
if t > 1 && t < numel(grid) && all(isfinite(res(t-1:t+1)))
    den = res(t-1) - 2 * res(t) + res(t+1);
    if den > 0, d = grid(t) + 0.5 * (res(t-1) - res(t+1)) / den * (grid(2) - grid(1)); end
end
end

function c = redblue_local()
x = linspace(0, 1, 128)'; c = [[x; ones(128, 1)] [x; flipud(x)] [ones(128, 1); flipud(x)]];
end
