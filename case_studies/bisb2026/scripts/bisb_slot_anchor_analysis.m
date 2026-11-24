function out = bisb_slot_anchor_analysis(result_dir)
%BISB_SLOT_ANCHOR_ANALYSIS Bragg-anchored momentum and robust regional metrics.
%   Uses the direct beam and the strongest symmetric off-axis elastic pair.
%   q/G_anchor is dimensionless: the pair is at -1/+1 by assignment. No
%   diffraction order, lattice constant or absolute reciprocal-space calibration.
arguments
    result_dir (1,:) char {mustBeFolder}
end
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root, 'src')));
a = load(fullfile(result_dir, 'analysis.mat'));
E = a.E(:); q = a.q(:).'; elastic = a.elastic(:).';
[tf, prominence] = islocalmax(elastic, 'MinProminence', max(elastic) * 0.0005);
[~, direct] = max(elastic);
left = find(tf & q < q(direct) - 20); right = find(tf & q > q(direct) + 20);
assert(~isempty(left) && ~isempty(right), 'bisb_slot_anchor_analysis:Anchors', 'Missing off-axis elastic anchors');
[~, il] = max(prominence(left)); [~, ir] = max(prominence(right));
idx = [left(il) direct right(ir)];
centres = arrayfun(@(k) local_log_parabola(q, elastic, k), idx); centres = centres(:);
assigned = [-1; 0; 1];
fit = [assigned ones(3, 1)] \ centres;
scale = fit(1); zero = fit(2); qn = (a.qb - zero) / scale;
A = table({'left'; 'direct'; 'right'}, q(idx).', centres, assigned, centres - (zero + scale * assigned), ...
    'VariableNames', {'anchor', 'integer_channel', 'fitted_channel', 'assigned_q_over_G', 'residual_channel'});
writetable(A, fullfile(result_dir, 'elastic_anchors.csv'));

nd = ceil(size(a.B, 2) / 2); Bd = zeros(numel(E), nd); qd = zeros(1, nd);
for k = 1:nd
    jj = (2*k-1):min(2*k, size(a.B, 2)); Bd(:, k) = sum(a.B(:, jj), 2); qd(k) = mean(qn(jj));
end
qe_plot_heatmap(E, qd, movmean(Bd, 25, 1), fullfile(result_dir, '08_anchor_dispersion'), ...
    q_window=[min(qd) max(qd)], e_window=[200 1800], norm_window=[400 1600], ...
    q_label='q / G_{anchor}', title='16-channel q bins; 100 meV display smoothing; pair at +/-1');

ranges = [65 80; 89 104; 225 232; 233 240; 369 376; 121 200; 265 344];
labels = {'left_outer'; 'left_spot'; 'central_left'; 'central_right'; 'right_spot'; 'between_left'; 'between_right'};
K = numel(labels); values = nan(K, 4); lo = nan(K, 4); hi = nan(K, 4); success = nan(K, 2); qc = nan(K, 1);
rng(a.cfg.seed);
for j = 1:K
    take = a.qb >= ranges(j, 1) & a.qb <= ranges(j, 2); qc(j) = mean(qn(take));
    X = squeeze(sum(a.S(:, take, :), 2)); y = sum(X, 2); values(j, :) = local_metrics(E, y, a.cfg);
    rows = squeeze(sum(reshape(X, numel(E), 19, 19), 2)); boot = nan(a.cfg.bootstrap, 4);
    for k = 1:a.cfg.bootstrap
        yy = sum(rows(:, randi(19, 1, 19)), 2); boot(k, :) = local_metrics(E, yy, a.cfg);
    end
    for m = 1:4
        ci = prctile(boot(:, m), [2.5 97.5]); lo(j, m) = ci(1); hi(j, m) = ci(2);
    end
    success(j, :) = mean(isfinite(boot(:, 1:2)), 1);
end
T = table(labels, ranges(:, 1), ranges(:, 2), qc, ...
    values(:, 1), lo(:, 1), hi(:, 1), values(:, 2), lo(:, 2), hi(:, 2), ...
    values(:, 3), lo(:, 3), hi(:, 3), values(:, 4), lo(:, 4), hi(:, 4), success(:, 1), success(:, 2), ...
    'VariableNames', {'region', 'channel_first', 'channel_last', 'q_over_G', ...
    'low_apex_meV', 'low_lo', 'low_hi', 'high_apex_meV', 'high_lo', 'high_hi', ...
    'centroid_meV', 'centroid_lo', 'centroid_hi', 'low_high_ratio', 'ratio_lo', 'ratio_hi', ...
    'low_success', 'high_success'});
writetable(T, fullfile(result_dir, 'region_metrics.csv'));

f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1200 820]); tl = tiledlayout(f, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
names = {'Low-band apex (0.4–0.8 eV search)', 'High-band apex (0.8–1.6 eV search)', ...
    'Spectral centroid (0.4–1.6 eV)', 'I(0.4–0.8) / I(1.2–1.6)'};
for m = 1:4
    ax = nexttile(tl); hold(ax, 'on');
    errorbar(ax, qc, values(:, m), values(:, m) - lo(:, m), hi(:, m) - values(:, m), 'o', 'LineWidth', 1.1, 'CapSize', 8);
    for u = [-1 0 1], xline(ax, u, ':', 'HandleVisibility', 'off'); end
    grid(ax, 'on'); box(ax, 'on'); xlabel(ax, 'q / G_{anchor}'); title(ax, names{m}, 'FontWeight', 'normal');
    if m < 4, ylabel(ax, 'Energy (meV)'); else, ylabel(ax, 'Integrated intensity ratio'); end
end
exportgraphics(f, fullfile(result_dir, '09_region_metrics.png'), 'Resolution', 200); close(f);

f = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1000 560]); ax = axes(f); hold(ax, 'on');
semilogy(ax, (q - zero) / scale, elastic, 'LineWidth', 1.2);
plot(ax, assigned, elastic(idx), 'o', 'MarkerSize', 9, 'LineWidth', 1.5);
set(ax, 'YScale', 'log'); ylim(ax, [1e3 2e9]);
for k = 1:3, text(ax, assigned(k), elastic(idx(k)) * 1.3, sprintf('%.2f ch', centres(k)), 'HorizontalAlignment', 'center'); end
grid(ax, 'on'); xlabel(ax, 'q / G_{anchor}'); ylabel(ax, 'Counts within \pm40 meV'); title(ax, 'Elastic anchors: direct beam and strongest off-axis pair', 'FontWeight', 'normal');
exportgraphics(f, fullfile(result_dir, '10_elastic_anchors.png'), 'Resolution', 200); close(f);

fid = fopen(fullfile(result_dir, 'anchor_summary.txt'), 'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Direct channel %.4f; left %.4f; right %.4f\n', centres(2), centres(1), centres(3));
fprintf(fid, 'Pair separation %.4f channels; half-separation %.4f channels\n', centres(3) - centres(1), (centres(3) - centres(1)) / 2);
fprintf(fid, 'Three-anchor fit: channel = %.6f + %.6f * (q/G_anchor)\n', zero, scale);
fprintf(fid, 'Pair midpoint minus direct = %.4f channels\n', (centres(1) + centres(3)) / 2 - centres(2));
fprintf(fid, 'If the off-axis anchors are -G/+G, |G| is one normalized momentum unit and dq = 2|G| / pair separation. No diffraction order or |G| assumed here.\n');
fprintf(fid, 'Peak windows are separate low/high bands with 100 meV smoothing. Intervals resample 19 spatial rows (200 draws); they are not independent-acquisition errors.\n');
disp(T);
out = struct('anchors', A, 'regions', T, 'zero_channel', zero, 'scale_channel', scale);
end

function c = local_log_parabola(q, y, i)
ly = log(max(y(i-1:i+1), 1)); dq = q(i+1) - q(i);
d = (ly(3) - ly(1)) / (2 * (2 * ly(2) - ly(1) - ly(3))); c = q(i) + d * dq;
end

function v = local_metrics(E, y, cfg)
v = [qe_peak_apex(E, y, [400 800], smooth=25, half_width=80), ...
    qe_peak_apex(E, y, [800 1600], smooth=25, half_width=80), NaN, NaN];
m = E >= 400 & E <= 1600;
v(3) = trapz(E(m), y(m) .* E(m)) / trapz(E(m), y(m));
l = E >= 400 & E <= 800; h = E >= 1200 & E <= 1600;
v(4) = trapz(E(l), y(l)) / trapz(E(h), y(h));
end
