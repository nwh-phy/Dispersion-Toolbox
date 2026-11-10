function out = b1_whole_band_check(run_root, session, stage)
% Whole-B1 band without the two-branch decomposition (2026-10-19).
% Measured loss function per high-|q| bin of the stored fits:
%   (observed - fitted ZLP tail and low-energy losses) / kinematic prefactor,
% 100 meV moving mean, 250-1800 meV. Prints apex, half-maximum points and
% LF(300 meV)/apex per bin. Then assumes NO dispersion: the mean normalised
% loss function of |q| = 0.0085-0.0155 times the point-Q 2D kinematic factor,
% and prints the maximum at small Q (what a peak tracker would see).
arguments
 run_root char = ''
 session char = '590_PL2_10w'
 stage char = 'stage2_2d_h5'
end
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
if isempty(run_root), run_root = fullfile(root, 'paper_results', 'b1_adaptq', '20261009T054142Z'); end
f = dir(fullfile(run_root, stage, session, 'bin_*.mat'));
kap = 1e-3 / 6.582119569e-16 / (sqrt(1 - 1 / (1 + 30e3 / 510998.95) ^ 2) * 299792458) * 1e-10;   % qE per meV, 30 kV
rows = {}; acc = []; n = 0;
for k = 1:numel(f)
    L = load(fullfile(f(k).folder, f(k).name), 'R'); R = L.R; m = R.fits([R.fits.main]);
    if abs(R.q) < 0.0054 || abs(R.q) > 0.0156 || ~m.success, continue; end
    c = m.curves; E = c.E(:);
    lf = movmean((c.observed(:) - c.background(:)) ./ c.prefactor(:), max(1, round(100 / median(diff(E)))));
    r = E >= 250 & E <= 1800; Er = E(r); y = lf(r); [ym, i] = max(y);
    lo = find(y(1:i) < ym / 2, 1, 'last'); hi = find(y(i:end) < ym / 2, 1, 'first');
    elo = NaN; ehi = NaN;
    if ~isempty(lo), elo = Er(lo); end
    if ~isempty(hi), ehi = Er(i + hi - 1); end
    rows(end + 1, :) = {R.q, R.N, Er(i), elo, ehi, interp1(Er, y, 300) / ym}; %#ok<AGROW>
    if abs(R.q) >= 0.0084
        if isempty(acc), acc = y / ym; else, acc = acc + y / ym; end
        n = n + 1;
    end
end
T = cell2table(rows, 'VariableNames', {'q', 'N', 'apex_meV', 'low_halfmax_meV', 'high_halfmax_meV', 'LF300_over_apex'});
disp(T); fprintf('high_halfmax NaN: the band stays above half maximum up to 1800 meV\n');
y = acc / n; Z = table();
for Q = [0.0005 0.001 0.002 0.003 0.005 0.0155]
    K = Q ./ (Q ^ 2 + (kap * Er) .^ 2) .^ 2; [~, i] = max(K .* y);
    Z = [Z; table(Q, Er(i), i == 1, 'VariableNames', {'Q', 'max_meV', 'at_window_edge'})]; %#ok<AGROW>
end
fprintf('zero dispersion, mean of %d bins (|q| 0.0085-0.0155), point-Q 2D kinematic factor:\n', n); disp(Z);
out = struct('bins', T, 'zero_dispersion', Z);
end
