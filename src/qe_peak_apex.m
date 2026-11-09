function [e, info] = qe_peak_apex(E, y, window, opts)
%QE_PEAK_APEX  Energy of the most prominent local maximum of a spectrum inside a window.
%   e = qe_peak_apex(E, y, [E1 E2]) finds the local maxima of a moving mean of y
%   (opts.smooth channels, default 11 = 44 meV at 4 meV per channel) strictly
%   inside [E1, E2] with prominence >= opts.min_prominence * max, takes the most
%   prominent one and refines it with a parabola fitted to the unsmoothed y
%   within +/- opts.half_width meV. e is NaN when there is no interior maximum
%   or the parabola opens upwards; info.reason says why.
%   This is a model-free peak position (no background, no kinematic factor);
%   for a peak on a ZLP tail subtract qe_powerlaw_tail first.
%
%   See also: qe_powerlaw_tail, qe_zlp_center, fit_loss_function
arguments
    E (:,1) double
    y (:,1) double
    window (1,2) double
    opts.smooth (1,1) double {mustBePositive, mustBeInteger} = 11
    opts.min_prominence (1,1) double = 0.005
    opts.half_width (1,1) double {mustBePositive} = 80
end
e = NaN; info = struct('reason', '', 'prominence', NaN, 'coarse', NaN);
ys = movmean(y, opts.smooth); w = find(E >= window(1) & E <= window(2));
[tf, pr] = islocalmax(ys(w), 'MinProminence', opts.min_prominence * max(ys(w)));
tf([1 end]) = false;
if ~any(tf), info.reason = 'no interior maximum'; return; end
c = find(tf); [~, m] = max(pr(c)); i0 = w(c(m));
info.prominence = pr(c(m)); info.coarse = E(i0);
id = find(abs(E - E(i0)) <= opts.half_width); p = polyfit(E(id) - E(i0), y(id), 2);
if p(1) >= 0, info.reason = 'parabola not concave'; return; end
e = E(i0) - p(2) / (2 * p(1));
if abs(e - E(i0)) > opts.half_width, e = NaN; info.reason = 'vertex outside the fit range'; end
end
