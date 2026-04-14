function peak_height = measure_peak_height(energy_meV, spectrum, peak_energy_meV, gamma_meV, opts)
%MEASURE_PEAK_HEIGHT  Estimate a measured peak height from a local spectrum window.
%
%   peak_height = measure_peak_height(E, S, E0, Gamma)
%   returns the maximum finite spectral intensity in a window centered on
%   E0 with half-width set by the fitted linewidth Gamma.
%
%   This helper is intentionally simple: it measures the local peak height
%   from the spectrum itself, rather than from the fitted model curve. That
%   keeps downstream amplitude analyses tied to the measured EELS signal.

arguments
    energy_meV (:,1) double
    spectrum   (:,1) double
    peak_energy_meV (1,1) double
    gamma_meV (1,1) double
    opts.min_half_window_pts (1,1) double {mustBePositive, mustBeInteger} = 3
end

E = double(energy_meV(:));
S = double(spectrum(:));

if isempty(E) || isempty(S) || numel(E) ~= numel(S)
    peak_height = NaN;
    return
end

[~, nearest_idx] = min(abs(E - peak_energy_meV));

if numel(E) > 1
    dE = median(abs(diff(E)), 'omitnan');
else
    dE = NaN;
end

if ~isfinite(dE) || dE <= 0
    half_window_pts = opts.min_half_window_pts;
else
    half_window_pts = max(opts.min_half_window_pts, round(abs(gamma_meV) / dE));
end

win = max(1, nearest_idx - half_window_pts):min(numel(S), nearest_idx + half_window_pts);
window_vals = S(win);
window_vals = window_vals(isfinite(window_vals));

if isempty(window_vals)
    peak_height = NaN;
else
    peak_height = max(window_vals);
end
end
