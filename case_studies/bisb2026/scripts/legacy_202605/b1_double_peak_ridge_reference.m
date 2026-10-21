function ref = b1_double_peak_ridge_reference(energy_axis, q_axis, intensity, options)
%B1_DOUBLE_PEAK_RIDGE_REFERENCE Build lower/upper B1 ridge references.
%
% The ridge positions are used only as fit seeds or window centers. They are
% not final extracted points.

arguments
    energy_axis (:,1) double
    q_axis (:,1) double
    intensity double
    options.energyWindowMeV (1,2) double = [300 2100]
    options.lowerBandMeV (1,2) double = [300 1050]
    options.upperBandMeV (1,2) double = [850 2100]
    options.smoothWindow (1,1) double = 9
    options.minSeparationMeV (1,1) double = 10
end

energy_axis = double(energy_axis(:));
q_axis = double(q_axis(:));
intensity = double(intensity);
if size(intensity, 1) ~= numel(energy_axis) && ...
        size(intensity, 2) == numel(energy_axis)
    intensity = intensity.';
end
if size(intensity, 1) ~= numel(energy_axis) || ...
        size(intensity, 2) ~= numel(q_axis)
    error('b1_double_peak_ridge_reference:InvalidSize', ...
        'Intensity must be energy-by-q.');
end

energy_window = sort(options.energyWindowMeV);
lower_band = local_clamp_band(sort(options.lowerBandMeV), energy_window);
upper_band = local_clamp_band(sort(options.upperBandMeV), energy_window);
smooth_window = local_odd_window(options.smoothWindow);

n = numel(q_axis);
lower_energy = NaN(n, 1);
upper_energy = NaN(n, 1);
lower_score = NaN(n, 1);
upper_score = NaN(n, 1);

for i = 1:n
    y = double(intensity(:, i));
    y = local_smooth_trace(y, smooth_window);
    y = y - median(y(isfinite(y)), 'omitnan');
    [lower_energy(i), lower_score(i)] = local_band_peak( ...
        energy_axis, y, lower_band);
    [upper_energy(i), upper_score(i)] = local_band_peak( ...
        energy_axis, y, upper_band);
    if ~isfinite(lower_energy(i)) || ~isfinite(upper_energy(i)) || ...
            upper_energy(i) - lower_energy(i) < options.minSeparationMeV
        [pair, score_pair] = local_two_peak_pair(energy_axis, y, ...
            energy_window, options.minSeparationMeV);
        if numel(pair) == 2
            lower_energy(i) = pair(1);
            upper_energy(i) = pair(2);
            lower_score(i) = score_pair(1);
            upper_score(i) = score_pair(2);
        end
    end
end

lower_energy = local_fill_missing_ridge(q_axis, lower_energy);
upper_energy = local_fill_missing_ridge(q_axis, upper_energy);
swap = isfinite(lower_energy) & isfinite(upper_energy) & ...
    lower_energy > upper_energy;
tmp = lower_energy(swap);
lower_energy(swap) = upper_energy(swap);
upper_energy(swap) = tmp;

lower_points = local_ridge_table(q_axis, lower_energy, lower_score, ...
    1, 'b1_double_peak_lower');
upper_points = local_ridge_table(q_axis, upper_energy, upper_score, ...
    2, 'b1_double_peak_upper');

ref = struct();
ref.lower_points = lower_points;
ref.upper_points = upper_points;
ref.reference_points = [lower_points; upper_points];
end


function band = local_clamp_band(band, window)
band(1) = max(band(1), window(1));
band(2) = min(band(2), window(2));
if band(2) <= band(1)
    band = window;
end
end


function y = local_smooth_trace(y, window)
if window < 3 || numel(y) < window
    return
end
try
    y = smoothdata(y, 'sgolay', window);
catch
    y = smoothdata(y, 'movmedian', window);
end
end


function [energy, score] = local_band_peak(energy_axis, y, band)
energy = NaN;
score = NaN;
mask = energy_axis >= band(1) & energy_axis <= band(2) & isfinite(y);
if ~any(mask)
    return
end
yy = y(mask);
ee = energy_axis(mask);
[score, idx] = max(yy, [], 'omitnan');
if isempty(idx) || ~isfinite(score)
    return
end
energy = ee(idx);
end


function [pair, score] = local_two_peak_pair(energy_axis, y, window, min_sep)
pair = [];
score = [];
mask = energy_axis >= window(1) & energy_axis <= window(2) & isfinite(y);
if nnz(mask) < 3
    return
end
ee = energy_axis(mask);
yy = max(y(mask), 0);
try
    [pks, locs] = findpeaks(yy, ee, ...
        'MinPeakDistance', min_sep, 'SortStr', 'descend');
catch
    pks = [];
    locs = [];
end
if numel(locs) < 2
    [sorted_y, order] = sort(yy, 'descend');
    locs = ee(order);
    pks = sorted_y;
end
if numel(locs) < 2
    return
end
for i = 1:numel(locs)
    for j = (i + 1):numel(locs)
        if abs(locs(i) - locs(j)) >= min_sep
            pair = sort([locs(i); locs(j)]);
            score = [pks(i); pks(j)];
            [~, order] = sort([locs(i); locs(j)]);
            score = score(order);
            return
        end
    end
end
end


function energy = local_fill_missing_ridge(q_axis, energy)
valid = isfinite(q_axis) & isfinite(energy);
if all(valid)
    return
end
if nnz(valid) >= 2
    energy(~valid) = interp1(q_axis(valid), energy(valid), ...
        q_axis(~valid), 'pchip', 'extrap');
elseif nnz(valid) == 1
    energy(~valid) = energy(valid);
end
end


function tbl = local_ridge_table(q_axis, energy, score, branch, label)
n = numel(q_axis);
tbl = table(q_axis(:), abs(q_axis(:)), energy(:), ...
    repmat(branch, n, 1), repmat({label}, n, 1), score(:), ...
    repmat({'ridge_reference'}, n, 1), ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'energy_meV', ...
    'branch', 'branch_label', 'ridge_score', 'ridge_source'});
end


function window = local_odd_window(value)
window = max(3, round(double(value)));
if mod(window, 2) == 0
    window = window + 1;
end
end
