function result = measure_peak_fwhm(energy_meV, peak_curve)
%MEASURE_PEAK_FWHM Numerically measure peak FWHM from a fitted peak curve.
%   The half-maximum is measured relative to the local curve floor. This is
%   more stable for asymmetric Fano-like curves than assuming a zero baseline.

arguments
    energy_meV (:,1) double
    peak_curve (:,1) double
end

E = double(energy_meV(:));
Y = double(peak_curve(:));
valid = isfinite(E) & isfinite(Y);
E = E(valid);
Y = Y(valid);

result = struct( ...
    'fwhm_meV', NaN, ...
    'left_meV', NaN, ...
    'right_meV', NaN, ...
    'apex_meV', NaN, ...
    'half_level', NaN, ...
    'baseline_level', NaN, ...
    'peak_level', NaN, ...
    'status', "invalid");

if numel(E) < 5
    return
end

[E, order] = sort(E);
Y = Y(order);

[peak_level, apex_idx] = max(Y);
baseline_level = min(Y);
if baseline_level >= 0 && peak_level > 0 && baseline_level / peak_level < 0.10
    baseline_level = 0;
end
contrast = peak_level - baseline_level;
if ~isfinite(contrast) || contrast <= 0
    result.status = "flat";
    return
end

half_level = baseline_level + 0.5 * contrast;
left = local_crossing(E, Y, half_level, apex_idx, -1);
right = local_crossing(E, Y, half_level, apex_idx, 1);

result.apex_meV = E(apex_idx);
result.half_level = half_level;
result.baseline_level = baseline_level;
result.peak_level = peak_level;

if ~isfinite(left) || ~isfinite(right) || right <= left
    result.status = "no_two_sided_crossing";
    return
end

result.left_meV = left;
result.right_meV = right;
result.fwhm_meV = right - left;
result.status = "ok";
end


function x_cross = local_crossing(E, Y, level, apex_idx, direction)
x_cross = NaN;
if direction < 0
    idx_range = apex_idx:-1:2;
    for k = idx_range
        y1 = Y(k);
        y0 = Y(k - 1);
        if (y1 - level) * (y0 - level) <= 0 && y1 ~= y0
            x_cross = interp1([y0 y1], [E(k - 1) E(k)], level, 'linear');
            return
        end
    end
else
    idx_range = apex_idx:(numel(E) - 1);
    for k = idx_range
        y0 = Y(k);
        y1 = Y(k + 1);
        if (y0 - level) * (y1 - level) <= 0 && y1 ~= y0
            x_cross = interp1([y0 y1], [E(k) E(k + 1)], level, 'linear');
            return
        end
    end
end
end
