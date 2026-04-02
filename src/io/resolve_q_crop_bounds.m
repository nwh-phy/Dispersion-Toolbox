function [q_lo, q_hi] = resolve_q_crop_bounds(q_crop, n_q, min_q_width, err_id)
%RESOLVE_Q_CROP_BOUNDS  Validate manual q-crop input and return integer bounds.
%
%   [q_lo, q_hi] = resolve_q_crop_bounds(q_crop, n_q, min_q_width, err_id)
%   validates q_crop=[lo,hi] against [1,n_q] and minimum width.

arguments
    q_crop (1,2) double
    n_q (1,1) double {mustBePositive, mustBeInteger}
    min_q_width (1,1) double {mustBePositive, mustBeInteger}
    err_id (1,1) string
end

q_lo_in = q_crop(1);
q_hi_in = q_crop(2);

if ~all(isfinite(q_crop))
    error(err_id, ...
        'Invalid q_crop = [%g, %g]: values must be finite.', q_lo_in, q_hi_in);
end

if any(abs(q_crop - round(q_crop)) > eps(max(abs(q_crop), 1)))
    error(err_id, ...
        ['Invalid q_crop = [%g, %g]: bounds must be integer channel indices. ' ...
         'Valid range is [1, %d].'], ...
        q_lo_in, q_hi_in, n_q);
end

q_lo = round(q_lo_in);
q_hi = round(q_hi_in);

if q_lo < 1 || q_hi > n_q
    error(err_id, ...
        ['Invalid q_crop = [%d, %d]: bounds are out of range. ' ...
         'Valid channel range is [1, %d].'], ...
        q_lo, q_hi, n_q);
end

if q_lo > q_hi
    error(err_id, ...
        ['Invalid q_crop = [%d, %d]: lo must be <= hi. ' ...
         'Valid channel range is [1, %d].'], ...
        q_lo, q_hi, n_q);
end

q_width = q_hi - q_lo + 1;
if q_width < min_q_width
    error(err_id, ...
        ['Invalid q_crop = [%d, %d]: width is %d channels. ' ...
         'Minimum required width is %d channels.'], ...
        q_lo, q_hi, q_width, min_q_width);
end
end
