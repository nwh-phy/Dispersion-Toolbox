function offset = qe_auto_stack_offset(traces, scale_factor)
%QE_AUTO_STACK_OFFSET Estimate a readable vertical spacing for stacked spectra.

if nargin < 2 || isempty(scale_factor)
    scale_factor = 1.2;
end

traces = double(traces);
scale_factor = max(double(scale_factor), eps);
spans = [];

for col = 1:size(traces, 2)
    y = traces(:, col);
    y = y(isfinite(y));
    if numel(y) < 2
        continue
    end

    lo = prctile(y, 5);
    hi = prctile(y, 95);
    span = hi - lo;
    if ~isfinite(span) || span <= eps
        span = max(y) - min(y);
    end
    if isfinite(span) && span > eps
        spans(end + 1) = span; %#ok<AGROW>
    end
end

if isempty(spans)
    finite_values = traces(isfinite(traces));
    if isempty(finite_values)
        offset = 1;
        return
    end

    span = max(finite_values) - min(finite_values);
    if isfinite(span) && span > eps
        offset = span * scale_factor;
    else
        offset = 1;
    end
    return
end

offset = median(spans) * scale_factor;
if ~isfinite(offset) || offset <= eps
    offset = 1;
end
end
