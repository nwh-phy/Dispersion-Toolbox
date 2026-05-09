function offsets = qe_adaptive_stack_offsets(traces, min_gap)
%QE_ADAPTIVE_STACK_OFFSETS Compute non-crossing offsets for stacked spectra.
%   offsets = qe_adaptive_stack_offsets(traces, min_gap) returns one offset
%   per trace column. Later traces are shifted upward just enough to stay
%   above the envelope of all earlier shifted traces over the shared energy
%   grid. min_gap is the requested minimum vertical clearance.

if nargin < 2 || isempty(min_gap)
    min_gap = 0;
end

traces = double(traces);
min_gap = max(double(min_gap), 0);
n_traces = size(traces, 2);
offsets = zeros(1, n_traces);

if n_traces <= 1
    return
end

shifted = traces;
for col = 2:n_traces
    previous_envelope = max(shifted(:, 1:col-1), [], 2);
    current_trace = traces(:, col);
    finite_mask = isfinite(previous_envelope) & isfinite(current_trace);

    if any(finite_mask)
        required = max(previous_envelope(finite_mask) + min_gap - current_trace(finite_mask));
    else
        required = offsets(col - 1) + min_gap;
    end

    if ~isfinite(required)
        required = offsets(col - 1) + min_gap;
    end

    offsets(col) = max(required, offsets(col - 1) + min_gap);
    shifted(:, col) = current_trace + offsets(col);
end
end
