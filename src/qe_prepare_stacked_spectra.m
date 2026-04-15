function stack = qe_prepare_stacked_spectra(energy_meV, spectra, q_axis, opts)
%QE_PREPARE_STACKED_SPECTRA Prepare stacked single-spectrum traces for plotting.
%   stack = qe_prepare_stacked_spectra(energy_meV, spectra, q_axis, opts)
%   selects representative q-channels over the requested range and applies
%   a fixed vertical offset to each trace.
%
%   Required opts fields:
%     q_start   start q value for target sampling
%     q_end     end q value for target sampling
%     q_step    target sampling spacing
%     offset    vertical offset between neighboring traces

if nargin < 4 || ~isstruct(opts)
    error('qe_prepare_stacked_spectra:InvalidOptions', ...
        'opts must be a struct with q_start, q_end, q_step, and offset fields.');
end
if ~isfield(opts, 'q_start') || ~isfield(opts, 'q_end') || ~isfield(opts, 'q_step')
    error('qe_prepare_stacked_spectra:MissingOptions', ...
        'opts.q_start, opts.q_end, and opts.q_step are required.');
end
if ~isfield(opts, 'offset')
    opts.offset = 0;
end

energy_meV = double(energy_meV(:));
spectra = double(spectra);
q_axis = double(q_axis(:)).';

if size(spectra, 1) ~= numel(energy_meV)
    error('qe_prepare_stacked_spectra:SizeMismatch', ...
        'spectra row count must match energy axis length.');
end
if size(spectra, 2) ~= numel(q_axis)
    error('qe_prepare_stacked_spectra:SizeMismatch', ...
        'spectra column count must match q axis length.');
end

q_axis = double(q_axis(:)).';
q_min = min(opts.q_start, opts.q_end);
q_max = max(opts.q_start, opts.q_end);
q_step = max(double(opts.q_step), eps);

targets = q_min:q_step:q_max;
if isempty(targets)
    targets = q_min;
end

q_indices = zeros(1, numel(targets));
for idx = 1:numel(targets)
    [~, q_indices(idx)] = min(abs(q_axis - targets(idx)));
end
q_indices = unique(q_indices, 'stable');

if isempty(q_indices)
    [~, nearest_idx] = min(abs(q_axis - mean([q_min, q_max])));
    q_indices = nearest_idx;
end

traces = double(spectra(:, q_indices));
offsets = (0:(numel(q_indices) - 1)) * double(opts.offset);
shifted_traces = traces + reshape(offsets, 1, []);

stack = struct();
stack.energy_meV = double(energy_meV(:));
stack.q_indices = q_indices;
stack.q_values = q_axis(q_indices);
stack.traces = traces;
stack.offsets = offsets;
stack.shifted_traces = shifted_traces;
end
