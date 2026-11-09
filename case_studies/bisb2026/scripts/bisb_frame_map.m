function S = bisb_frame_map(id, opts)
%BISB_FRAME_MAP  Frame-aligned, drift-following q-E map of a registered dataset (cached).
%   S = bisb_frame_map(id) takes dataset id from bisb_datasets, aligns every
%   frame with qe_align_sequence_frames (the raw-frame preparation of
%   b1_low_branch_matlab: one integer ZLP shift per frame, beam q-drift from the
%   elastic centroid), bins the native channels |q| <= q_max with
%   qe_transport_frames following the drift, and sums the frames, also per
%   contiguous time block. Raw counts: no normalisation, background or
%   kinematic factor.
%
%   The result is cached in paper_results/frame_maps/<id>_<hash8>.mat, the hash
%   covering the registry entry and the options; refresh = true rebuilds it.
%
%   Fields of S: dataset (registry entry), opts, E (meV, relative to the
%   reference ZLP pixel), q (1/A, channel centres), M (E x q summed counts,
%   double), Mb (E x q x block), block_frames, offsets_meV (energy shift removed
%   from each frame), drift_Ainv (per frame), sigma_probe, beam, file.
%
%   See also: bisb_datasets, bisb_quicklook, qe_align_sequence_frames, qe_transport_frames
arguments
    id (1,:) char
    opts.q_max (1,1) double {mustBePositive} = 0.03
    opts.n_blocks (1,1) double {mustBePositive, mustBeInteger} = 6
    opts.refresh (1,1) logical = false
end
d = bisb_datasets(id);
root = bisb_find_project_root(fileparts(mfilename('fullpath')));
o = rmfield(opts, 'refresh');
key = b1_sha256(struct('dataset', rmfield(d, 'path'), 'opts', o, 'version', 'bisb_frame_map v1'));
folder = fullfile(root, 'paper_results', 'frame_maps');
file = fullfile(folder, sprintf('%s_%s.mat', d.id, key(1:8)));
if isfile(file) && ~opts.refresh
    L = load(file, 'S'); S = L.S; return
end
if ~isfolder(folder), mkdir(folder); end
assert(strcmp(b1_sha256({'file', d.path}), d.sha256), 'bisb_frame_map:Source', ...
    'Raw file of %s differs from the hash in bisb_datasets', d.id);
A = qe_align_sequence_frames(d.path, d.zero, d.dq, q_select=opts.q_max + 0.002);
keep = abs(A.Q) <= opts.q_max + 1e-12; centres = A.Q(keep);
frames = qe_transport_frames(A.F, A.Q, d.dq, centres, d.dq * ones(size(centres)), A.drift_Ainv);
T = size(frames, 2); blk = ceil((1:T) * opts.n_blocks / T);
M = reshape(sum(double(frames), 2), size(frames, 1), []);
Mb = zeros(size(M, 1), size(M, 2), opts.n_blocks);
for b = 1:opts.n_blocks
    Mb(:, :, b) = reshape(sum(double(frames(:, blk == b, :)), 2), size(frames, 1), []);
end
dE = median(diff(A.E));
S = struct('dataset', d, 'opts', o, 'E', A.E, 'q', centres, 'M', M, 'Mb', Mb, ...
    'block_frames', accumarray(blk(:), 1).', 'offsets_meV', A.alignment.measured_offset_pixels * dE, ...
    'drift_Ainv', A.drift_Ainv, 'sigma_probe', A.sigma_probe, 'beam', A.beam, 'file', file);
save(file, 'S', '-v7.3');
fprintf('bisb_frame_map %s: %d frames, energy shifts %.0f..%+.0f meV, drift %.2f..%+.2f channels -> %s\n', ...
    d.id, T, min(S.offsets_meV), max(S.offsets_meV), min(A.drift_Ainv)/d.dq, max(A.drift_Ainv)/d.dq, file);
end
