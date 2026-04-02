function dataset = load_npy_session(npy_path, options)
%LOAD_NPY_SESSION  Load raw Nion .npy + .json and produce a qeData struct.
%
%   dataset = load_npy_session('path/to/file.npy')
%   dataset = load_npy_session('path/to/file.npy', q_crop=[51 408])
%
%   Processing pipeline:
%     1. Read .npy (4D float32) + companion .json (metadata)
%     2. Sum frames along x-dimension (300 frames → 1)
%     3. Crop q-dimension (auto-detect or manual)
%     4. ZLP alignment (cross-correlation, ported from Nion Align4DEELS)
%     5. Build calibrated energy axis from JSON
%     6. Output standardized qeData struct via make_qe_struct
%
%   Options:
%     q_crop       - [lo, hi] pixel indices for q cropping (default: auto)
%                    Example: [51 408]
%     align_sigma  - Gaussian weight width for iterative alignment (default: 200)
%     align_target - Convergence threshold: var(shifts) (default: 0.5)
%     max_iter     - Maximum alignment iterations (default: 10)
%     dq_Ainv      - Momentum step in 1/Å per pixel (default: auto from path)
%     show_progress - Show alignment convergence figure (default: true)
%
%   See also: read_npy, load_qe_dataset, make_qe_struct

arguments
    npy_path {mustBeFile}
    options.q_crop (1,2) double = [NaN NaN]
    options.align_sigma (1,1) double {mustBePositive} = 200
    options.align_target (1,1) double {mustBePositive} = 0.5
    options.max_iter (1,1) double {mustBePositive, mustBeInteger} = 10
    options.dq_Ainv (1,1) double = NaN
    options.show_progress (1,1) logical = true
end

%% 1. Read .npy + .json
fprintf('Loading NPY file: %s\n', npy_path);
tic;
raw_data = double(read_npy(npy_path));
fprintf('  Loaded in %.1f s, shape = [%s], %.0f MB\n', ...
    toc, num2str(size(raw_data)), numel(raw_data)*8/1e6);

% Read companion JSON
[fdir, fname, ~] = fileparts(npy_path);
json_path = fullfile(fdir, [fname, '.json']);
if exist(json_path, 'file') ~= 2
    error('load_npy_session:NoJSON', ...
        'Companion JSON not found: %s', json_path);
end
json = jsondecode(fileread(json_path));
fprintf('  JSON metadata loaded: %s\n', json_path);

%% 2. Frame summation
ndims_data = ndims(raw_data);
fprintf('  Original dimensions: %dD\n', ndims_data);

if ndims_data == 4
    % Shape: y × x × E × q  (Nion convention after readNPYtoEELS permute)
    % But raw .npy from Nion Swift is: y × x × E × q (already)
    % Sum over first two dimensions (y, x = frames)
    summed = squeeze(sum(raw_data, [1 2]));
    fprintf('  Frames summed: %d × %d frames → 2D [%s]\n', ...
        size(raw_data, 1), size(raw_data, 2), num2str(size(summed)));
elseif ndims_data == 3
    % Shape: x × E × q  (sequence without y)
    summed = squeeze(sum(raw_data, 1));
    fprintf('  Frames summed: %d frames → 2D [%s]\n', ...
        size(raw_data, 1), num2str(size(summed)));
elseif ndims_data == 2
    summed = raw_data;
    fprintf('  Already 2D, no summation needed.\n');
else
    error('load_npy_session:UnsupportedDims', ...
        'Expected 2D-4D data, got %dD.', ndims_data);
end
clear raw_data  % free memory

% Ensure summed is [E × q]
% From JSON: spatial_calibrations order matches dimensions.
% Last calibration with units 'eV' is the energy axis.
% Heuristic: energy dim has more pixels than q dim for typical Nion data.
if size(summed, 1) < size(summed, 2)
    summed = summed.';
    fprintf('  Transposed to [E × q] = [%d × %d]\n', size(summed, 1), size(summed, 2));
end
n_E = size(summed, 1);
n_q = size(summed, 2);
fprintf('  Data matrix: %d energy × %d momentum channels\n', n_E, n_q);

%% 3. Q cropping
if all(isfinite(options.q_crop))
    % Manual crop
    min_q_width = 10;
    [q_lo, q_hi] = resolve_q_crop_bounds( ...
        options.q_crop, n_q, min_q_width, "load_npy_session:InvalidQCrop");
    fprintf('  Q crop (manual): channels %d : %d\n', q_lo, q_hi);
else
    % Auto-detect: find q range where total intensity > 5% of peak
    q_total = sum(summed, 1, 'omitnan');
    q_total_smooth = smoothdata(q_total, 'gaussian', 20);
    threshold = max(q_total_smooth) * 0.05;
    valid_q = find(q_total_smooth > threshold);
    if isempty(valid_q)
        q_lo = 1; q_hi = n_q;
    else
        q_lo = valid_q(1);
        q_hi = valid_q(end);
    end
    fprintf('  Q crop (auto): channels %d : %d (of %d)\n', q_lo, q_hi, n_q);
end
summed = summed(:, q_lo:q_hi);
n_q_cropped = size(summed, 2);

%% 4. ZLP alignment (ported from Nion Align4DEELS)
fprintf('  ZLP alignment...\n');
aligned = local_align_zlp(summed, ...
    options.align_sigma, options.align_target, options.max_iter, ...
    options.show_progress);

%% 5. Build energy axis from JSON
ax = json.spatial_calibrations;

% Find energy calibration — it's the last axis with units containing 'eV'
e_scale_eV = NaN;
e_offset_eV = NaN;
for ai = 1:numel(ax)
    if contains(ax(ai).units, 'eV', 'IgnoreCase', true)
        e_scale_eV = ax(ai).scale;
        e_offset_eV = ax(ai).offset;
    end
end

if ~isfinite(e_scale_eV)
    % Fallback: use the last calibration axis
    e_scale_eV = ax(end).scale;
    e_offset_eV = ax(end).offset;
    fprintf('  Warning: no eV axis found in JSON, using last axis.\n');
end

% Build energy axis in meV
e_scale_meV = e_scale_eV * 1e3;
energy_meV = ((1:n_E) - 1) * e_scale_meV + e_offset_eV * 1e3;

% Center on ZLP maximum
sum_spectrum = sum(aligned, 2, 'omitnan');
[~, zlp_idx] = max(sum_spectrum);
energy_meV = energy_meV - energy_meV(zlp_idx);
energy_meV = energy_meV(:);

fprintf('  Energy axis: [%.1f, %.1f] meV, dE = %.1f meV, ZLP at pixel %d\n', ...
    min(energy_meV), max(energy_meV), e_scale_meV, zlp_idx);

%% 6. Resolve dq
dq_Ainv = options.dq_Ainv;
if ~isfinite(dq_Ainv)
    % Try to infer from path
    lower_path = lower(npy_path);
    if contains(lower_path, '20w')
        dq_Ainv = 0.0025;
    elseif contains(lower_path, '10w')
        dq_Ainv = 0.005;
    else
        % Prompt user
        answer = inputdlg({'Enter dq in 1/Å per pixel:'}, ...
            'Momentum Step', [1 40], {'0.005'});
        if isempty(answer)
            error('load_npy_session:DqCancelled', 'dq entry was cancelled.');
        end
        dq_Ainv = str2double(answer{1});
    end
end
fprintf('  dq = %.4f Å⁻¹/pixel\n', dq_Ainv);

%% 7. Build output struct
q_channel = 1:size(aligned, 2);  % re-index from 1 (data is already cropped)

% Get label from path
[~, parent_name] = fileparts(fdir);
label = sprintf('%s/%s [npy]', parent_name, fname);

dataset = struct();
dataset.source_kind = 'npy';
dataset.source_path = char(npy_path);
dataset.label = label;
dataset.dq_Ainv = dq_Ainv;
dataset.qe = make_qe_struct( ...
    aligned, energy_meV, dq_Ainv, q_channel, ...
    source_path=string(npy_path), source_kind="npy", ...
    label=string(label), view_kind="physical", stage_name="npy_import");

fprintf('  Done. qeData struct created: %d E × %d q\n', ...
    size(aligned, 1), size(aligned, 2));

% Optionally save eq3D.mat for future fast loading
eq3d_path = fullfile(fdir, 'eq3D_from_npy.mat');
a3 = aligned;
e = energy_meV(:).';
q = q_channel;
save(eq3d_path, 'a3', 'e', 'q', '-v7.3');
fprintf('  Saved processed data to: %s\n', eq3d_path);
end


%% ====================================================================
function aligned = local_align_zlp(data, sigma, target, max_iter, show_progress)
%LOCAL_ALIGN_ZLP  ZLP energy alignment via cross-correlation.
%   Ported from Nion 4D-EELS Toolbox Align4DEELS.m
%   Adapted for 2D (E × q) data.

[n_E, n_q] = size(data);
aligned = data;

% --- Phase 1: Sequential cross-correlation from center outward ---
% Find center channel (brightest ZLP)
totint = sum(aligned, 1, 'omitnan');
[~, center] = max(totint);

% Align rightward from center
for qi = center+1 : n_q
    specref = aligned(:, qi-1);
    spec    = aligned(:, qi);
    [cor, lags] = xcorr(spec, specref);
    [~, idx] = max(cor);
    zlpshift = lags(idx);
    aligned(:, qi) = circshift(spec, -zlpshift);
end

% Align leftward from center
for qi = center-1 : -1 : 1
    specref = aligned(:, qi+1);
    spec    = aligned(:, qi);
    [cor, lags] = xcorr(spec, specref);
    [~, idx] = max(cor);
    zlpshift = lags(idx);
    aligned(:, qi) = circshift(spec, -zlpshift);
end

fprintf('    Phase 1: Sequential cross-correlation done (center = %d)\n', center);

% --- Phase 2: Iterative global alignment ---
if show_progress
    fig_h = figure('Name', 'ZLP Alignment', 'NumberTitle', 'off', ...
        'Position', [100 100 900 300]);
end

gaussian_wt = exp(-((-n_E+1:n_E-1).^2) / sigma^2);

residual = inf;
iter = 1;

while residual > target && iter <= max_iter
    % Build reference: sum spectrum excluding dim q channels
    diffref = sum(aligned, 1, 'omitnan');
    bright_threshold = 0.8 * max(diffref);
    bright_mask = diffref >= bright_threshold;
    specref = sum(aligned(:, bright_mask), 2, 'omitnan');

    shift_arr = zeros(1, n_q);

    for qi = 1:n_q
        spec = aligned(:, qi);
        cor = xcorr(spec, specref);
        weighted_cor = cor(:).' .* gaussian_wt;
        [~, idx] = max(weighted_cor);
        zlpshift = n_E - idx;
        shift_arr(qi) = zlpshift;
        aligned(:, qi) = circshift(spec, zlpshift);
    end

    residual = var(shift_arr);

    if show_progress && ishandle(fig_h)
        subplot(1, 3, 1);
        imagesc([], [], log(max(aligned, 0) + 1));
        axis xy; colormap hot;
        title('Aligned q-E map'); xlabel('q channel'); ylabel('E pixel');

        subplot(1, 3, 2);
        histogram(shift_arr, 'FaceColor', [0.3 0.5 0.8]);
        title(sprintf('Shift histogram (iter %d)', iter));
        xlabel('Shift (pixels)'); ylabel('Count');

        subplot(1, 3, 3); hold on;
        scatter(iter, residual, 50, 'bo', 'filled');
        set(gca, 'YScale', 'log');
        xlabel('Iteration'); ylabel('Var(shifts)');
        title('Convergence');
        xlim([0 iter + 1]);
        box on;
        drawnow;
    end

    fprintf('    Iteration #%d: residual = %.2g\n', iter, residual);
    iter = iter + 1;
end

if residual < target
    fprintf('    Phase 2: Converged (residual %.2g < target %.2g)\n', residual, target);
else
    fprintf('    Phase 2: Max iterations reached (residual %.2g)\n', residual);
end
end
