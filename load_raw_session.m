function dataset = load_raw_session(file_path, options)
%LOAD_RAW_SESSION  Load raw EELS data (.npy or 4D .mat) and produce a qeData struct.
%
%   dataset = load_raw_session('path/to/file.npy')
%   dataset = load_raw_session('path/to/file.npy', q_crop=[51 408])
%   dataset = load_raw_session('path/to/4d_aligned.mat')
%
%   Supports:
%     - .npy files with companion .json (Nion Swift format)
%     - .mat files containing EELS4D objects or 3D/4D numeric arrays
%
%   Processing pipeline:
%     1. Read data (detect format automatically)
%     2. Sum frames along spatial/time dimensions
%     3. Crop q-dimension (auto-detect or manual)
%     4. ZLP alignment (cross-correlation, ported from Nion Align4DEELS)
%     5. Build calibrated energy axis
%     6. Output standardized qeData struct via make_qe_struct
%
%   Options:
%     q_crop       - [lo, hi] pixel indices for q cropping (default: auto)
%     align_sigma  - Gaussian weight width for iterative alignment (default: 200)
%     align_target - Convergence threshold: var(shifts) (default: 0.5)
%     max_iter     - Maximum alignment iterations (default: 10)
%     dq_Ainv      - Momentum step in 1/Å per pixel (default: auto from path)
%     show_progress - Show alignment convergence figure (default: true)
%
%   See also: read_npy, load_qe_dataset, make_qe_struct

arguments
    file_path {mustBeFile}
    options.q_crop (1,2) double = [NaN NaN]
    options.align_sigma (1,1) double {mustBePositive} = 200
    options.align_target (1,1) double {mustBePositive} = 0.5
    options.max_iter (1,1) double {mustBePositive, mustBeInteger} = 10
    options.dq_Ainv (1,1) double = NaN
    options.show_progress (1,1) logical = true
end

[fdir, fname, ext] = fileparts(file_path);

%% 1. Load data — detect format
energy_meV_from_file = [];  % may be provided by .mat

if strcmpi(ext, '.npy')
    fprintf('Loading NPY file: %s\n', file_path);
    tic;
    raw_data = double(read_npy(file_path));
    fprintf('  Loaded in %.1f s, shape = [%s], %.0f MB\n', ...
        toc, num2str(size(raw_data)), numel(raw_data)*8/1e6);

    % Read companion JSON for energy calibration
    json_path = fullfile(fdir, [fname, '.json']);
    if exist(json_path, 'file') == 2
        json = jsondecode(fileread(json_path));
        fprintf('  JSON metadata loaded.\n');
    else
        json = [];
        fprintf('  Warning: no companion JSON found.\n');
    end
    source_kind = 'npy';

elseif strcmpi(ext, '.mat')
    fprintf('Loading MAT file: %s\n', file_path);
    
    % Auto-add Nion toolbox to path if EELS4D class is not available
    local_ensure_nion_toolbox_on_path();
    
    tic;
    lastwarn('');  % clear last warning
    loaded = load(file_path);
    [warnMsg, ~] = lastwarn();
    fprintf('  Loaded in %.1f s\n', toc);
    json = [];

    % If EELS4D couldn't be instantiated, try reloading with toolbox
    if contains(warnMsg, 'uint32') || contains(warnMsg, '无法实例化')
        fprintf('  Warning: EELS4D class not found, attempting struct extraction...\n');
    end
    
    % Detect data format in .mat
    [raw_data, energy_meV_from_file] = local_extract_mat_data(loaded, file_path);
    source_kind = 'mat4d';
else
    error('load_raw_session:UnsupportedFormat', ...
        'Unsupported file format: %s', ext);
end

%% 2. Frame summation
ndims_data = ndims(raw_data);
fprintf('  Original dimensions: %dD, shape = [%s]\n', ndims_data, num2str(size(raw_data)));

if ndims_data == 4
    summed = squeeze(sum(raw_data, [1 2]));
    fprintf('  Frames summed: %d × %d → 2D [%s]\n', ...
        size(raw_data, 1), size(raw_data, 2), num2str(size(summed)));
elseif ndims_data == 3
    summed = squeeze(sum(raw_data, 1));
    fprintf('  Frames summed: %d frames → 2D [%s]\n', ...
        size(raw_data, 1), num2str(size(summed)));
elseif ndims_data == 2
    summed = raw_data;
    fprintf('  Already 2D, no summation needed.\n');
else
    error('load_raw_session:UnsupportedDims', ...
        'Expected 2D-4D data, got %dD.', ndims_data);
end
clear raw_data

% Ensure summed is [E × q]
if size(summed, 1) < size(summed, 2)
    summed = summed.';
    fprintf('  Transposed to [E × q] = [%d × %d]\n', size(summed, 1), size(summed, 2));
end
n_E = size(summed, 1);
n_q = size(summed, 2);
fprintf('  Data matrix: %d energy × %d momentum channels\n', n_E, n_q);

%% 3. Q cropping
if all(isfinite(options.q_crop))
    q_lo = max(1, round(options.q_crop(1)));
    q_hi = min(n_q, round(options.q_crop(2)));
    fprintf('  Q crop (manual): channels %d : %d\n', q_lo, q_hi);
else
    % Auto-detect usable q range
    % Strategy: find where the intensity profile drops off at the edges.
    % Use a generous threshold (1% of median) to handle pre-aligned data
    % where the intensity is relatively uniform across q.
    q_total = sum(summed, 1, 'omitnan');
    q_total_smooth = smoothdata(q_total, 'gaussian', 20);
    
    % Use median-based threshold (robust to central ZLP peak)
    med_val = median(q_total_smooth(q_total_smooth > 0));
    threshold = med_val * 0.01;  % 1% of median
    valid_q = find(q_total_smooth > threshold);
    
    if isempty(valid_q) || (valid_q(end) - valid_q(1) + 1) < 10
        % Fallback: keep everything
        q_lo = 1; q_hi = n_q;
    else
        q_lo = valid_q(1);
        q_hi = valid_q(end);
    end
    
    % Sanity check: auto-crop should keep at least 50% of channels
    if (q_hi - q_lo + 1) < n_q * 0.5
        fprintf('  Q crop auto too aggressive (%d of %d), keeping all.\n', ...
            q_hi - q_lo + 1, n_q);
        q_lo = 1; q_hi = n_q;
    end
    
    fprintf('  Q crop (auto): channels %d : %d (of %d)\n', q_lo, q_hi, n_q);
end
summed = summed(:, q_lo:q_hi);

%% 4. ZLP alignment
fprintf('  ZLP alignment...\n');
aligned = local_align_zlp(summed, ...
    options.align_sigma, options.align_target, options.max_iter, ...
    options.show_progress);

%% 5. Build energy axis
if ~isempty(energy_meV_from_file)
    % Energy axis from .mat file (EELS4D.ene)
    energy_meV = energy_meV_from_file(:);
    % Re-center on ZLP
    sum_spectrum = sum(aligned, 2, 'omitnan');
    [~, zlp_idx] = max(sum_spectrum);
    energy_meV = energy_meV - energy_meV(zlp_idx);
    fprintf('  Energy axis from file: [%.1f, %.1f] meV, %d points\n', ...
        min(energy_meV), max(energy_meV), numel(energy_meV));
elseif ~isempty(json)
    % Energy axis from JSON
    ax = json.spatial_calibrations;
    e_scale_eV = NaN; e_offset_eV = NaN;
    for ai = 1:numel(ax)
        if contains(ax(ai).units, 'eV', 'IgnoreCase', true)
            e_scale_eV = ax(ai).scale;
            e_offset_eV = ax(ai).offset;
        end
    end
    if ~isfinite(e_scale_eV)
        e_scale_eV = ax(end).scale;
        e_offset_eV = ax(end).offset;
    end
    e_scale_meV = e_scale_eV * 1e3;
    energy_meV = ((1:n_E)' - 1) * e_scale_meV + e_offset_eV * 1e3;
    sum_spectrum = sum(aligned, 2, 'omitnan');
    [~, zlp_idx] = max(sum_spectrum);
    energy_meV = energy_meV - energy_meV(zlp_idx);
    fprintf('  Energy axis from JSON: [%.1f, %.1f] meV, dE = %.1f meV\n', ...
        min(energy_meV), max(energy_meV), e_scale_meV);
else
    % Fallback: guess dE from path or use default
    lower_path = lower(file_path);
    if contains(lower_path, '0.004')
        dE = 4;
    elseif contains(lower_path, '0.002')
        dE = 2;
    else
        dE = 4;  % default
    end
    energy_meV = ((1:n_E)' - 1) * dE;
    sum_spectrum = sum(aligned, 2, 'omitnan');
    [~, zlp_idx] = max(sum_spectrum);
    energy_meV = energy_meV - energy_meV(zlp_idx);
    fprintf('  Energy axis (estimated): [%.1f, %.1f] meV, dE = %.1f meV\n', ...
        min(energy_meV), max(energy_meV), dE);
end

%% 6. Resolve dq
dq_Ainv = options.dq_Ainv;
if ~isfinite(dq_Ainv)
    lower_path = lower(file_path);
    if contains(lower_path, '20w')
        dq_Ainv = 0.0025;
    elseif contains(lower_path, '10w')
        dq_Ainv = 0.005;
    else
        answer = inputdlg({'Enter dq in 1/Å per pixel:'}, ...
            'Momentum Step', [1 40], {'0.005'});
        if isempty(answer)
            error('load_raw_session:DqCancelled', 'dq entry was cancelled.');
        end
        dq_Ainv = str2double(answer{1});
    end
end
fprintf('  dq = %.4f Å⁻¹/pixel\n', dq_Ainv);

%% 7. Build output struct
q_channel = 1:size(aligned, 2);  % re-index from 1 (data is already cropped)
[~, parent_name] = fileparts(fdir);
label = sprintf('%s/%s [%s]', parent_name, fname, source_kind);

dataset = struct();
dataset.source_kind = source_kind;
dataset.source_path = char(file_path);
dataset.label = label;
dataset.dq_Ainv = dq_Ainv;
dataset.qe = make_qe_struct( ...
    aligned, energy_meV, dq_Ainv, q_channel, ...
    source_path=string(file_path), source_kind=string(source_kind), ...
    label=string(label), view_kind="physical", stage_name="raw_import");

fprintf('  Done. qeData struct created: %d E × %d q\n', ...
    size(aligned, 1), size(aligned, 2));

% Save eq3D cache for future fast loading
eq3d_path = fullfile(fdir, 'eq3D_processed.mat');
a3 = aligned; %#ok<NASGU>
e = energy_meV(:).'; %#ok<NASGU>
q = q_channel; %#ok<NASGU>
save(eq3d_path, 'a3', 'e', 'q', '-v7.3');
fprintf('  Saved processed data to: %s\n', eq3d_path);
end


%% ====================================================================
function local_ensure_nion_toolbox_on_path()
%LOCAL_ENSURE_NION_TOOLBOX_ON_PATH  Add Nion 4D toolbox to path if needed.
if exist('EELS4D', 'class') == 8
    return  % already available
end

% Search common locations
desktop = fullfile(getenv('USERPROFILE'), 'Desktop');
candidates = dir(fullfile(desktop, 'Nion-EELS*toolbox'));
for i = 1:numel(candidates)
    toolbox_dir = fullfile(candidates(i).folder, candidates(i).name, '4D-EELS TOOLBOX');
    if isfolder(toolbox_dir)
        addpath(genpath(toolbox_dir));
        fprintf('  Auto-added Nion 4D toolbox: %s\n', toolbox_dir);
        return
    end
end
end


%% ====================================================================
function [raw_data, energy_meV] = local_extract_mat_data(loaded, file_path)
%LOCAL_EXTRACT_MAT_DATA  Extract 3D/4D data and energy axis from a .mat file.
%   Supports:
%     - EELS4D objects (with .data and .ene properties)
%     - Degraded EELS4D (class not on path, loaded as struct or uint32)
%     - Plain numeric arrays (largest 3D/4D variable)

energy_meV = [];
fields = fieldnames(loaded);

% === Case 1: Proper EELS4D object ===
for i = 1:numel(fields)
    val = loaded.(fields{i});
    if isobject(val)
        try
            raw_data = double(val.data);
            if isprop(val, 'ene')
                energy_meV = double(val.ene(:));
            end
            fprintf('  Found EELS4D object in field "%s", data size = [%s]\n', ...
                fields{i}, num2str(size(raw_data)));
            return
        catch
        end
    end
end

% === Case 2: Degraded EELS4D (struct fallback) ===
% When EELS4D class is not on path, MATLAB may load it as a struct
for i = 1:numel(fields)
    val = loaded.(fields{i});
    if isstruct(val) && isfield(val, 'data') && isfield(val, 'ene')
        raw_data = double(val.data);
        energy_meV = double(val.ene(:));
        fprintf('  Found degraded EELS4D struct in field "%s", data size = [%s]\n', ...
            fields{i}, num2str(size(raw_data)));
        return
    end
end

% === Case 3: EELS4D loaded as uint32 (total degradation) ===
% Use whos to peek at the original class, then use matfile for raw access
if nargin >= 2
    try
        info = whos('-file', file_path);
        for i = 1:numel(info)
            if strcmp(info(i).class, 'EELS4D')
                fprintf('  EELS4D detected in file but class unavailable.\n');
                fprintf('  Attempting direct HDF5/v7.3 extraction...\n');
                % v7.3 .mat files are HDF5, try h5read
                try
                    raw_data = double(h5read(file_path, '/data/data'));
                    fprintf('  Extracted data via HDF5: size = [%s]\n', ...
                        num2str(size(raw_data)));
                catch
                end
                try
                    energy_meV = double(h5read(file_path, '/data/ene'));
                    energy_meV = energy_meV(:);
                    fprintf('  Extracted energy via HDF5: %d points\n', ...
                        numel(energy_meV));
                catch
                end
                if exist('raw_data', 'var') && ~isempty(raw_data)
                    return
                end
                % If HDF5 failed, try matfile object access
                try
                    mf = matfile(file_path);
                    raw_data = double(mf.data);
                    fprintf('  Extracted data via matfile: size = [%s]\n', ...
                        num2str(size(raw_data)));
                catch
                end
                break
            end
        end
    catch
    end
end

% === Case 4: Plain numeric array (largest 3D+ field) ===
best_field = '';
best_numel = 0;
for i = 1:numel(fields)
    val = loaded.(fields{i});
    if isnumeric(val) && ndims(val) >= 3 && numel(val) > best_numel
        best_field = fields{i};
        best_numel = numel(val);
    end
end

if ~isempty(best_field)
    raw_data = double(loaded.(best_field));
    fprintf('  Using largest numeric field "%s", size = [%s]\n', ...
        best_field, num2str(size(raw_data)));

    % Try to find energy axis in other fields
    for i = 1:numel(fields)
        if strcmp(fields{i}, best_field); continue; end
        val = loaded.(fields{i});
        if isnumeric(val) && isvector(val) && numel(val) > 100
            energy_meV = double(val(:));
            fprintf('  Found energy axis in field "%s", %d points\n', ...
                fields{i}, numel(energy_meV));
            break
        end
    end
    return
end

error('load_raw_session:NoData', ...
    ['Could not find 3D/4D data in .mat file.\n' ...
     'If this file contains EELS4D objects, ensure the Nion 4D-EELS ' ...
     'toolbox is on your MATLAB path:\n' ...
     '  addpath(genpath(''path/to/Nion-EELS toolbox/4D-EELS TOOLBOX''))']);
end


%% ====================================================================
function aligned = local_align_zlp(data, sigma, target, max_iter, show_progress)
%LOCAL_ALIGN_ZLP  ZLP energy alignment via cross-correlation.
%   Ported from Nion 4D-EELS Toolbox Align4DEELS.m
%   Adapted for 2D (E × q) data.

[n_E, n_q] = size(data);
aligned = data;

% --- Phase 1: Sequential cross-correlation from center outward ---
totint = sum(aligned, 1, 'omitnan');
[~, center] = max(totint);

for qi = center+1 : n_q
    specref = aligned(:, qi-1);
    spec    = aligned(:, qi);
    [cor, lags] = xcorr(spec, specref);
    [~, idx] = max(cor);
    zlpshift = lags(idx);
    aligned(:, qi) = circshift(spec, -zlpshift);
end

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
        subplot(1,3,1);
        imagesc([],[],log(max(aligned,0)+1)); axis xy; colormap hot;
        title('Aligned q-E map'); xlabel('q channel'); ylabel('E pixel');
        subplot(1,3,2);
        histogram(shift_arr,'FaceColor',[0.3 0.5 0.8]);
        title(sprintf('Shift histogram (iter %d)',iter));
        xlabel('Shift (pixels)'); ylabel('Count');
        subplot(1,3,3); hold on;
        scatter(iter,residual,50,'bo','filled');
        set(gca,'YScale','log');
        xlabel('Iteration'); ylabel('Var(shifts)');
        title('Convergence'); xlim([0 iter+1]); box on;
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
