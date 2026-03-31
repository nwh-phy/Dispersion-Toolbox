function dataset = load_qe_dataset(path, dqOverride, options)
%LOAD_QE_DATASET  Load q-E data from eq3D.mat or raw .npy + .json.
%
%   dataset = load_qe_dataset('path/to/folder')
%   dataset = load_qe_dataset('path/to/eq3D.mat')
%   dataset = load_qe_dataset('path/to/file.npy')
%   dataset = load_qe_dataset('path/to/file.npy', NaN, q_crop=[51 408])
%
%   Supports:
%     - eq3D.mat files (preprocessed by dealTBG or similar)
%     - Raw .npy files with companion .json (Nion Swift format)
%
%   For .npy files, the full preprocessing pipeline runs:
%   read → frame sum → q crop → ZLP alignment → energy calibration
%
%   See also: load_npy_session, make_qe_struct, read_npy
arguments
    path {mustBeTextScalar}
    dqOverride (1,1) double = NaN
    options.q_crop (1,2) double = [NaN NaN]
end

[source_path, file_type] = local_resolve_source(string(path));

if strcmp(file_type, 'npy') || strcmp(file_type, 'mat_raw')
    % Check for cached processed result first
    [raw_dir, ~, ~] = fileparts(char(source_path));
    cache_path = fullfile(raw_dir, 'eq3D_processed.mat');
    if exist(cache_path, 'file') == 2
        cache_info = dir(cache_path);
        raw_info = dir(char(source_path));
        if cache_info.datenum >= raw_info.datenum
            fprintf('  Using cached eq3D_processed.mat (newer than raw)\n');
            dataset = load_qe_dataset(cache_path, dqOverride);
            return
        end
    end
    % No cache or stale — process raw data
    dataset = load_raw_session(char(source_path), ...
        q_crop=options.q_crop, dq_Ainv=dqOverride);
    return
end

% eq3D.mat loading path (file_type == 'mat_eq3d')
dq_Ainv = local_resolve_dq(source_path, dqOverride);

loaded = load(char(source_path));
if ~isfield(loaded, "a3") || ~isfield(loaded, "e")
    error("load_qe_dataset:InvalidEq3D", ...
        "File %s must contain variables 'a3' and 'e'.", source_path);
end

intensity = double(loaded.a3);
energy_meV = double(loaded.e(:));
intensity = local_match_energy_dimension(intensity, energy_meV, source_path);

if isfield(loaded, "q") && ~isempty(loaded.q) && numel(loaded.q) == size(intensity, 2)
    q_channel = double(loaded.q(:)).';
else
    q_channel = 1:size(intensity, 2);
end

label = local_make_label(source_path);

dataset = struct();
dataset.source_kind = 'eq3d';
dataset.source_path = char(source_path);
dataset.label = char(label);
dataset.dq_Ainv = dq_Ainv;
dataset.qe = make_qe_struct( ...
    intensity, energy_meV, dq_Ainv, q_channel, ...
    source_path=source_path, source_kind="eq3d", ...
    label=label, view_kind="physical", stage_name="eq3d");
end


function [source_path, file_type] = local_resolve_source(path)
%LOCAL_RESOLVE_SOURCE  Find eq3D.mat or .npy in the given path.
%   Returns the resolved path and type ('mat' or 'npy').
candidate = char(path);

% Direct file paths
if exist(candidate, "file") == 2
    [fdir, fname, ext] = fileparts(candidate);
    if strcmpi(ext, ".mat")
        % Fast path: if a processed cache exists alongside, skip classification
        cache_path = fullfile(fdir, 'eq3D_processed.mat');
        if exist(cache_path, 'file') == 2 && ~strcmp(fname, 'eq3D_processed')
            cache_info = dir(cache_path);
            raw_info = dir(candidate);
            if cache_info.datenum >= raw_info.datenum
                source_path = string(cache_path);
                file_type = 'mat_eq3d';
                return
            end
        end
        % Detect: eq3D (2D) vs raw 4D .mat
        mat_type = local_classify_mat(candidate);
        source_path = string(candidate);
        file_type = mat_type;  % 'mat_eq3d' or 'mat_raw'
        return
    elseif strcmpi(ext, ".npy")
        % Fast path: check for cache
        cache_path = fullfile(fdir, 'eq3D_processed.mat');
        if exist(cache_path, 'file') == 2
            npy_info = dir(candidate);
            cache_info = dir(cache_path);
            if cache_info.datenum >= npy_info.datenum
                source_path = string(cache_path);
                file_type = 'mat_eq3d';
                return
            end
        end
        source_path = string(candidate);
        file_type = 'npy';
        return
    end
end

% Folder: look for processed cache first, then eq3D.mat, then .npy
if isfolder(candidate)
    % Cached result from raw import (highest priority)
    processed_path = fullfile(candidate, "eq3D_processed.mat");
    if exist(processed_path, "file") == 2
        source_path = string(processed_path);
        file_type = 'mat_eq3d';
        fprintf('  Found cached processed data: eq3D_processed.mat\n');
        return
    end

    eq3d_path = fullfile(candidate, "eq3D.mat");
    if exist(eq3d_path, "file") == 2
        source_path = string(eq3d_path);
        file_type = 'mat_eq3d';  % known eq3D
        return
    end

    % Look for .npy files
    npy_files = dir(fullfile(candidate, '*.npy'));
    if ~isempty(npy_files)
        % Pick the largest .npy file (most likely the EELS data)
        [~, biggest] = max([npy_files.bytes]);
        source_path = string(fullfile(candidate, npy_files(biggest).name));
        file_type = 'npy';
        fprintf('  No eq3D.mat found, using NPY: %s\n', npy_files(biggest).name);
        return
    end

    error("load_qe_dataset:NoDataFound", ...
        "Folder %s contains neither eq3D.mat nor .npy files.", candidate);
end

error("load_qe_dataset:PathNotFound", ...
    "Cannot resolve %s as a data file.", candidate);
end


function mat_type = local_classify_mat(mat_path)
%LOCAL_CLASSIFY_MAT  Determine if .mat is eq3D (2D) or raw 4D data.
%   Returns 'mat_eq3d' or 'mat_raw'.
try
    info = whos('-file', mat_path);
    names = {info.name};
    sizes = {info.size};
    classes = {info.class};

    % eq3D signature: has 'a3' (2D) and 'e'
    has_a3 = ismember('a3', names);
    has_e = ismember('e', names);
    if has_a3 && has_e
        a3_idx = find(strcmp(names, 'a3'), 1);
        if numel(sizes{a3_idx}) == 2
            mat_type = 'mat_eq3d';
            return
        end
    end

    % Check for EELS4D objects
    for i = 1:numel(classes)
        if strcmp(classes{i}, 'EELS4D')
            mat_type = 'mat_raw';
            return
        end
    end

    % Check for any 3D/4D numeric array
    for i = 1:numel(info)
        if numel(info(i).size) >= 3 && info(i).bytes > 1e6
            mat_type = 'mat_raw';
            return
        end
    end

    % Default: treat as eq3D attempt
    mat_type = 'mat_eq3d';
catch
    mat_type = 'mat_eq3d';
end
end


function intensity = local_match_energy_dimension(intensity, energy_meV, source_path)
if ~ismatrix(intensity) || isempty(intensity)
    error("load_qe_dataset:InvalidIntensity", ...
        "Intensity data loaded from %s must be a non-empty 2D matrix.", source_path);
end

energy_count = numel(energy_meV);
if size(intensity, 1) == energy_count
    return
end
if size(intensity, 2) == energy_count
    intensity = intensity.';
    return
end

error("load_qe_dataset:EnergyDimensionMismatch", ...
    "Energy axis length (%d) does not match any dimension of intensity from %s.", ...
    energy_count, source_path);
end


function dq_Ainv = local_resolve_dq(source_path, dqOverride)
persistent session_dq_Ainv

if isfinite(dqOverride) && dqOverride > 0
    dq_Ainv = dqOverride;
    session_dq_Ainv = dq_Ainv;
    return
end

lower_path = lower(char(source_path));
if contains(lower_path, "20w")
    dq_Ainv = 0.0025;
    return
end
if contains(lower_path, "10w")
    dq_Ainv = 0.005;
    return
end

if ~isempty(session_dq_Ainv) && isfinite(session_dq_Ainv) && session_dq_Ainv > 0
    dq_Ainv = session_dq_Ainv;
    return
end

if ~usejava("desktop")
    error("load_qe_dataset:DqRequired", ...
        "Unable to infer dq from %s and no dqOverride was provided.", source_path);
end

answer = inputdlg( ...
    {"Enter dq in 1/A per pixel:"}, ...
    "Set Momentum Step", ...
    [1 40], ...
    {"0.005"});

if isempty(answer)
    error("load_qe_dataset:DqCancelled", ...
        "dq entry was cancelled for %s.", source_path);
end

dq_Ainv = str2double(answer{1});
if ~isfinite(dq_Ainv) || dq_Ainv <= 0
    error("load_qe_dataset:DqInvalid", ...
        "dq must be a positive numeric value.");
end

session_dq_Ainv = dq_Ainv;
end


function label = local_make_label(source_path)
[parent_dir, file_name, ext] = fileparts(char(source_path));
[~, parent_name] = fileparts(parent_dir);

if strlength(parent_name) == 0
    label = string(file_name) + string(ext);
else
    label = string(parent_name) + filesep + string(file_name) + string(ext);
end
label = label + " [eq3d]";
end
