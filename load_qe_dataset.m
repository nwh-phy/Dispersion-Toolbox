function dataset = load_qe_dataset(path, dqOverride)
arguments
    path {mustBeTextScalar}
    dqOverride (1,1) double = NaN
end

source_path = local_resolve_eq3d(string(path));
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


function source_path = local_resolve_eq3d(path)
candidate = char(path);

if isfolder(candidate)
    eq3d_path = fullfile(candidate, "eq3D.mat");
    if exist(eq3d_path, "file") == 2
        source_path = string(eq3d_path);
        return
    end
    error("load_qe_dataset:NoEq3D", ...
        "Folder %s does not contain eq3D.mat.", candidate);
end

if exist(candidate, "file") == 2
    [~, ~, ext] = fileparts(candidate);
    if strcmpi(ext, ".mat")
        source_path = string(candidate);
        return
    end
end

error("load_qe_dataset:PathNotFound", ...
    "Cannot resolve %s as an eq3D.mat file.", candidate);
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
