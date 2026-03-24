function qe = make_qe_struct(intensity, energy_meV, dq_Ainv, q_channel, metadata)
arguments
    intensity {mustBeNumeric}
    energy_meV {mustBeNumeric}
    dq_Ainv (1,1) double {mustBePositive}
    q_channel {mustBeNumeric} = []
    metadata.q_zero_index (1,1) double = NaN
    metadata.source_path string = ""
    metadata.source_kind string = ""
    metadata.label string = ""
    metadata.view_kind string = "physical"
    metadata.stage_name string = ""
    metadata.note string = ""
end

energy_meV = double(energy_meV(:));
intensity = double(intensity);

if size(intensity, 1) ~= numel(energy_meV)
    if size(intensity, 2) == numel(energy_meV)
        intensity = intensity.';
    else
        error("make_qe_struct:EnergySizeMismatch", ...
            "Intensity dimensions do not match the energy axis.");
    end
end

if isempty(q_channel) || numel(q_channel) ~= size(intensity, 2)
    q_channel = 1:size(intensity, 2);
end
q_channel = double(q_channel(:)).';
q_zero_index = local_resolve_q_zero_index(intensity, energy_meV, metadata.q_zero_index);

qe = struct();
qe.intensity = intensity;
qe.energy_meV = energy_meV;
qe.q_channel = q_channel;
qe.q_zero_index = q_zero_index;
qe.q_Ainv = ((1:size(intensity, 2)) - q_zero_index) * dq_Ainv;
qe.q_abs_Ainv = abs(qe.q_Ainv);
qe.dq_Ainv = dq_Ainv;
qe.source_path = char(metadata.source_path);
qe.source_kind = char(metadata.source_kind);
qe.label = char(metadata.label);
qe.view_kind = char(metadata.view_kind);
qe.stage_name = char(metadata.stage_name);
qe.note = char(metadata.note);
end


function q_zero_index = local_resolve_q_zero_index(intensity, energy_meV, q_zero_index)
nQ = size(intensity, 2);
if isfinite(q_zero_index) && q_zero_index >= 1 && q_zero_index <= nQ
    q_zero_index = round(q_zero_index);
    return
end

energy_meV = double(energy_meV(:));
if numel(energy_meV) > 1
    dE = median(abs(diff(energy_meV)), "omitnan");
else
    dE = NaN;
end
if ~isfinite(dE) || dE <= 0
    dE = 4;
end

window_half_width = max(20, 3 * dE);
mask = energy_meV >= -window_half_width & energy_meV <= window_half_width;
if ~any(mask)
    [~, nearest_zero] = min(abs(energy_meV));
    mask(nearest_zero) = true;
end

window_data = max(double(intensity(mask, :)), 0);
metric = sum(window_data, 1, "omitnan");
if isempty(metric) || ~any(isfinite(metric))
    q_zero_index = round((nQ + 1) / 2);
    return
end

[~, q_zero_index] = max(metric);
if isempty(q_zero_index) || ~isfinite(q_zero_index)
    q_zero_index = round((nQ + 1) / 2);
end
q_zero_index = min(max(round(q_zero_index), 1), nQ);
end
