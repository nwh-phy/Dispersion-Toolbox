function comparison = build_comparison_qe(qeData, energyWindow, normalizationMode, smoothingMode, smoothingWidth, refAbsQMin, refAbsQMax)
arguments
    qeData struct
    energyWindow (1,2) double = [-20 20]
    normalizationMode {mustBeTextScalar} = "zlp-area"
    smoothingMode {mustBeTextScalar} = "gaussian"
    smoothingWidth (1,1) double {mustBePositive} = 3
    refAbsQMin (1,1) double = NaN
    refAbsQMax (1,1) double = NaN
end

local_validate_qe(qeData);

energyWindow = sort(double(energyWindow(:))).';
energy_axis = double(qeData.energy_meV(:));
fit_window_mask = energy_axis >= energyWindow(1) & energy_axis <= energyWindow(2);
if ~any(fit_window_mask)
    [~, nearest_idx] = min(abs(energy_axis - mean(energyWindow)));
    fit_window_mask(nearest_idx) = true;
end

[ref_abs_q_min, ref_abs_q_max] = local_effective_reference_range(qeData, refAbsQMin, refAbsQMax);
[reference_q_mask_left, reference_q_mask_right, center_target_q_mask] = local_reference_masks(qeData, ref_abs_q_min, ref_abs_q_max);

normalized_original = double(qeData.intensity);
raw_window_metric = zeros(1, size(normalized_original, 2));
normalization_factor = ones(1, size(normalized_original, 2));

for qi = 1:size(normalized_original, 2)
    spectrum = double(qeData.intensity(:, qi));
    scale = local_normalization_scale(spectrum, energy_axis, fit_window_mask, normalizationMode);
    raw_window_metric(qi) = scale;

    if ~isfinite(scale) || scale <= eps
        positive_values = spectrum(spectrum > 0);
        if isempty(positive_values)
            scale = 1;
        else
            scale = max(max(positive_values), eps);
        end
    end

    normalization_factor(qi) = scale;
    normalized_original(:, qi) = max(spectrum, 0) ./ scale;
end

left_reference_spectrum = mean(normalized_original(:, reference_q_mask_left), 2);
right_reference_spectrum = mean(normalized_original(:, reference_q_mask_right), 2);

off_axis_component = normalized_original;
fit_coeff_left = nan(1, size(normalized_original, 2));
fit_coeff_right = nan(1, size(normalized_original, 2));

fit_matrix = [left_reference_spectrum(fit_window_mask), right_reference_spectrum(fit_window_mask)];
for qi = find(center_target_q_mask)
    target_vector = normalized_original(fit_window_mask, qi);
    coeff = local_nonnegative_reference_fit(fit_matrix, target_vector);
    fit_coeff_left(qi) = coeff(1);
    fit_coeff_right(qi) = coeff(2);
    off_axis_component(:, qi) = max(coeff(1) * left_reference_spectrum + coeff(2) * right_reference_spectrum, 0);
end

axis_excess = normalized_original - off_axis_component;
axis_excess(:, ~center_target_q_mask) = 0;
smoothed = local_apply_smoothing(off_axis_component, smoothingMode, smoothingWidth);

comparison = make_qe_struct( ...
    off_axis_component, ...
    qeData.energy_meV, ...
    qeData.dq_Ainv, ...
    qeData.q_channel, ...
    q_zero_index=local_get_q_zero_index(qeData), ...
    source_path=string(qeData.source_path), ...
    source_kind=string(qeData.source_kind), ...
    label=string(qeData.label) + " | normalized off-axis component, center reconstructed from manual symmetric off-axis reference bands, not physical intensity", ...
    view_kind="comparison", ...
    stage_name=string(qeData.stage_name), ...
    note="normalized off-axis component, center reconstructed from manual symmetric off-axis reference bands, not physical intensity");
comparison.smoothed_intensity = smoothed;
comparison.off_axis_component = off_axis_component;
comparison.normalized_original = normalized_original;
comparison.axis_excess = axis_excess;
comparison.normalization_mode = char(normalizationMode);
comparison.smoothing_mode = char(smoothingMode);
comparison.smoothing_width = smoothingWidth;
comparison.zlp_window = energyWindow;
comparison.reference_abs_q_min_Ainv = ref_abs_q_min;
comparison.reference_abs_q_max_Ainv = ref_abs_q_max;
comparison.left_reference_spectrum = left_reference_spectrum;
comparison.right_reference_spectrum = right_reference_spectrum;
comparison.reference_q_mask_left = reference_q_mask_left;
comparison.reference_q_mask_right = reference_q_mask_right;
comparison.center_target_q_mask = center_target_q_mask;
comparison.fit_coeff_left = fit_coeff_left;
comparison.fit_coeff_right = fit_coeff_right;
comparison.raw_window_metric = raw_window_metric;
comparison.normalization_factor = normalization_factor;
comparison.total_scale_factor = normalization_factor;
end


function smoothed = local_apply_smoothing(data, smoothingMode, smoothingWidth)
smoothed = double(data);
mode_name = lower(char(smoothingMode));
if strcmp(mode_name, "off")
    return
end

width = max(1, round(smoothingWidth));
for qi = 1:size(smoothed, 2)
    spectrum = smoothed(:, qi);
    switch mode_name
        case "gaussian"
            window = min(numel(spectrum), max(1, width));
            spectrum = smoothdata(spectrum, "gaussian", window);

        case "sgolay"
            window = min(numel(spectrum), max(3, width));
            if mod(window, 2) == 0
                window = max(3, window - 1);
            end
            if window < 3
                spectrum = smoothdata(spectrum, "movmean", max(1, width));
            else
                spectrum = smoothdata(spectrum, "sgolay", window);
            end

        otherwise
            error("build_comparison_qe:UnknownSmoothing", ...
                "Unknown smoothing mode '%s'.", smoothingMode);
    end
    smoothed(:, qi) = spectrum;
end
end


function local_validate_qe(qeData)
required_fields = {"intensity", "energy_meV", "q_channel", "q_Ainv", "dq_Ainv"};
for idx = 1:numel(required_fields)
    if ~isfield(qeData, required_fields{idx})
        error("build_comparison_qe:MissingField", ...
            "Input qeData is missing field '%s'.", required_fields{idx});
    end
end
end


function q_zero_index = local_get_q_zero_index(qeData)
if isfield(qeData, "q_zero_index") && isfinite(qeData.q_zero_index)
    q_zero_index = double(qeData.q_zero_index);
else
    q_zero_index = NaN;
end
end


function [ref_abs_q_min, ref_abs_q_max] = local_effective_reference_range(qeData, refAbsQMin, refAbsQMax)
q_abs_max = max(abs(double(qeData.q_Ainv(:))));
ref_abs_q_min = double(refAbsQMin);
ref_abs_q_max = double(refAbsQMax);

if ~isfinite(ref_abs_q_min) || ref_abs_q_min <= 0
    ref_abs_q_min = max(2 * double(qeData.dq_Ainv), double(qeData.dq_Ainv));
end
if ~isfinite(ref_abs_q_max) || ref_abs_q_max <= ref_abs_q_min
    ref_abs_q_max = max(2 * ref_abs_q_min, ref_abs_q_min + double(qeData.dq_Ainv));
end

ref_abs_q_max = min(ref_abs_q_max, q_abs_max);
if ref_abs_q_max <= ref_abs_q_min
    error("build_comparison_qe:InvalidReferenceRange", ...
        "Reference q range is invalid after clipping. Increase Ref |q| max or decrease Ref |q| min.");
end
end


function [left_mask, right_mask, center_mask] = local_reference_masks(qeData, ref_abs_q_min, ref_abs_q_max)
q_axis = double(qeData.q_Ainv(:)).';
left_mask = q_axis >= -ref_abs_q_max - eps & q_axis <= -ref_abs_q_min + eps;
right_mask = q_axis >= ref_abs_q_min - eps & q_axis <= ref_abs_q_max + eps;
center_mask = abs(q_axis) < ref_abs_q_min - eps;

if ~any(center_mask)
    center_mask = abs(q_axis) < ref_abs_q_min + eps;
end

if ~any(left_mask) || ~any(right_mask)
    error("build_comparison_qe:MissingReferenceColumns", ...
        "The selected reference q bands do not contain columns on both sides. Adjust Ref |q| min/max.");
end
if ~any(center_mask)
    error("build_comparison_qe:MissingCenterColumns", ...
        "The selected Ref |q| min does not cover any center columns. Increase Ref |q| min.");
end
end


function scale = local_normalization_scale(spectrum, energy_axis, fit_window_mask, normalizationMode)
switch lower(char(normalizationMode))
    case "zlp-area"
        normalization_window = max(spectrum(fit_window_mask), 0);
        if nnz(fit_window_mask) > 1
            scale = trapz(energy_axis(fit_window_mask), normalization_window);
        else
            scale = sum(normalization_window);
        end

    case "zlp-peak"
        normalization_window = spectrum(fit_window_mask);
        scale = max(normalization_window);

    otherwise
        error("build_comparison_qe:UnknownNormalization", ...
            "Unknown normalization mode '%s'.", normalizationMode);
end
end


function coeff = local_nonnegative_reference_fit(fit_matrix, target_vector)
fit_matrix = double(fit_matrix);
target_vector = double(target_vector(:));

valid = all(isfinite(fit_matrix), 2) & isfinite(target_vector);
fit_matrix = fit_matrix(valid, :);
target_vector = target_vector(valid);

if isempty(target_vector) || size(fit_matrix, 2) ~= 2
    coeff = [0; 0];
    return
end

if exist("lsqnonneg", "file") == 2
    coeff = lsqnonneg(fit_matrix, target_vector);
else
    coeff = max(fit_matrix \ target_vector, 0);
end

coeff = coeff(:);
if numel(coeff) < 2
    coeff(2, 1) = 0;
end
end
