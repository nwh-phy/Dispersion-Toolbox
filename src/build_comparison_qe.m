function comparison = build_comparison_qe(qeData, energyWindow, normalizationMode, smoothingMode, smoothingWidth, refAbsQMin, refAbsQMax, opts)
arguments
    qeData struct
    energyWindow (1,2) double = [-20 20]
    normalizationMode {mustBeTextScalar} = "zlp-area"
    smoothingMode {mustBeTextScalar} = "gaussian"
    smoothingWidth (1,1) double {mustBePositive} = 3
    refAbsQMin (1,1) double = NaN
    refAbsQMax (1,1) double = NaN
    opts.method {mustBeMember(opts.method, {'symmetric-avg', 'q-interp', 'q-profile'})} = 'q-profile'
    opts.smooth_scale (1,1) logical = true
    opts.idw_power (1,1) double {mustBePositive} = 2
end

local_validate_qe(qeData);

energyWindow = sort(double(energyWindow(:))).';
energy_axis = double(qeData.energy_meV(:));
q_axis = double(qeData.q_Ainv(:)).';
fit_window_mask = energy_axis >= energyWindow(1) & energy_axis <= energyWindow(2);
if ~any(fit_window_mask)
    [~, nearest_idx] = min(abs(energy_axis - mean(energyWindow)));
    fit_window_mask(nearest_idx) = true;
end

[ref_abs_q_min, ref_abs_q_max] = local_effective_reference_range(qeData, refAbsQMin, refAbsQMax);
[reference_q_mask_left, reference_q_mask_right, center_target_q_mask] = local_reference_masks(qeData, ref_abs_q_min, ref_abs_q_max);

% === Normalization (shared by both methods) ===
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

% === Reference spectra (still computed for diagnostics) ===
left_reference_spectrum = mean(normalized_original(:, reference_q_mask_left), 2);
right_reference_spectrum = mean(normalized_original(:, reference_q_mask_right), 2);

% === On-axis removal ===
off_axis_component = normalized_original;
fit_coeff_left = nan(1, size(normalized_original, 2));
fit_coeff_right = nan(1, size(normalized_original, 2));
scale_factors = ones(1, size(normalized_original, 2));
center_indices = find(center_target_q_mask);

switch opts.method
    case 'symmetric-avg'
        % Original method: fit α·S_left + β·S_right in ZLP window
        fit_matrix = [left_reference_spectrum(fit_window_mask), ...
                      right_reference_spectrum(fit_window_mask)];
        for qi = center_indices
            target_vector = normalized_original(fit_window_mask, qi);
            coeff = local_nonnegative_reference_fit(fit_matrix, target_vector);
            fit_coeff_left(qi) = coeff(1);
            fit_coeff_right(qi) = coeff(2);
            off_axis_component(:, qi) = max(coeff(1) * left_reference_spectrum ...
                                          + coeff(2) * right_reference_spectrum, 0);
        end

    case 'q-interp'
        % Plan A: IDW q-interpolation from individual reference channels
        % Plan C: Wider fit window (ZLP + sub-plasmon tail, exclude loss features)
        % Plan D: Smooth scale factors across q
        ref_mask = reference_q_mask_left | reference_q_mask_right;
        ref_indices = find(ref_mask);
        ref_q_values = q_axis(ref_indices);
        ref_spectra = normalized_original(:, ref_indices);

        % Build wider fit window: ZLP region + low-energy tail below plasmon onset
        % Auto-detect plasmon onset as the energy where off-axis spectra show
        % a prominent peak (use averaged reference as diagnostic)
        avg_ref = mean(ref_spectra, 2);
        avg_ref_smooth = smoothdata(avg_ref, 'gaussian', max(5, round(numel(avg_ref)*0.02)));
        positive_E_mask = energy_axis > max(energyWindow(2), 50);
        if any(positive_E_mask)
            [~, peak_idx_local] = max(avg_ref_smooth(positive_E_mask));
            peak_E_indices = find(positive_E_mask);
            plasmon_onset = energy_axis(peak_E_indices(max(1, peak_idx_local - 5)));
            plasmon_onset = max(plasmon_onset, 200);  % at least 200 meV
        else
            plasmon_onset = 300;
        end

        % Wider fit mask: ZLP window + pre-plasmon region [50, plasmon_onset]
        % Also include far tail > 3500 meV if available
        wide_fit_mask = fit_window_mask | ...
            (energy_axis >= 50 & energy_axis <= plasmon_onset * 0.8) | ...
            (energy_axis >= 3500);
        if nnz(wide_fit_mask) < 10
            wide_fit_mask = fit_window_mask;  % fallback to ZLP-only
        end

        for ki = 1:numel(center_indices)
            qi = center_indices(ki);
            q_c = q_axis(qi);

            % IDW weights: w_j = 1 / |q_c - q_j|^p
            distances = abs(q_c - ref_q_values);
            distances(distances < 1e-10) = 1e-10;
            weights = 1 ./ distances.^opts.idw_power;
            weights = weights / sum(weights);

            % Weighted reference spectrum for this specific q
            predicted_ref = ref_spectra * weights(:);

            % Fit single scale factor in wide window
            target_w = normalized_original(wide_fit_mask, qi);
            pred_w = predicted_ref(wide_fit_mask);
            valid = isfinite(target_w) & isfinite(pred_w) & pred_w > 0;
            if nnz(valid) >= 3
                alpha = max(0, sum(target_w(valid) .* pred_w(valid)) / ...
                             max(sum(pred_w(valid).^2), eps));
            else
                alpha = 1.0;
            end

            scale_factors(qi) = alpha;
            off_axis_component(:, qi) = max(alpha * predicted_ref, 0);

            % Store left/right decomposition for diagnostics
            left_w = weights(ref_q_values < 0);
            right_w = weights(ref_q_values >= 0);
            fit_coeff_left(qi) = sum(left_w);
            fit_coeff_right(qi) = sum(right_w);
        end

        % Plan D: smooth scale factors across center channels
        if opts.smooth_scale && numel(center_indices) >= 5
            raw_scales = scale_factors(center_indices);
            smooth_width = max(3, round(numel(center_indices) * 0.2));
            if mod(smooth_width, 2) == 0; smooth_width = smooth_width + 1; end
            smoothed_scales = smoothdata(raw_scales, 'sgolay', smooth_width);
            smoothed_scales = max(smoothed_scales, 0.01);

            for ki = 1:numel(center_indices)
                qi = center_indices(ki);
                old_scale = scale_factors(qi);
                new_scale = smoothed_scales(ki);
                if old_scale > eps
                    off_axis_component(:, qi) = off_axis_component(:, qi) * ...
                        (new_scale / old_scale);
                end
                scale_factors(qi) = new_scale;
            end
        end

    case 'q-profile'
        % === q-Profile Decomposition ===
        % Instead of matching spectra at different q (E-space fitting),
        % learn the on-axis beam q-profile shape from sub-plasmon energies,
        % then subtract a scaled copy at each energy.
        %
        % Physics: on-axis signal has FIXED q-shape f(q) (beam profile)
        %          but VARYING amplitude α(E).
        %          At each E: I(q,E) = α(E)·f(q) + I_offaxis(q,E)

        ref_mask = reference_q_mask_left | reference_q_mask_right;

        % Step 1: Auto-detect plasmon onset from reference channels
        avg_ref = mean(normalized_original(:, ref_mask), 2);
        avg_ref_s = smoothdata(avg_ref, 'gaussian', max(5, round(numel(avg_ref)*0.02)));
        pos_E = energy_axis > max(energyWindow(2), 50);
        if any(pos_E)
            [~, pk_loc] = max(avg_ref_s(pos_E));
            pk_idx = find(pos_E);
            plasmon_peak_E = energy_axis(pk_idx(pk_loc));
            template_E_max = max(min(plasmon_peak_E * 0.4, 400), 100);
        else
            template_E_max = 300;
        end

        % Step 2: Build q-profile template from sub-plasmon energies
        % In this range, the spectrum is dominated by ZLP tail (= on-axis)
        template_E_min = max(energyWindow(2) + 10, 30);
        tmask = energy_axis >= template_E_min & energy_axis <= template_E_max;
        if nnz(tmask) < 5
            tmask = energy_axis >= 20 & energy_axis <= 200;
        end

        q_profile_raw = mean(normalized_original(tmask, :), 1);
        q_profile_raw = max(q_profile_raw, 0);

        % Enforce q-symmetry: f(q) = f(-q)
        nq = numel(q_profile_raw);
        q_zero_idx = round((nq + 1) / 2);
        if isfield(qeData, 'q_zero_index') && isfinite(qeData.q_zero_index)
            q_zero_idx = round(qeData.q_zero_index);
        end
        for ii = 1:min(q_zero_idx-1, nq-q_zero_idx)
            left_val = q_profile_raw(q_zero_idx - ii);
            right_val = q_profile_raw(q_zero_idx + ii);
            avg_val = (left_val + right_val) / 2;
            q_profile_raw(q_zero_idx - ii) = avg_val;
            q_profile_raw(q_zero_idx + ii) = avg_val;
        end

        % Smooth the q-profile
        q_profile_smooth = smoothdata(q_profile_raw, 'gaussian', max(3, round(nq*0.03)));
        q_profile_smooth = max(q_profile_smooth, 0);

        % Step 3: For each energy, fit α(E) using center + near-reference channels
        % Use |q| <= ref_abs_q_max for fitting (on-axis tail is visible there)
        fit_q_mask = abs(q_axis) <= ref_abs_q_max;
        f_fit = q_profile_smooth(fit_q_mask);
        f_fit_sq = sum(f_fit.^2);

        alpha_E = zeros(numel(energy_axis), 1);
        for ei = 1:numel(energy_axis)
            row = normalized_original(ei, fit_q_mask);
            alpha_E(ei) = max(0, sum(row .* f_fit) / max(f_fit_sq, eps));
        end

        % Smooth α(E) to suppress noise in the subtraction
        smooth_span = max(7, round(numel(alpha_E) * 0.03));
        if mod(smooth_span, 2) == 0; smooth_span = smooth_span + 1; end
        alpha_smooth = smoothdata(alpha_E, 'sgolay', smooth_span);
        alpha_smooth = max(alpha_smooth, 0);

        % Zero out α in the ZLP core to avoid subtracting the elastic peak itself
        zlp_core = energy_axis >= energyWindow(1) & energy_axis <= energyWindow(2);
        alpha_smooth(zlp_core) = 0;

        % Step 4: Subtract on-axis component from ALL channels
        on_axis_2d = alpha_smooth * q_profile_smooth;
        off_axis_component = max(normalized_original - on_axis_2d, 0);

        % Store diagnostics
        scale_factors = q_profile_smooth;
        fit_coeff_left(:) = NaN;
        fit_coeff_right(:) = NaN;
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
    label=string(qeData.label) + " | off-axis component (" + opts.method + ")", ...
    view_kind="comparison", ...
    stage_name=string(qeData.stage_name), ...
    note="off-axis component reconstructed via " + opts.method);
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
comparison.method = opts.method;
comparison.scale_factors = scale_factors;
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
