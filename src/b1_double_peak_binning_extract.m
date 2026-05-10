function result = b1_double_peak_binning_extract(qe_pp, qe_raw, opts)
%B1_DOUBLE_PEAK_BINNING_EXTRACT Mandatory B1 double-peak extraction.
%
% Low-SNR q channels are combined before fitting. Every extraction unit is
% then fit with exactly two B1 components. Failed double-peak fits are
% recorded instead of falling back to a single peak.

if nargin < 3 || isempty(opts)
    opts = struct();
end
opts = local_defaults(opts);

energy_axis = double(qe_pp.energy_meV(:));
q_axis = double(qe_pp.q_Ainv(:));
intensity = double(qe_pp.intensity);
raw_intensity = double(qe_raw.intensity);
if isempty(raw_intensity) || ~isequal(size(raw_intensity), size(intensity))
    raw_intensity = intensity;
end

energy_window = sort(double(opts.energy_window_meV(:)).');
energy_mask = energy_axis >= energy_window(1) & ...
    energy_axis <= energy_window(2);
if nnz(energy_mask) < 10
    error('b1_double_peak_binning_extract:InvalidEnergyWindow', ...
        'B1 energy window must contain at least 10 energy samples.');
end

q_bounds = sort(double(opts.q_range_Ainv(:)).');
q_candidates = find(q_axis >= q_bounds(1) & q_axis <= q_bounds(2) & ...
    abs(q_axis) >= opts.q_skip_Ainv & isfinite(q_axis));

noise_profile = local_noise_profile(energy_axis, intensity, ...
    q_axis, q_candidates, energy_mask, opts);
units = local_extraction_units(q_axis, noise_profile, opts);
unit_data = local_prepare_unit_data(energy_axis, q_axis, intensity, ...
    raw_intensity, units, energy_mask, opts);

switch opts.tracking_mode
    case 'independent_double_peak'
        result = local_extract_independent(unit_data, energy_axis, ...
            energy_window, opts);
    case 'propagated_double_peak'
        result = local_extract_propagated(unit_data, energy_axis, ...
            energy_window, noise_profile, opts);
    case 'windowed_branch_tracking'
        result = local_extract_windowed(unit_data, energy_axis, ...
            energy_window, opts);
    otherwise
        error('b1_double_peak_binning_extract:UnknownTrackingMode', ...
            'Unknown tracking mode "%s".', opts.tracking_mode);
end
result.noise_profile = noise_profile;
end


function unit_data = local_prepare_unit_data(energy_axis, q_axis, intensity, ...
    raw_intensity, units, energy_mask, opts)
unit_data = repmat(local_empty_unit_data(), 1, numel(units));
for ui = 1:numel(units)
    unit = units(ui);
    spectrum = mean(intensity(:, unit.q_indices), 2, 'omitnan');
    raw_spectrum = mean(raw_intensity(:, unit.q_indices), 2, 'omitnan');
    [fit_spectrum, fit_meta] = local_fit_spectrum(spectrum, opts, unit);
    fit_meta.peak_model = opts.peak_model;
    fit_meta.tracking_mode = opts.tracking_mode;
    guesses = local_initial_guesses(energy_axis, fit_spectrum, q_axis, ...
        unit, energy_mask, opts);

    unit_data(ui).unit = unit;
    unit_data(ui).spectrum = spectrum;
    unit_data(ui).raw_spectrum = raw_spectrum;
    unit_data(ui).fit_spectrum = fit_spectrum;
    unit_data(ui).fit_meta = fit_meta;
    unit_data(ui).guesses = guesses(:);
end
end


function data = local_empty_unit_data()
data = struct('unit', local_empty_unit(), 'spectrum', [], 'raw_spectrum', [], ...
    'fit_spectrum', [], 'fit_meta', struct(), 'guesses', []);
end


function result = local_extract_independent(unit_data, energy_axis, ...
    energy_window, opts)
result = local_empty_extraction_result(numel(unit_data));
for ui = 1:numel(unit_data)
    outcome = local_fit_double_unit(unit_data(ui), energy_axis, ...
        energy_window, opts, unit_data(ui).guesses);
    result = local_append_outcome(result, ui, outcome);
end
result = local_sort_extraction_result(result);
end


function result = local_extract_propagated(unit_data, energy_axis, ...
    energy_window, noise_profile, opts)
result = local_empty_extraction_result(numel(unit_data));
done = false(1, numel(unit_data));
q_values = arrayfun(@(x) x.unit.q_Ainv, unit_data);

side_sets = {find(q_values < 0), find(q_values > 0), find(q_values == 0)};
for si = 1:numel(side_sets)
    side_idx = side_sets{si};
    if isempty(side_idx)
        continue
    end
    [~, order] = sort(q_values(side_idx), 'ascend');
    side_idx = side_idx(order);
    [result, done] = local_propagate_one_side(result, done, unit_data, ...
        side_idx, energy_axis, energy_window, noise_profile, opts);
end

remaining = find(~done);
for i = remaining(:).'
    outcome = local_fit_double_unit(unit_data(i), energy_axis, ...
        energy_window, opts, unit_data(i).guesses);
    result = local_append_outcome(result, i, outcome);
end
result = local_sort_extraction_result(result);
end


function [result, done] = local_propagate_one_side(result, done, unit_data, ...
    side_idx, energy_axis, energy_window, noise_profile, opts)
seed_candidates = local_tracking_seed_candidates(unit_data, side_idx, ...
    noise_profile);
seed_idx = NaN;
seed_outcome = struct();
for ci = 1:numel(seed_candidates)
    candidate = seed_candidates(ci);
    outcome = local_fit_double_unit(unit_data(candidate), energy_axis, ...
        energy_window, opts, unit_data(candidate).guesses);
    if outcome.success
        seed_idx = candidate;
        seed_outcome = outcome;
        break
    end
end

if ~isfinite(seed_idx)
    for idx = side_idx(:).'
        outcome = local_fit_double_unit(unit_data(idx), energy_axis, ...
            energy_window, opts, unit_data(idx).guesses);
        result = local_append_outcome(result, idx, outcome);
        done(idx) = true;
    end
    return
end

result = local_append_outcome(result, seed_idx, seed_outcome);
done(seed_idx) = true;
seed_pos = find(side_idx == seed_idx, 1, 'first');
seed_energy = local_outcome_peak_energy(seed_outcome);

before = side_idx((seed_pos - 1):-1:1);
[result, done] = local_propagate_sequence(result, done, unit_data, before, ...
    energy_axis, energy_window, opts, seed_energy);

after = side_idx((seed_pos + 1):end);
[result, done] = local_propagate_sequence(result, done, unit_data, after, ...
    energy_axis, energy_window, opts, seed_energy);
end


function [result, done] = local_propagate_sequence(result, done, unit_data, ...
    indices, energy_axis, energy_window, opts, seed_energy)
previous_energy = seed_energy(:);
for idx = indices(:).'
    outcome = local_fit_double_unit(unit_data(idx), energy_axis, ...
        energy_window, opts, previous_energy);
    if outcome.success
        energy = local_outcome_peak_energy(outcome);
        shift = max(abs(energy(:) - previous_energy(:)));
        if shift > opts.max_tracking_shift_meV
            outcome = local_tracking_failure(unit_data(idx), previous_energy, ...
                sprintf('tracking shift %.3g meV', shift), opts);
        else
            previous_energy = energy(:);
        end
    end
    result = local_append_outcome(result, idx, outcome);
    done(idx) = true;
end
end


function candidates = local_tracking_seed_candidates(unit_data, side_idx, ...
    noise_profile)
scores = NaN(numel(side_idx), 1);
for i = 1:numel(side_idx)
    scores(i) = local_unit_seed_score(unit_data(side_idx(i)).unit, ...
        noise_profile);
end
if all(~isfinite(scores))
    q_abs = abs(arrayfun(@(x) x.unit.q_Ainv, unit_data(side_idx)));
    target = median(q_abs(isfinite(q_abs) & q_abs > 0), 'omitnan');
    if ~isfinite(target)
        target = median(q_abs, 'omitnan');
    end
    scores = -abs(q_abs(:) - target);
end
[~, order] = sort(scores, 'descend', 'MissingPlacement', 'last');
candidates = side_idx(order);
end


function score = local_unit_seed_score(unit, noise_profile)
score = NaN;
if isempty(noise_profile) || height(noise_profile) == 0
    return
end
rows = ismember(noise_profile.q_index, unit.q_indices(:));
if ~any(rows)
    return
end
signal = noise_profile.b1_signal(rows);
noise = noise_profile.b1_noise_mad(rows);
snr = signal ./ max(noise, eps);
score = median(snr(isfinite(snr)), 'omitnan');
if ~isfinite(score)
    score = NaN;
end
end


function result = local_extract_windowed(unit_data, energy_axis, energy_window, ...
    opts)
result = local_empty_extraction_result(numel(unit_data));
for ui = 1:numel(unit_data)
    outcome = local_fit_windowed_unit(unit_data(ui), energy_axis, ...
        energy_window, opts);
    result = local_append_outcome(result, ui, outcome);
end
result = local_sort_extraction_result(result);
end


function outcome = local_fit_windowed_unit(data, energy_axis, energy_window, opts)
unit = data.unit;
pred_lower = local_reference_energy(opts.reference_lower_points, unit.q_Ainv);
pred_upper = local_reference_energy(opts.reference_upper_points, unit.q_Ainv);
if ~isfinite(pred_lower) || ~isfinite(pred_upper)
    guesses = data.guesses(:);
    pred_lower = guesses(1);
    pred_upper = guesses(2);
end
if pred_lower > pred_upper
    tmp = pred_lower;
    pred_lower = pred_upper;
    pred_upper = tmp;
end

lower_half_width = local_tracking_window_half_width(unit, opts);
upper_half_width = local_upper_tracking_window_half_width(unit, opts);
lower_window = local_clamp_window(pred_lower + [-lower_half_width lower_half_width], ...
    energy_window);
upper_window = local_clamp_window(pred_upper + [-upper_half_width upper_half_width], ...
    energy_window);
guesses = [pred_lower; pred_upper];
if isempty(lower_window) || isempty(upper_window)
    if strcmp(opts.tracking_window_invalid_fallback, 'independent_double_peak')
        fallback = local_fit_double_unit(data, energy_axis, energy_window, ...
            opts, data.guesses);
        if fallback.success
            outcome = local_mark_outcome_repair(fallback, unit, ...
                'tracking_window_invalid_independent_double_peak', ...
                'tracking window outside B1 range; repaired with double-peak fit', ...
                guesses);
            return
        end
    end
    outcome = local_failed_double_outcome(unit, guesses, ...
        'tracking_window_invalid', 'tracking window outside B1 range', ...
        data.fit_meta);
    return
end

[lower_ok, lower_fit, lower_detail] = local_fit_single_branch(data, ...
    energy_axis, lower_window, pred_lower, opts);
[upper_ok, upper_fit, upper_detail] = local_fit_single_branch(data, ...
    energy_axis, upper_window, pred_upper, opts);

if ~lower_ok || ~upper_ok
    outcome = local_failed_double_outcome(unit, guesses, ...
        'windowed_branch_fit_failed', ...
        sprintf('lower: %s; upper: %s', lower_detail, upper_detail), ...
        data.fit_meta);
    return
end

[lower_energy, lower_ci] = local_peak_energy_and_ci(lower_fit);
[upper_energy, upper_ci] = local_peak_energy_and_ci(upper_fit);
lower_energy = lower_energy(1);
upper_energy = upper_energy(1);
if upper_energy < lower_energy
    tmp_fit = lower_fit;
    tmp_energy = lower_energy;
    tmp_ci = lower_ci;
    lower_fit = upper_fit;
    lower_energy = upper_energy;
    lower_ci = upper_ci;
    upper_fit = tmp_fit;
    upper_energy = tmp_energy;
    upper_ci = tmp_ci;
end
if upper_energy - lower_energy < opts.min_peak_separation_meV
    outcome = local_failed_double_outcome(unit, guesses, ...
        'collapsed_double_peak', ...
        sprintf('peak separation %.3g meV', upper_energy - lower_energy), ...
        data.fit_meta);
    return
end

upper_repair_source = 'raw_fit';
upper_repair_detail = '';
upper_original_energy = upper_energy;
if opts.upper_quality_retry
    [quality_ok, quality_detail] = local_upper_fit_quality(upper_fit, 1, ...
        upper_energy, opts);
    if ~quality_ok
        [retry_ok, retry_fit, retry_energy, retry_ci, retry_detail] = ...
            local_retry_upper_branch(data, energy_axis, energy_window, ...
            pred_upper, opts);
        if retry_ok
            [retry_quality_ok, retry_quality_detail] = ...
                local_upper_fit_quality(retry_fit, 1, retry_energy, opts);
            if retry_quality_ok && retry_energy - lower_energy >= ...
                    opts.min_peak_separation_meV
                upper_fit = retry_fit;
                upper_energy = retry_energy;
                upper_ci = retry_ci;
                upper_repair_source = 'upper_quality_retry';
                upper_repair_detail = sprintf( ...
                    'original upper rejected (%s); retry ok (%s)', ...
                    quality_detail, retry_detail);
            else
                if retry_quality_ok
                    retry_quality_detail = sprintf( ...
                        'retry separation %.3g meV', ...
                        retry_energy - lower_energy);
                end
                outcome = local_failed_double_outcome(unit, guesses, ...
                    'overbroad_upper_peak', ...
                    sprintf(['original upper rejected (%s); retry ', ...
                    'rejected (%s); %s'], quality_detail, ...
                    retry_quality_detail, retry_detail), data.fit_meta);
                return
            end
        else
            outcome = local_failed_double_outcome(unit, guesses, ...
                'overbroad_upper_peak', ...
                sprintf('original upper rejected (%s); retry failed (%s)', ...
                quality_detail, retry_detail), data.fit_meta);
            return
        end
    end
end

lower_row = local_point_row(unit, 1, 'b1_double_peak_lower', ...
    lower_fit, 1, lower_energy, lower_ci(1, :), energy_axis, ...
    data.raw_spectrum);
upper_row = local_point_row(unit, 2, 'b1_double_peak_upper', ...
    upper_fit, 1, upper_energy, upper_ci(1, :), energy_axis, ...
    data.raw_spectrum);
binning_row = local_binning_map_row(unit, guesses, true, 'ok', ...
    data.fit_meta);
repair_row = local_empty_repair_log_table();
if strcmp(upper_repair_source, 'upper_quality_retry')
    upper_row.repair_source(:) = {upper_repair_source};
    upper_row.original_energy_meV(:) = upper_original_energy;
    upper_row.repair_detail(:) = {upper_repair_detail};
    binning_row.repair_source(:) = {upper_repair_source};
    binning_row.repair_detail(:) = {upper_repair_detail};
    repair_row = local_repair_log_row(unit, 2, ...
        'b1_double_peak_upper', upper_repair_source, ...
        upper_original_energy, upper_energy, upper_repair_detail);
end

outcome = struct('success', true, 'lower_row', lower_row, ...
    'upper_row', upper_row, 'failure_row', local_empty_failures_table(), ...
    'binning_row', binning_row, 'fit', struct('lower_fit', lower_fit, ...
    'upper_fit', upper_fit), 'peak_energy', [lower_energy; upper_energy], ...
    'repair_row', repair_row);
end


function [success, fit, detail] = local_fit_single_branch(data, energy_axis, ...
    branch_window, guess, opts)
success = false;
fit = struct();
detail = '';
try
    fit = fit_loss_function(energy_axis, data.fit_spectrum, ...
        'E_min', branch_window(1), ...
        'E_max', branch_window(2), ...
        'max_peaks', 1, ...
        'min_prominence', opts.min_prominence, ...
        'smooth_width', opts.smooth_width, ...
        'initial_guesses', guess, ...
        'peak_model', opts.peak_model, ...
        'pre_subtracted', opts.pre_subtracted, ...
        'bootstrap_ci_samples', opts.bootstrap_ci_samples);
    fit = local_attach_fit_spectrum_metadata(fit, data.fit_meta, ...
        data.spectrum, data.fit_spectrum);
catch ME
    detail = ME.message;
    return
end
if fit.n_peaks ~= 1
    detail = sprintf('fit returned %d peaks', fit.n_peaks);
    return
end
success = true;
detail = 'ok';
end


function half_width = local_tracking_window_half_width(unit, opts)
if abs(unit.q_Ainv) >= opts.tracking_window_highq_abs_Ainv
    half_width = opts.tracking_window_highq_half_width_meV;
else
    half_width = opts.tracking_window_half_width_meV;
end
end


function half_width = local_upper_tracking_window_half_width(unit, opts)
if abs(unit.q_Ainv) >= opts.upper_tracking_window_highq_abs_Ainv
    half_width = opts.upper_tracking_window_highq_half_width_meV;
else
    half_width = opts.upper_tracking_window_half_width_meV;
end
end


function half_width = local_upper_retry_half_width(unit, opts)
if abs(unit.q_Ainv) >= opts.upper_retry_highq_abs_Ainv
    half_width = opts.upper_retry_highq_half_width_meV;
else
    half_width = opts.upper_retry_window_half_width_meV;
end
end


function [quality_ok, detail] = local_upper_fit_quality(fit, peak_idx, ...
    peak_energy, opts)
quality_ok = true;
reasons = {};
gamma = NaN;
if isfield(fit, 'gamma') && numel(fit.gamma) >= peak_idx
    gamma = fit.gamma(peak_idx);
end
if ~isfinite(gamma)
    reasons{end + 1} = 'nonfinite upper gamma'; %#ok<AGROW>
else
    gamma_over_E = gamma ./ max(abs(peak_energy), eps);
    if isfinite(opts.upper_max_gamma_over_E) && ...
            gamma_over_E > opts.upper_max_gamma_over_E
        reasons{end + 1} = sprintf('gamma/E %.3g > %.3g', ...
            gamma_over_E, opts.upper_max_gamma_over_E); %#ok<AGROW>
    end
    if isfinite(opts.upper_max_gamma_meV) && ...
            gamma > opts.upper_max_gamma_meV
        reasons{end + 1} = sprintf('gamma %.3g > %.3g meV', ...
            gamma, opts.upper_max_gamma_meV); %#ok<AGROW>
    end
end
if ~isempty(reasons)
    quality_ok = false;
    detail = strjoin(reasons, '; ');
else
    detail = 'ok';
end
end


function [success, fit, energy, ci, detail] = local_retry_upper_branch( ...
    data, energy_axis, energy_window, pred_upper, opts)
success = false;
fit = struct();
energy = NaN;
ci = [NaN NaN];
half_width = local_upper_retry_half_width(data.unit, opts);
retry_window = local_clamp_window(pred_upper + [-half_width half_width], ...
    energy_window);
if isempty(retry_window)
    detail = 'upper retry window outside B1 range';
    return
end
[success, fit, detail] = local_fit_single_branch(data, energy_axis, ...
    retry_window, pred_upper, opts);
if ~success
    return
end
[energy, ci] = local_peak_energy_and_ci(fit);
energy = energy(1);
ci = ci(1, :);
detail = sprintf('upper retry window %.3g-%.3g meV', ...
    retry_window(1), retry_window(2));
end


function win = local_clamp_window(win, energy_window)
win = sort(double(win(:)).');
win(1) = max(win(1), energy_window(1));
win(2) = min(win(2), energy_window(2));
if numel(win) ~= 2 || win(2) - win(1) < 20
    win = [];
end
end


function energy = local_reference_energy(points, q_value)
energy = NaN;
if isempty(points) || ~istable(points) || ...
        ~all(ismember({'q_Ainv', 'energy_meV'}, points.Properties.VariableNames))
    return
end
q = double(points.q_Ainv(:));
E = double(points.energy_meV(:));
valid = isfinite(q) & isfinite(E);
q = q(valid);
E = E(valid);
if isempty(q)
    return
end
[q_unique, ~, group] = unique(q);
E_unique = splitapply(@median, E, group);
if numel(E_unique) >= 3
    E_unique = smoothdata(E_unique, 'movmedian', 3);
end
if numel(q_unique) == 1
    energy = E_unique(1);
else
    energy = interp1(q_unique, E_unique, q_value, 'pchip', 'extrap');
end
end


function result = local_empty_extraction_result(n_units)
result = struct();
result.lower_points = local_empty_points_table();
result.upper_points = local_empty_points_table();
result.combined_points = local_empty_points_table();
result.binning_map = local_empty_binning_map();
result.fit_failures = local_empty_failures_table();
result.repair_log = local_empty_repair_log_table();
result.fit_details = cell(n_units, 1);
end


function result = local_append_outcome(result, ui, outcome)
if outcome.success
    result.lower_points = [result.lower_points; outcome.lower_row]; %#ok<AGROW>
    result.upper_points = [result.upper_points; outcome.upper_row]; %#ok<AGROW>
    result.combined_points = [result.combined_points; outcome.lower_row; ...
        outcome.upper_row]; %#ok<AGROW>
    result.fit_details{ui} = outcome.fit;
else
    result.fit_failures = [result.fit_failures; outcome.failure_row]; %#ok<AGROW>
end
if isfield(outcome, 'repair_row') && ~isempty(outcome.repair_row)
    result.repair_log = [result.repair_log; outcome.repair_row]; %#ok<AGROW>
end
result.binning_map = [result.binning_map; outcome.binning_row]; %#ok<AGROW>
end


function result = local_sort_extraction_result(result)
result.lower_points = sortrows(result.lower_points, {'q_Ainv', 'energy_meV'});
result.upper_points = sortrows(result.upper_points, {'q_Ainv', 'energy_meV'});
result.combined_points = sortrows(result.combined_points, ...
    {'q_Ainv', 'branch', 'energy_meV'});
end


function outcome = local_fit_double_unit(data, energy_axis, energy_window, ...
    opts, guesses)
unit = data.unit;
guesses = local_two_guess_vector(guesses, data.guesses);
try
    fit = fit_loss_function(energy_axis, data.fit_spectrum, ...
        'E_min', energy_window(1), ...
        'E_max', energy_window(2), ...
        'max_peaks', 2, ...
        'min_prominence', opts.min_prominence, ...
        'smooth_width', opts.smooth_width, ...
        'initial_guesses', guesses(:), ...
        'peak_model', opts.peak_model, ...
        'pre_subtracted', opts.pre_subtracted, ...
        'bootstrap_ci_samples', opts.bootstrap_ci_samples);
    fit = local_attach_fit_spectrum_metadata(fit, data.fit_meta, ...
        data.spectrum, data.fit_spectrum);
catch ME
    outcome = local_failed_double_outcome(unit, guesses, ...
        'double_peak_fit_failed', ME.message, data.fit_meta);
    return
end

if fit.n_peaks ~= 2
    outcome = local_failed_double_outcome(unit, guesses, 'not_two_peaks', ...
        sprintf('fit returned %d peaks', fit.n_peaks), data.fit_meta);
    return
end

[peak_energy, peak_ci] = local_peak_energy_and_ci(fit);
[peak_energy_sorted, sort_idx] = sort(peak_energy(:));
if numel(peak_energy_sorted) ~= 2 || any(~isfinite(peak_energy_sorted))
    outcome = local_failed_double_outcome(unit, guesses, ...
        'nonfinite_double_peak_energy', ...
        'one or both peak energies are non-finite', data.fit_meta);
    return
end
if diff(peak_energy_sorted) < opts.min_peak_separation_meV
    outcome = local_failed_double_outcome(unit, guesses, ...
        'collapsed_double_peak', sprintf('peak separation %.3g meV', ...
        diff(peak_energy_sorted)), data.fit_meta);
    return
end

lower_idx = sort_idx(1);
upper_idx = sort_idx(2);
lower_row = local_point_row(unit, 1, 'b1_double_peak_lower', ...
    fit, lower_idx, peak_energy(lower_idx), peak_ci(lower_idx, :), ...
    energy_axis, data.raw_spectrum);
upper_row = local_point_row(unit, 2, 'b1_double_peak_upper', ...
    fit, upper_idx, peak_energy(upper_idx), peak_ci(upper_idx, :), ...
    energy_axis, data.raw_spectrum);

outcome = struct('success', true, 'lower_row', lower_row, ...
    'upper_row', upper_row, 'failure_row', local_empty_failures_table(), ...
    'binning_row', local_binning_map_row(unit, guesses, true, 'ok', ...
    data.fit_meta), 'fit', fit, 'peak_energy', peak_energy_sorted(:), ...
    'repair_row', local_empty_repair_log_table());
end


function guesses = local_two_guess_vector(guesses, fallback)
guesses = double(guesses(:));
guesses = guesses(isfinite(guesses));
if numel(guesses) < 2
    guesses = double(fallback(:));
end
if numel(guesses) < 2
    guesses = [NaN; NaN];
else
    guesses = sort(guesses(1:2));
end
end


function energy = local_outcome_peak_energy(outcome)
energy = outcome.peak_energy(:);
if numel(energy) > 2
    energy = energy(1:2);
end
energy = sort(energy(:));
end


function outcome = local_tracking_failure(data, guesses, detail, opts)
outcome = local_failed_double_outcome(data.unit, guesses, ...
    'tracking_shift_exceeded', detail, data.fit_meta);
outcome.binning_row.fit_status{1} = 'tracking_shift_exceeded';
outcome.binning_row.initial_guess_lower_meV(1) = guesses(1);
outcome.binning_row.initial_guess_upper_meV(1) = guesses(2);
if isfield(opts, 'tracking_mode')
    outcome.binning_row.tracking_mode{1} = opts.tracking_mode;
end
end


function outcome = local_failed_double_outcome(unit, guesses, status, detail, ...
    fit_meta)
guesses = local_two_guess_vector(guesses, [NaN; NaN]);
outcome = struct('success', false, ...
    'lower_row', local_empty_points_table(), ...
    'upper_row', local_empty_points_table(), ...
    'failure_row', local_failure_row(unit, status, detail, fit_meta), ...
    'binning_row', local_binning_map_row(unit, guesses, false, status, ...
    fit_meta), 'fit', [], 'peak_energy', [NaN; NaN], ...
    'repair_row', local_empty_repair_log_table());
end


function outcome = local_mark_outcome_repair(outcome, unit, source, detail, ...
    guesses)
if ~outcome.success
    return
end
outcome.lower_row.repair_source(:) = {source};
outcome.upper_row.repair_source(:) = {source};
outcome.lower_row.repair_detail(:) = {detail};
outcome.upper_row.repair_detail(:) = {detail};
outcome.lower_row.original_energy_meV(:) = outcome.lower_row.energy_meV;
outcome.upper_row.original_energy_meV(:) = outcome.upper_row.energy_meV;
outcome.binning_row.fit_status(:) = {'ok_repaired_tracking_window_invalid'};
outcome.binning_row.repair_source(:) = {source};
outcome.binning_row.repair_detail(:) = {detail};
outcome.binning_row.initial_guess_lower_meV(:) = guesses(1);
outcome.binning_row.initial_guess_upper_meV(:) = guesses(2);
outcome.repair_row = [ ...
    local_repair_log_row(unit, 1, 'b1_double_peak_lower', source, NaN, ...
    outcome.lower_row.energy_meV(1), detail); ...
    local_repair_log_row(unit, 2, 'b1_double_peak_upper', source, NaN, ...
    outcome.upper_row.energy_meV(1), detail)];
end


function opts = local_defaults(opts)
opts = local_set_default(opts, 'energy_window_meV', [300 2100]);
opts = local_set_default(opts, 'q_range_Ainv', [-0.15 0.15]);
opts = local_set_default(opts, 'q_skip_Ainv', 0.005);
opts = local_set_default(opts, 'bin_size', 3);
opts = local_set_default(opts, 'noise_threshold', NaN);
opts = local_set_default(opts, 'low_q_no_bin_abs_Ainv', 0.05);
opts = local_set_default(opts, 'high_q_force_bin_abs_Ainv', Inf);
opts = local_set_default(opts, 'edge_fraction', 0.18);
opts = local_set_default(opts, 'peak_model', 'fano');
opts = local_set_default(opts, 'pre_subtracted', false);
opts = local_set_default(opts, 'min_prominence', 0.02);
opts = local_set_default(opts, 'smooth_width', 1);
opts = local_set_default(opts, 'bootstrap_ci_samples', 0);
opts = local_set_default(opts, 'old_branch_points', table());
opts = local_set_default(opts, 'fallback_split_meV', 180);
opts = local_set_default(opts, 'fallback_split_candidates_meV', opts.fallback_split_meV);
opts = local_set_default(opts, 'min_peak_separation_meV', 10);
opts = local_set_default(opts, 'initial_peak_margin_meV', 25);
opts = local_set_default(opts, 'tracking_mode', 'independent_double_peak');
opts = local_set_default(opts, 'max_tracking_shift_meV', 180);
opts = local_set_default(opts, 'tracking_window_half_width_meV', 220);
opts = local_set_default(opts, 'tracking_window_highq_half_width_meV', 300);
opts = local_set_default(opts, 'tracking_window_highq_abs_Ainv', 0.09);
opts = local_set_default(opts, 'tracking_window_invalid_fallback', 'fail');
opts = local_set_default(opts, 'upper_tracking_window_half_width_meV', ...
    opts.tracking_window_half_width_meV);
opts = local_set_default(opts, 'upper_tracking_window_highq_half_width_meV', ...
    opts.tracking_window_highq_half_width_meV);
opts = local_set_default(opts, 'upper_tracking_window_highq_abs_Ainv', ...
    opts.tracking_window_highq_abs_Ainv);
opts = local_set_default(opts, 'upper_quality_retry', false);
opts = local_set_default(opts, 'upper_max_gamma_over_E', Inf);
opts = local_set_default(opts, 'upper_max_gamma_meV', Inf);
opts = local_set_default(opts, 'upper_retry_window_half_width_meV', 180);
opts = local_set_default(opts, 'upper_retry_highq_half_width_meV', 240);
opts = local_set_default(opts, 'upper_retry_highq_abs_Ainv', 0.09);
opts = local_set_default(opts, 'reference_lower_points', table());
opts = local_set_default(opts, 'reference_upper_points', table());
opts = local_set_default(opts, 'fit_denoise_method', 'none');
opts = local_set_default(opts, 'fit_denoise_profile', 'global');
opts = local_set_default(opts, 'fit_denoise_window', 11);
opts = local_set_default(opts, 'fit_denoise_low_window', opts.fit_denoise_window);
opts = local_set_default(opts, 'fit_denoise_high_window', opts.fit_denoise_window);
opts = local_set_default(opts, 'fit_denoise_q_start_Ainv', 0.07);
opts = local_set_default(opts, 'fit_denoise_q_end_Ainv', 0.15);
opts = local_set_default(opts, 'fit_denoise_order', 3);
opts.bin_size = max(1, round(double(opts.bin_size)));
opts.smooth_width = max(1, round(double(opts.smooth_width)));
opts.low_q_no_bin_abs_Ainv = max(0, double(opts.low_q_no_bin_abs_Ainv));
opts.high_q_force_bin_abs_Ainv = double(opts.high_q_force_bin_abs_Ainv);
opts.fit_denoise_method = lower(strtrim(char(string(opts.fit_denoise_method))));
opts.fit_denoise_profile = lower(strtrim(char(string(opts.fit_denoise_profile))));
opts.peak_model = lower(strtrim(char(string(opts.peak_model))));
opts.tracking_mode = lower(strtrim(char(string(opts.tracking_mode))));
opts.tracking_window_invalid_fallback = lower(strtrim( ...
    char(string(opts.tracking_window_invalid_fallback))));
opts.fallback_split_candidates_meV = double(opts.fallback_split_candidates_meV(:).');
opts.fallback_split_candidates_meV = opts.fallback_split_candidates_meV( ...
    isfinite(opts.fallback_split_candidates_meV) & opts.fallback_split_candidates_meV > 0);
if isempty(opts.fallback_split_candidates_meV)
    opts.fallback_split_candidates_meV = opts.fallback_split_meV;
end
opts.max_tracking_shift_meV = max(1, double(opts.max_tracking_shift_meV));
opts.tracking_window_half_width_meV = max(1, double(opts.tracking_window_half_width_meV));
opts.tracking_window_highq_half_width_meV = max(opts.tracking_window_half_width_meV, ...
    double(opts.tracking_window_highq_half_width_meV));
opts.tracking_window_highq_abs_Ainv = max(0, double(opts.tracking_window_highq_abs_Ainv));
opts.upper_tracking_window_half_width_meV = max(1, ...
    double(opts.upper_tracking_window_half_width_meV));
opts.upper_tracking_window_highq_half_width_meV = max( ...
    opts.upper_tracking_window_half_width_meV, ...
    double(opts.upper_tracking_window_highq_half_width_meV));
opts.upper_tracking_window_highq_abs_Ainv = max(0, ...
    double(opts.upper_tracking_window_highq_abs_Ainv));
opts.upper_quality_retry = logical(opts.upper_quality_retry);
opts.upper_max_gamma_over_E = double(opts.upper_max_gamma_over_E);
opts.upper_max_gamma_meV = double(opts.upper_max_gamma_meV);
opts.upper_retry_window_half_width_meV = max(1, ...
    double(opts.upper_retry_window_half_width_meV));
opts.upper_retry_highq_half_width_meV = max( ...
    opts.upper_retry_window_half_width_meV, ...
    double(opts.upper_retry_highq_half_width_meV));
opts.upper_retry_highq_abs_Ainv = max(0, ...
    double(opts.upper_retry_highq_abs_Ainv));
opts.fit_denoise_window = local_odd_window(opts.fit_denoise_window);
opts.fit_denoise_low_window = local_odd_window(opts.fit_denoise_low_window);
opts.fit_denoise_high_window = local_odd_window(opts.fit_denoise_high_window);
opts.fit_denoise_q_start_Ainv = max(0, double(opts.fit_denoise_q_start_Ainv));
opts.fit_denoise_q_end_Ainv = max(opts.fit_denoise_q_start_Ainv, ...
    double(opts.fit_denoise_q_end_Ainv));
opts.fit_denoise_order = max(1, round(double(opts.fit_denoise_order)));
min_window = min([opts.fit_denoise_window, opts.fit_denoise_low_window, ...
    opts.fit_denoise_high_window]);
if opts.fit_denoise_order >= min_window
    opts.fit_denoise_order = max(1, min_window - 2);
end
end


function window = local_odd_window(value)
window = max(3, round(double(value)));
if mod(window, 2) == 0
    window = window + 1;
end
end


function s = local_set_default(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end


function profile = local_noise_profile(energy_axis, intensity, q_axis, ...
    q_candidates, energy_mask, opts)
energy_win = energy_axis(energy_mask);
n_edge = max(3, round(numel(energy_win) * opts.edge_fraction));
edge_local = false(size(energy_win));
edge_local(1:n_edge) = true;
edge_local(end-n_edge+1:end) = true;

n = numel(q_candidates);
q_index = q_candidates(:);
q = q_axis(q_index);
signal = NaN(n, 1);
noise = NaN(n, 1);
noise_ratio = NaN(n, 1);

for i = 1:n
    y = double(intensity(energy_mask, q_index(i)));
    edge = y(edge_local);
    baseline = median(edge(isfinite(edge)), 'omitnan');
    if ~isfinite(baseline)
        baseline = 0;
    end
    centered_edge = edge - baseline;
    noise(i) = 1.4826 * median(abs(centered_edge(isfinite(centered_edge))), ...
        'omitnan');
    signal(i) = max(y, [], 'omitnan') - baseline;
    if ~isfinite(noise(i))
        noise(i) = NaN;
    end
    if isfinite(signal(i)) && signal(i) > eps
        noise_ratio(i) = noise(i) ./ signal(i);
    end
end

threshold = opts.noise_threshold;
if isnan(threshold)
    threshold = local_auto_noise_threshold(noise_ratio);
    threshold_source = 'auto_per_dataset';
else
    threshold_source = 'manual_override';
end
[use_binning, noise_exceeds_threshold, low_q_protected, high_q_force_binning] = ...
    local_threshold_binning_mask(q, noise_ratio, threshold, ...
    opts.low_q_no_bin_abs_Ainv, opts.high_q_force_bin_abs_Ainv);

profile = table(q_index, q, abs(q), signal, noise, noise_ratio, ...
    repmat(threshold, n, 1), ...
    repmat(opts.low_q_no_bin_abs_Ainv, n, 1), ...
    repmat(opts.high_q_force_bin_abs_Ainv, n, 1), ...
    repmat({threshold_source}, n, 1), ...
    noise_exceeds_threshold, low_q_protected, high_q_force_binning, ...
    use_binning, ...
    'VariableNames', {'q_index', 'q_Ainv', 'q_abs_Ainv', ...
    'b1_signal', 'b1_noise_mad', 'noise_ratio', ...
    'noise_threshold', 'low_q_no_bin_abs_Ainv', ...
    'high_q_force_bin_abs_Ainv', ...
    'noise_threshold_source', 'noise_exceeds_threshold', ...
    'low_q_no_bin_protected', 'high_q_force_binning', ...
    'use_binning'});
end


function threshold = local_auto_noise_threshold(noise_ratio)
valid = noise_ratio(isfinite(noise_ratio));
if isempty(valid)
    threshold = Inf;
    return
end
center = median(valid, 'omitnan');
spread = 1.4826 * median(abs(valid - center), 'omitnan');
if ~isfinite(spread) || spread <= 0
    threshold = center * 1.5;
else
    threshold = center + 2.5 * spread;
end
threshold = max(threshold, prctile(valid, 70));
end


function [use_binning, noise_exceeds_threshold, low_q_protected, ...
    high_q_force_binning] = local_threshold_binning_mask(q, noise_ratio, ...
    threshold, low_q_no_bin_abs, high_q_force_bin_abs)
use_binning = false(size(q));
noise_exceeds_threshold = false(size(q));
low_q_protected = false(size(q));
high_q_force_binning = false(size(q));
if isfinite(threshold)
    noise_exceeds_threshold = isfinite(noise_ratio) & noise_ratio > threshold;
end
if isfinite(low_q_no_bin_abs) && low_q_no_bin_abs > 0
    low_q_protected = isfinite(q) & abs(q) <= low_q_no_bin_abs;
end
if isfinite(high_q_force_bin_abs) && high_q_force_bin_abs > 0
    high_q_force_binning = isfinite(q) & abs(q) >= high_q_force_bin_abs;
end
use_binning = (noise_exceeds_threshold | high_q_force_binning) & ...
    ~low_q_protected;
end


function units = local_extraction_units(q_axis, noise_profile, opts)
units = local_empty_unit();
units = units([]);

q_index = noise_profile.q_index;
use_binning = logical(noise_profile.use_binning);
direct_idx = q_index(~use_binning);
for i = 1:numel(direct_idx)
    idx = direct_idx(i);
    units(end+1) = local_make_unit(q_axis, idx, ... %#ok<AGROW>
        'single_q_direct', opts.bin_size);
end

for sign_value = [-1 1]
    if sign_value < 0
        side = noise_profile.q_Ainv < 0;
    else
        side = noise_profile.q_Ainv > 0;
    end
    side_rows = find(side(:));
    if isempty(side_rows)
        continue
    end
    [~, order] = sort(q_axis(q_index(side_rows)), 'ascend');
    side_rows = side_rows(order);
    run_idx = [];
    for ri = 1:numel(side_rows)
        row_idx = side_rows(ri);
        if use_binning(row_idx)
            run_idx(end+1) = q_index(row_idx); %#ok<AGROW>
        else
            units = local_append_binned_run(units, q_axis, run_idx, opts);
            run_idx = [];
        end
    end
    units = local_append_binned_run(units, q_axis, run_idx, opts);
end

[~, order] = sort([units.q_Ainv]);
units = units(order);
end


function units = local_append_binned_run(units, q_axis, run_idx, opts)
if isempty(run_idx)
    return
end
start_idx = 1;
while start_idx <= numel(run_idx)
    stop_idx = min(start_idx + opts.bin_size - 1, numel(run_idx));
    group = run_idx(start_idx:stop_idx);
    units(end+1) = local_make_unit(q_axis, group, ... %#ok<AGROW>
        sprintf('combined_q_binning_%d', opts.bin_size), opts.bin_size);
    start_idx = stop_idx + 1;
end
end


function unit = local_empty_unit()
unit = struct('q_Ainv', NaN, 'q_abs_Ainv', NaN, 'q_indices', [], ...
    'source_mode', '', 'source_q_count', 0, 'source_q_Ainv', '', ...
    'source_q_index', '', 'bin_size_requested', NaN);
end


function unit = local_make_unit(q_axis, q_indices, source_mode, bin_size)
q_indices = q_indices(:).';
q_values = q_axis(q_indices);
unit = local_empty_unit();
unit.q_Ainv = mean(q_values, 'omitnan');
unit.q_abs_Ainv = abs(unit.q_Ainv);
unit.q_indices = q_indices;
unit.source_mode = source_mode;
unit.source_q_count = numel(q_indices);
unit.source_q_Ainv = local_join_numbers(q_values);
unit.source_q_index = local_join_numbers(q_indices);
unit.bin_size_requested = bin_size;
end


function [fit_spectrum, meta] = local_fit_spectrum(spectrum, opts, unit)
method = lower(strtrim(char(string(opts.fit_denoise_method))));
window = local_fit_denoise_window_for_unit(opts, unit);
fit_spectrum = double(spectrum(:));
meta = struct();
meta.fit_denoise_method = method;
meta.fit_denoise_window = window;
meta.fit_denoise_order = opts.fit_denoise_order;
meta.fit_denoise_profile = opts.fit_denoise_profile;
meta.fit_denoise_q_start_Ainv = opts.fit_denoise_q_start_Ainv;
meta.fit_denoise_q_end_Ainv = opts.fit_denoise_q_end_Ainv;

switch method
    case {'none', 'off', ''}
        method = 'none';
    case {'sgolay', 'savgol', 'savitzky-golay'}
        fit_spectrum = local_sgolay_denoise(fit_spectrum, ...
            window, opts.fit_denoise_order);
        method = 'sgolay';
    case 'gaussian'
        fit_spectrum = smoothdata(fit_spectrum, 'gaussian', ...
            window);
        method = 'gaussian';
    otherwise
        error('b1_double_peak_binning_extract:UnknownFitDenoiseMethod', ...
            'Unknown fit denoise method "%s".', opts.fit_denoise_method);
end

meta.fit_denoise_method = method;
if strcmp(method, 'none')
    meta.fit_spectrum_source = 'original';
else
    meta.fit_spectrum_source = 'denoised';
end


function window = local_fit_denoise_window_for_unit(opts, unit)
profile = lower(strtrim(char(string(opts.fit_denoise_profile))));
switch profile
    case {'global', 'none', ''}
        window = opts.fit_denoise_window;
    case {'adaptive_absq', 'adaptive-q', 'adaptive_q', 'absq'}
        q_abs = abs(double(unit.q_Ainv));
        q0 = opts.fit_denoise_q_start_Ainv;
        q1 = opts.fit_denoise_q_end_Ainv;
        if ~isfinite(q_abs)
            q_abs = 0;
        end
        if q1 <= q0
            t = double(q_abs >= q0);
        else
            t = min(max((q_abs - q0) / (q1 - q0), 0), 1);
            t = t * t * (3 - 2 * t);
        end
        window = opts.fit_denoise_low_window + ...
            t * (opts.fit_denoise_high_window - opts.fit_denoise_low_window);
        window = local_odd_window(window);
    otherwise
        error('b1_double_peak_binning_extract:UnknownFitDenoiseProfile', ...
            'Unknown fit denoise profile "%s".', opts.fit_denoise_profile);
end
end
delta = fit_spectrum - double(spectrum(:));
meta.fit_spectrum_delta_rms = sqrt(mean(delta(isfinite(delta)).^2, ...
    'omitnan'));
if ~isfinite(meta.fit_spectrum_delta_rms)
    meta.fit_spectrum_delta_rms = 0;
end
end


function y = local_sgolay_denoise(y, window, order)
window = min(window, numel(y));
if mod(window, 2) == 0
    window = window - 1;
end
if window <= order || window < 3
    return
end
try
    y = sgolayfilt(y, order, window);
catch
    y = smoothdata(y, 'sgolay', window);
end
end


function fit = local_attach_fit_spectrum_metadata(fit, meta, raw_spectrum, ...
    fit_spectrum)
fit.fit_spectrum_source = meta.fit_spectrum_source;
fit.fit_denoise_method = meta.fit_denoise_method;
fit.fit_denoise_window = meta.fit_denoise_window;
fit.fit_denoise_order = meta.fit_denoise_order;
fit.fit_denoise_profile = meta.fit_denoise_profile;
fit.fit_denoise_q_start_Ainv = meta.fit_denoise_q_start_Ainv;
fit.fit_denoise_q_end_Ainv = meta.fit_denoise_q_end_Ainv;
fit.fit_spectrum_delta_rms = meta.fit_spectrum_delta_rms;
fit.peak_model = meta.peak_model;
fit.tracking_mode = meta.tracking_mode;
fit.raw_unit_spectrum = double(raw_spectrum(:));
fit.fit_input_spectrum = double(fit_spectrum(:));
end


function guesses = local_initial_guesses(energy_axis, spectrum, q_axis, ...
    unit, energy_mask, opts)
energy_win = energy_axis(energy_mask);
y = double(spectrum(energy_mask));
baseline = median(y, 'omitnan');
if ~isfinite(baseline)
    baseline = 0;
end
y_det = y - baseline;

guesses = local_two_local_peak_guesses(energy_win, y_det, opts);
if numel(guesses) == 2
    return
end

center = local_old_branch_center(opts.old_branch_points, unit.q_Ainv);
if ~isfinite(center)
    [~, max_idx] = max(y_det);
    if ~isempty(max_idx)
        center = energy_win(max_idx);
    end
end
if ~isfinite(center)
    center = mean(opts.energy_window_meV);
end

fallback_split = median(opts.fallback_split_candidates_meV, 'omitnan');
if ~isfinite(fallback_split)
    fallback_split = opts.fallback_split_meV;
end
half_split = fallback_split / 2;
win = sort(double(opts.energy_window_meV(:)).');
margin = opts.initial_peak_margin_meV;
lower = min(max(center - half_split, win(1) + margin), win(2) - margin);
upper = min(max(center + half_split, win(1) + margin), win(2) - margin);
if upper - lower < opts.min_peak_separation_meV
    lower = max(win(1) + margin, center - opts.min_peak_separation_meV);
    upper = min(win(2) - margin, center + opts.min_peak_separation_meV);
end
guesses = sort([lower; upper]);

% Keep the variable referenced for MATLAB Code Analyzer clarity.
if isempty(q_axis)
    guesses = guesses(:);
end
end


function guesses = local_two_local_peak_guesses(energy_win, y, opts)
guesses = [];
if numel(energy_win) < 5 || all(~isfinite(y))
    return
end
y(~isfinite(y)) = 0;
try
    [pks, locs] = findpeaks(max(y, 0), energy_win, ...
        'MinPeakProminence', max(max(y, [], 'omitnan') * opts.min_prominence, 0), ...
        'MinPeakDistance', opts.min_peak_separation_meV, ...
        'SortStr', 'descend');
catch
    pks = [];
    locs = [];
end
if numel(locs) < 2
    return
end
[~, order] = sort(pks, 'descend');
locs = locs(order(1:2));
guesses = sort(locs(:));
end


function center = local_old_branch_center(old_points, q_value)
center = NaN;
if isempty(old_points) || ~istable(old_points) || ...
        ~all(ismember({'q_Ainv', 'energy_meV'}, old_points.Properties.VariableNames))
    return
end
q = double(old_points.q_Ainv);
E = double(old_points.energy_meV);
valid = isfinite(q) & isfinite(E);
if ~any(valid)
    return
end
[~, idx] = min(abs(q(valid) - q_value));
valid_idx = find(valid);
center = E(valid_idx(idx));
end


function [peak_energy, peak_ci] = local_peak_energy_and_ci(fit)
peak_energy = fit.omega_p(:);
peak_ci = local_result_ci(fit, 'omega_p_ci', numel(peak_energy));
if isfield(fit, 'apex_energy_meV')
    apex = fit.apex_energy_meV(:);
    use_apex = isfinite(apex);
    peak_energy(use_apex) = apex(use_apex);
    if isfield(fit, 'apex_energy_ci')
        apex_ci = fit.apex_energy_ci;
        for i = 1:numel(use_apex)
            if use_apex(i) && i <= size(apex_ci, 1) && size(apex_ci, 2) >= 2
                peak_ci(i, :) = apex_ci(i, 1:2);
            end
        end
    end
end
peak_ci = local_fill_energy_ci(peak_energy, peak_ci);
end


function ci = local_result_ci(fit, field_name, n_rows)
ci = NaN(n_rows, 2);
if isfield(fit, field_name)
    source = fit.(field_name);
    n = min(n_rows, size(source, 1));
    if size(source, 2) >= 2
        ci(1:n, :) = source(1:n, 1:2);
    end
end
end


function ci = local_fill_energy_ci(energy, ci)
for i = 1:numel(energy)
    center = energy(i);
    if ~isfinite(center)
        continue
    end
    half_width = max(1, 0.001 * abs(center));
    if ~all(isfinite(ci(i, :))) || ci(i, 1) >= center || ci(i, 2) <= center
        ci(i, :) = [center - half_width, center + half_width];
    else
        ci(i, 1) = min(ci(i, 1), center - half_width);
        ci(i, 2) = max(ci(i, 2), center + half_width);
    end
end
end


function row = local_point_row(unit, branch_id, branch_label, fit, ...
    peak_idx, peak_energy, peak_ci, energy_axis, raw_spectrum)
gamma = fit.gamma(peak_idx);
amplitude = fit.amplitude(peak_idx);
gamma_ci = local_ci_row(fit.gamma_ci, peak_idx, gamma);
amp_ci = local_ci_row(fit.amplitude_ci, peak_idx, amplitude);
raw_height = measure_peak_height(energy_axis, raw_spectrum, ...
    peak_energy, gamma);
E_ci_half = 0.5 * (peak_ci(2) - peak_ci(1));

row = table(unit.q_Ainv, unit.q_abs_Ainv, peak_energy, gamma, ...
    fit.R_squared, amplitude, peak_ci(1), peak_ci(2), ...
    gamma_ci(1), gamma_ci(2), amp_ci(1), amp_ci(2), raw_height, ...
    branch_id, {branch_label}, E_ci_half, {unit.source_mode}, ...
    unit.source_q_count, {unit.source_q_Ainv}, {unit.source_q_index}, ...
    unit.bin_size_requested, {fit.peak_model}, {fit.tracking_mode}, ...
    {fit.fit_spectrum_source}, ...
    {fit.fit_denoise_method}, fit.fit_denoise_window, ...
    fit.fit_denoise_order, {fit.fit_denoise_profile}, ...
    fit.fit_denoise_q_start_Ainv, fit.fit_denoise_q_end_Ainv, ...
    fit.fit_spectrum_delta_rms, {'raw_fit'}, peak_energy, {''}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'energy_meV', ...
    'gamma_meV', 'R2', 'amplitude_fit', 'E_ci_lo', 'E_ci_hi', ...
    'gamma_ci_lo', 'gamma_ci_hi', 'A_ci_lo', 'A_ci_hi', ...
    'raw_height', 'branch', 'branch_label', 'E_ci_half_meV', ...
    'source_mode', 'source_q_count', 'source_q_Ainv', ...
    'source_q_index', 'bin_size_requested', 'peak_model', ...
    'tracking_mode', 'fit_spectrum_source', 'fit_denoise_method', ...
    'fit_denoise_window', 'fit_denoise_order', 'fit_denoise_profile', ...
    'fit_denoise_q_start_Ainv', 'fit_denoise_q_end_Ainv', ...
    'fit_spectrum_delta_rms', 'repair_source', 'original_energy_meV', ...
    'repair_detail'});
end


function ci = local_ci_row(ci_matrix, idx, center)
ci = [NaN NaN];
if idx <= size(ci_matrix, 1) && size(ci_matrix, 2) >= 2
    ci = ci_matrix(idx, 1:2);
end
if ~all(isfinite(ci)) || ci(1) >= center || ci(2) <= center
    half_width = max(1, 0.001 * abs(center));
    ci = [center - half_width, center + half_width];
end
end


function row = local_failure_row(unit, status, detail, fit_meta)
row = table(unit.q_Ainv, unit.q_abs_Ainv, {unit.source_mode}, ...
    unit.source_q_count, {unit.source_q_Ainv}, {unit.source_q_index}, ...
    {status}, {detail}, {fit_meta.peak_model}, {fit_meta.tracking_mode}, ...
    {fit_meta.fit_spectrum_source}, ...
    {fit_meta.fit_denoise_method}, fit_meta.fit_denoise_window, ...
    fit_meta.fit_denoise_order, {fit_meta.fit_denoise_profile}, ...
    fit_meta.fit_denoise_q_start_Ainv, fit_meta.fit_denoise_q_end_Ainv, ...
    fit_meta.fit_spectrum_delta_rms, {'unrepaired_failure'}, {''}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'status', 'detail', 'peak_model', 'tracking_mode', ...
    'fit_spectrum_source', 'fit_denoise_method', 'fit_denoise_window', ...
    'fit_denoise_order', 'fit_denoise_profile', ...
    'fit_denoise_q_start_Ainv', 'fit_denoise_q_end_Ainv', ...
    'fit_spectrum_delta_rms', 'repair_source', 'repair_detail'});
end


function row = local_binning_map_row(unit, guesses, success, status, fit_meta)
row = table(unit.q_Ainv, unit.q_abs_Ainv, {unit.source_mode}, ...
    unit.source_q_count, {unit.source_q_Ainv}, {unit.source_q_index}, ...
    unit.bin_size_requested, guesses(1), guesses(2), success, {status}, ...
    {fit_meta.peak_model}, {fit_meta.tracking_mode}, ...
    {fit_meta.fit_spectrum_source}, {fit_meta.fit_denoise_method}, ...
    fit_meta.fit_denoise_window, fit_meta.fit_denoise_order, ...
    {fit_meta.fit_denoise_profile}, fit_meta.fit_denoise_q_start_Ainv, ...
    fit_meta.fit_denoise_q_end_Ainv, fit_meta.fit_spectrum_delta_rms, ...
    {'raw_fit'}, {''}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'bin_size_requested', 'initial_guess_lower_meV', ...
    'initial_guess_upper_meV', 'fit_success', 'fit_status', ...
    'peak_model', 'tracking_mode', 'fit_spectrum_source', ...
    'fit_denoise_method', 'fit_denoise_window', ...
    'fit_denoise_order', 'fit_denoise_profile', ...
    'fit_denoise_q_start_Ainv', 'fit_denoise_q_end_Ainv', ...
    'fit_spectrum_delta_rms', 'repair_source', 'repair_detail'});
end


function tbl = local_empty_points_table()
names = {'q_Ainv', 'q_abs_Ainv', 'energy_meV', ...
    'gamma_meV', 'R2', 'amplitude_fit', 'E_ci_lo', 'E_ci_hi', ...
    'gamma_ci_lo', 'gamma_ci_hi', 'A_ci_lo', 'A_ci_hi', ...
    'raw_height', 'branch', 'branch_label', 'E_ci_half_meV', ...
    'source_mode', 'source_q_count', 'source_q_Ainv', ...
    'source_q_index', 'bin_size_requested', 'peak_model', ...
    'tracking_mode', 'fit_spectrum_source', 'fit_denoise_method', ...
    'fit_denoise_window', 'fit_denoise_order', 'fit_denoise_profile', ...
    'fit_denoise_q_start_Ainv', 'fit_denoise_q_end_Ainv', ...
    'fit_spectrum_delta_rms', 'repair_source', 'original_energy_meV', ...
    'repair_detail'};
types = {'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'double', 'double', 'double', 'double', ...
    'double', 'double', 'cell', 'double', 'cell', 'double', 'cell', ...
    'cell', 'double', 'cell', 'cell', 'cell', 'cell', 'double', ...
    'double', 'cell', 'double', 'double', 'double', 'cell', 'double', ...
    'cell'};
tbl = table('Size', [0 numel(names)], 'VariableTypes', types, ...
    'VariableNames', names);
end


function tbl = local_empty_failures_table()
names = {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'status', 'detail', 'peak_model', 'tracking_mode', ...
    'fit_spectrum_source', 'fit_denoise_method', 'fit_denoise_window', ...
    'fit_denoise_order', 'fit_denoise_profile', ...
    'fit_denoise_q_start_Ainv', 'fit_denoise_q_end_Ainv', ...
    'fit_spectrum_delta_rms', 'repair_source', 'repair_detail'};
types = {'double', 'double', 'cell', 'double', 'cell', 'cell', ...
    'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'double', ...
    'double', 'cell', 'double', 'double', 'double', 'cell', 'cell'};
tbl = table('Size', [0 numel(names)], 'VariableTypes', types, ...
    'VariableNames', names);
end


function tbl = local_empty_binning_map()
names = {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'bin_size_requested', 'initial_guess_lower_meV', ...
    'initial_guess_upper_meV', 'fit_success', 'fit_status', ...
    'peak_model', 'tracking_mode', 'fit_spectrum_source', ...
    'fit_denoise_method', 'fit_denoise_window', ...
    'fit_denoise_order', 'fit_denoise_profile', ...
    'fit_denoise_q_start_Ainv', 'fit_denoise_q_end_Ainv', ...
    'fit_spectrum_delta_rms', 'repair_source', 'repair_detail'};
types = {'double', 'double', 'cell', 'double', 'cell', 'cell', ...
    'double', 'double', 'double', 'logical', 'cell', 'cell', 'cell', ...
    'cell', 'cell', 'double', 'double', 'cell', 'double', 'double', ...
    'double', 'cell', 'cell'};
tbl = table('Size', [0 numel(names)], 'VariableTypes', types, ...
    'VariableNames', names);
end


function tbl = local_empty_repair_log_table()
names = {'q_Ainv', 'q_abs_Ainv', 'branch', 'branch_label', ...
    'repair_source', 'original_energy_meV', 'repaired_energy_meV', ...
    'source_mode', 'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'detail'};
types = {'double', 'double', 'double', 'cell', 'cell', 'double', ...
    'double', 'cell', 'double', 'cell', 'cell', 'cell'};
tbl = table('Size', [0 numel(names)], 'VariableTypes', types, ...
    'VariableNames', names);
end


function row = local_repair_log_row(unit, branch_id, branch_label, ...
    source, original_energy, repaired_energy, detail)
row = table(unit.q_Ainv, unit.q_abs_Ainv, branch_id, {branch_label}, ...
    {source}, original_energy, repaired_energy, {unit.source_mode}, ...
    unit.source_q_count, {unit.source_q_Ainv}, {unit.source_q_index}, ...
    {detail}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'branch', ...
    'branch_label', 'repair_source', 'original_energy_meV', ...
    'repaired_energy_meV', 'source_mode', 'source_q_count', ...
    'source_q_Ainv', 'source_q_index', 'detail'});
end


function text = local_join_numbers(values)
parts = arrayfun(@(x) sprintf('%.12g', x), values(:).', ...
    'UniformOutput', false);
text = strjoin(parts, '|');
end
