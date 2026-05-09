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

lower_points = local_empty_points_table();
upper_points = local_empty_points_table();
combined_points = local_empty_points_table();
fit_failures = local_empty_failures_table();
binning_map = local_empty_binning_map();
fit_details = cell(numel(units), 1);

for ui = 1:numel(units)
    unit = units(ui);
    spectrum = mean(intensity(:, unit.q_indices), 2, 'omitnan');
    raw_spectrum = mean(raw_intensity(:, unit.q_indices), 2, 'omitnan');
    guesses = local_initial_guesses(energy_axis, spectrum, q_axis, ...
        unit, energy_mask, opts);

    try
        fit = fit_loss_function(energy_axis, spectrum, ...
            'E_min', energy_window(1), ...
            'E_max', energy_window(2), ...
            'max_peaks', 2, ...
            'min_prominence', opts.min_prominence, ...
            'smooth_width', opts.smooth_width, ...
            'initial_guesses', guesses(:), ...
            'peak_model', opts.peak_model, ...
            'pre_subtracted', opts.pre_subtracted, ...
            'bootstrap_ci_samples', opts.bootstrap_ci_samples);
    catch ME
        fit_failures = [fit_failures; local_failure_row(unit, ...
            'double_peak_fit_failed', ME.message)]; %#ok<AGROW>
        binning_map = [binning_map; local_binning_map_row(unit, ...
            guesses, false, 'double_peak_fit_failed')]; %#ok<AGROW>
        continue
    end

    if fit.n_peaks ~= 2
        fit_failures = [fit_failures; local_failure_row(unit, ...
            'not_two_peaks', sprintf('fit returned %d peaks', fit.n_peaks))]; %#ok<AGROW>
        binning_map = [binning_map; local_binning_map_row(unit, ...
            guesses, false, 'not_two_peaks')]; %#ok<AGROW>
        continue
    end

    [peak_energy, peak_ci] = local_peak_energy_and_ci(fit);
    [peak_energy_sorted, sort_idx] = sort(peak_energy(:));
    if numel(peak_energy_sorted) ~= 2 || any(~isfinite(peak_energy_sorted))
        fit_failures = [fit_failures; local_failure_row(unit, ...
            'nonfinite_double_peak_energy', 'one or both peak energies are non-finite')]; %#ok<AGROW>
        binning_map = [binning_map; local_binning_map_row(unit, ...
            guesses, false, 'nonfinite_double_peak_energy')]; %#ok<AGROW>
        continue
    end
    if diff(peak_energy_sorted) < opts.min_peak_separation_meV
        fit_failures = [fit_failures; local_failure_row(unit, ...
            'collapsed_double_peak', sprintf('peak separation %.3g meV', ...
            diff(peak_energy_sorted)))]; %#ok<AGROW>
        binning_map = [binning_map; local_binning_map_row(unit, ...
            guesses, false, 'collapsed_double_peak')]; %#ok<AGROW>
        continue
    end

    lower_idx = sort_idx(1);
    upper_idx = sort_idx(2);
    lower_row = local_point_row(unit, 1, 'b1_double_peak_lower', ...
        fit, lower_idx, peak_energy(lower_idx), peak_ci(lower_idx, :), ...
        energy_axis, raw_spectrum);
    upper_row = local_point_row(unit, 2, 'b1_double_peak_upper', ...
        fit, upper_idx, peak_energy(upper_idx), peak_ci(upper_idx, :), ...
        energy_axis, raw_spectrum);

    lower_points = [lower_points; lower_row]; %#ok<AGROW>
    upper_points = [upper_points; upper_row]; %#ok<AGROW>
    combined_points = [combined_points; lower_row; upper_row]; %#ok<AGROW>
    fit_details{ui} = fit;
    binning_map = [binning_map; local_binning_map_row(unit, ...
        guesses, true, 'ok')]; %#ok<AGROW>
end

result = struct();
result.lower_points = sortrows(lower_points, {'q_Ainv', 'energy_meV'});
result.upper_points = sortrows(upper_points, {'q_Ainv', 'energy_meV'});
result.combined_points = sortrows(combined_points, ...
    {'q_Ainv', 'branch', 'energy_meV'});
result.binning_map = binning_map;
result.noise_profile = noise_profile;
result.fit_failures = fit_failures;
result.fit_details = fit_details;
end


function opts = local_defaults(opts)
opts = local_set_default(opts, 'energy_window_meV', [500 2100]);
opts = local_set_default(opts, 'q_range_Ainv', [-0.15 0.15]);
opts = local_set_default(opts, 'q_skip_Ainv', 0.005);
opts = local_set_default(opts, 'bin_size', 3);
opts = local_set_default(opts, 'noise_threshold', NaN);
opts = local_set_default(opts, 'low_q_no_bin_abs_Ainv', 0.03);
opts = local_set_default(opts, 'edge_fraction', 0.18);
opts = local_set_default(opts, 'peak_model', 'fano');
opts = local_set_default(opts, 'pre_subtracted', false);
opts = local_set_default(opts, 'min_prominence', 0.02);
opts = local_set_default(opts, 'smooth_width', 1);
opts = local_set_default(opts, 'bootstrap_ci_samples', 0);
opts = local_set_default(opts, 'old_branch_points', table());
opts = local_set_default(opts, 'fallback_split_meV', 180);
opts = local_set_default(opts, 'min_peak_separation_meV', 10);
opts = local_set_default(opts, 'initial_peak_margin_meV', 25);
opts.bin_size = max(1, round(double(opts.bin_size)));
opts.smooth_width = max(1, round(double(opts.smooth_width)));
opts.low_q_no_bin_abs_Ainv = max(0, double(opts.low_q_no_bin_abs_Ainv));
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
[use_binning, noise_exceeds_threshold, low_q_protected] = ...
    local_threshold_binning_mask(q, noise_ratio, threshold, ...
    opts.low_q_no_bin_abs_Ainv);

profile = table(q_index, q, abs(q), signal, noise, noise_ratio, ...
    repmat(threshold, n, 1), ...
    repmat(opts.low_q_no_bin_abs_Ainv, n, 1), ...
    repmat({threshold_source}, n, 1), ...
    noise_exceeds_threshold, low_q_protected, use_binning, ...
    'VariableNames', {'q_index', 'q_Ainv', 'q_abs_Ainv', ...
    'b1_signal', 'b1_noise_mad', 'noise_ratio', ...
    'noise_threshold', 'low_q_no_bin_abs_Ainv', ...
    'noise_threshold_source', 'noise_exceeds_threshold', ...
    'low_q_no_bin_protected', 'use_binning'});
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


function [use_binning, noise_exceeds_threshold, low_q_protected] = ...
    local_threshold_binning_mask(q, noise_ratio, threshold, low_q_no_bin_abs)
use_binning = false(size(q));
noise_exceeds_threshold = false(size(q));
low_q_protected = false(size(q));
if ~isfinite(threshold)
    return
end
noise_exceeds_threshold = isfinite(noise_ratio) & noise_ratio > threshold;
if isfinite(low_q_no_bin_abs) && low_q_no_bin_abs > 0
    low_q_protected = isfinite(q) & abs(q) <= low_q_no_bin_abs;
end
use_binning = noise_exceeds_threshold & ~low_q_protected;
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

half_split = opts.fallback_split_meV / 2;
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
    unit.bin_size_requested, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'energy_meV', ...
    'gamma_meV', 'R2', 'amplitude_fit', 'E_ci_lo', 'E_ci_hi', ...
    'gamma_ci_lo', 'gamma_ci_hi', 'A_ci_lo', 'A_ci_hi', ...
    'raw_height', 'branch', 'branch_label', 'E_ci_half_meV', ...
    'source_mode', 'source_q_count', 'source_q_Ainv', ...
    'source_q_index', 'bin_size_requested'});
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


function row = local_failure_row(unit, status, detail)
row = table(unit.q_Ainv, unit.q_abs_Ainv, {unit.source_mode}, ...
    unit.source_q_count, {unit.source_q_Ainv}, {unit.source_q_index}, ...
    {status}, {detail}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'status', 'detail'});
end


function row = local_binning_map_row(unit, guesses, success, status)
row = table(unit.q_Ainv, unit.q_abs_Ainv, {unit.source_mode}, ...
    unit.source_q_count, {unit.source_q_Ainv}, {unit.source_q_index}, ...
    unit.bin_size_requested, guesses(1), guesses(2), success, {status}, ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'bin_size_requested', 'initial_guess_lower_meV', ...
    'initial_guess_upper_meV', 'fit_success', 'fit_status'});
end


function tbl = local_empty_points_table()
tbl = table(zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    cell(0,1), zeros(0,1), cell(0,1), zeros(0,1), cell(0,1), ...
    cell(0,1), zeros(0,1), ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'energy_meV', ...
    'gamma_meV', 'R2', 'amplitude_fit', 'E_ci_lo', 'E_ci_hi', ...
    'gamma_ci_lo', 'gamma_ci_hi', 'A_ci_lo', 'A_ci_hi', ...
    'raw_height', 'branch', 'branch_label', 'E_ci_half_meV', ...
    'source_mode', 'source_q_count', 'source_q_Ainv', ...
    'source_q_index', 'bin_size_requested'});
end


function tbl = local_empty_failures_table()
tbl = table(zeros(0,1), zeros(0,1), cell(0,1), zeros(0,1), ...
    cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'status', 'detail'});
end


function tbl = local_empty_binning_map()
tbl = table(zeros(0,1), zeros(0,1), cell(0,1), zeros(0,1), ...
    cell(0,1), cell(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), false(0,1), cell(0,1), ...
    'VariableNames', {'q_Ainv', 'q_abs_Ainv', 'source_mode', ...
    'source_q_count', 'source_q_Ainv', 'source_q_index', ...
    'bin_size_requested', 'initial_guess_lower_meV', ...
    'initial_guess_upper_meV', 'fit_success', 'fit_status'});
end


function text = local_join_numbers(values)
parts = arrayfun(@(x) sprintf('%.12g', x), values(:).', ...
    'UniformOutput', false);
text = strjoin(parts, '|');
end
