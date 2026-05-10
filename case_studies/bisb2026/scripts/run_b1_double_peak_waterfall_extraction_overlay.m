function out = run_b1_double_peak_waterfall_extraction_overlay(options)
%RUN_B1_DOUBLE_PEAK_WATERFALL_EXTRACTION_OVERLAY Overlay B1 double peaks.
%
% This diagnostic uses the saved spectra/extraction results for one
% b1_double_peak_binning run. It does not refit; it places the already
% extracted lower/upper peak centers on the same signed-q waterfall traces.

arguments
    options.dateTag {mustBeTextScalar} = "260509_highqbin5_strong_adaptive_sgdenoise"
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
results_root = fullfile(project_root, 'paper_results');
date_tag = char(string(options.dateTag));

sessions = local_sessions(date_tag);
panel_paths = strings(1, numel(sessions));
session_outputs = cell(1, numel(sessions));
for i = 1:numel(sessions)
    folder = fullfile(results_root, sessions(i).folder);
    session_outputs{i} = local_overlay_one_session(folder);
    panel_paths(i) = string(session_outputs{i}.png);
end

comparison_dir = fullfile(results_root, ...
    sprintf('b1_double_peak_waterfall_extraction_overlay_three_dataset_comparison_%s', date_tag));
if ~isfolder(comparison_dir)
    mkdir(comparison_dir);
end

images = cell(1, numel(panel_paths));
for i = 1:numel(panel_paths)
    images{i} = local_read_rgb(panel_paths(i));
end
canvas = local_concat_horizontal(images, 55);

out = struct();
out.sessions = session_outputs;
out.output_dir = comparison_dir;
out.png = fullfile(comparison_dir, ...
    'b1_double_peak_waterfall_extraction_overlay_three_dataset_comparison.png');
out.pdf = fullfile(comparison_dir, ...
    'b1_double_peak_waterfall_extraction_overlay_three_dataset_comparison.pdf');
if strlength(string(out.pdf)) > 235
    out.pdf = fullfile(comparison_dir, ...
        'b1_double_peak_overlay_three_dataset_comparison.pdf');
end
imwrite(canvas, out.png);

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [50 50 min(size(canvas, 2), 5600) min(size(canvas, 1), 2800)]);
ax = axes(fig);
image(ax, canvas);
axis(ax, 'image');
axis(ax, 'off');
title(ax, sprintf('B1 double-peak extraction overlay | %s', date_tag), ...
    'Interpreter', 'none');
exportgraphics(fig, out.pdf, 'ContentType', 'vector');
if ~isfile(out.pdf)
    print(fig, out.pdf, '-dpdf', '-bestfit');
end
close(fig);

fprintf('B1 double-peak waterfall extraction overlay written:\n  %s\n', out.png);
end


function sessions = local_sessions(date_tag)
sessions = struct( ...
    'folder', { ...
        sprintf('590_gui_history_area_260506_b1_double_peak_binning_%s', date_tag), ...
        sprintf('n0_PL2_10w_gui_history_area_260506_b1_double_peak_binning_%s', date_tag), ...
        sprintf('no_PL2_20w_2film_gui_history_area_260506_highq_refined_b1_double_peak_binning_%s', date_tag)});
end


function out = local_overlay_one_session(folder)
mat_path = fullfile(folder, 'b1_double_peak_binning_results.mat');
if ~isfile(mat_path)
    error('run_b1_double_peak_waterfall_extraction_overlay:MissingMat', ...
        'Missing extraction MAT file: %s', mat_path);
end
saved = load(mat_path, 'extract', 'extract_opts', 'input_dir', 'session');
analysis = load(fullfile(saved.input_dir, 'analysis_results.mat'), 'output');
qe = analysis.output.qe_pp;
extract = saved.extract;
opts = saved.extract_opts;
session = saved.session;

[energy_axis, energy_mask] = local_waterfall_energy_window(qe, opts);
[q_values, traces, groups] = local_waterfall_trace_set(qe, extract, ...
    energy_mask);
[norm_energy_axis, norm_energy_mask] = local_waterfall_norm_window(qe, opts);
[~, norm_traces] = local_waterfall_trace_set(qe, extract, norm_energy_mask);
normalized = local_visual_normalize_traces(energy_axis, traces, ...
    norm_energy_axis, norm_traces, opts);

offset = 0.55;
offsets = (0:(size(normalized, 2) - 1)) .* offset;

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 40 760 1120]);
ax = axes(fig);
hold(ax, 'on');
for i = 1:size(normalized, 2)
    plot(ax, energy_axis, normalized(:, i) + offsets(i), '-', ...
        'Color', local_waterfall_color(i, size(normalized, 2)), ...
        'LineWidth', 1.05);
end

overlay_rows = local_overlay_branch(ax, extract.lower_points, energy_axis, ...
    normalized, offsets, q_values, groups, [0.00 0.40 1.00], 'lower');
overlay_rows = [overlay_rows; local_overlay_branch(ax, extract.upper_points, ...
    energy_axis, normalized, offsets, q_values, groups, [1.00 0.10 0.60], ...
    'upper')]; %#ok<AGROW>

hold(ax, 'off');
xlim(ax, [max(0, min(energy_axis)), min(max(energy_axis), ...
    max(opts.energy_window_meV(2), 2000))]);
ylim(ax, [-0.2, offsets(end) + 1.45]);
ax.YTick = [];
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.LineWidth = 1.35;
ax.FontSize = 18;
ax.TickDir = 'out';
ax.XTick = [0 1000 2000];
xlabel(ax, 'Energy loss (meV)', 'FontSize', 24);
title(ax, {'B1 double-peak waterfall extraction overlay', ...
    char(session.session_label), ...
    local_overlay_title_line(opts)}, ...
    'FontSize', 12, 'Interpreter', 'none');
box(ax, 'off');

out = struct();
out.output_dir = folder;
out.png = fullfile(folder, ...
    'b1_double_peak_waterfall_extraction_overlay_signed_q.png');
out.pdf = fullfile(folder, ...
    'b1_double_peak_waterfall_extraction_overlay_signed_q.pdf');
out.csv = fullfile(folder, ...
    'b1_double_peak_waterfall_extraction_overlay_points.csv');

exportgraphics(fig, out.png, 'Resolution', 300);
exportgraphics(fig, out.pdf, 'ContentType', 'vector');
close(fig);

if isempty(overlay_rows)
    overlay_table = table();
else
    overlay_table = cell2table(overlay_rows, 'VariableNames', { ...
        'branch', 'q_Ainv', 'q_abs_Ainv', 'energy_meV', ...
        'plot_trace_index', 'plot_q_Ainv', 'plot_y', ...
        'source_mode', 'source_q_count', 'source_q_Ainv'});
end
writetable(overlay_table, out.csv);
end


function line = local_overlay_title_line(opts)
parts = {'signed q'};
if isfield(opts, 'high_q_force_bin_abs_Ainv') && ...
        isfinite(opts.high_q_force_bin_abs_Ainv) && ...
        isfield(opts, 'bin_size')
    parts{end + 1} = sprintf('high-q bin%d >= %.2f', ...
        opts.bin_size, opts.high_q_force_bin_abs_Ainv); %#ok<AGROW>
end
if isfield(opts, 'fit_denoise_method') && ...
        strcmpi(char(opts.fit_denoise_method), 'sgolay')
    if isfield(opts, 'fit_denoise_profile') && ...
            contains(char(opts.fit_denoise_profile), 'adaptive')
        parts{end + 1} = sprintf('adaptive SG %d-%d', ...
            opts.fit_denoise_low_window, ...
            opts.fit_denoise_high_window); %#ok<AGROW>
    else
        parts{end + 1} = sprintf('SG %d', ...
            opts.fit_denoise_window); %#ok<AGROW>
    end
end
parts{end + 1} = 'lower blue'; %#ok<AGROW>
parts{end + 1} = 'upper magenta'; %#ok<AGROW>
line = strjoin(parts, ' | ');
end


function rows = local_overlay_branch(ax, points, energy_axis, normalized, ...
    offsets, q_values, groups, color, branch_label)
rows = {};
if isempty(points) || height(points) == 0
    return
end
for i = 1:height(points)
    q = double(points.q_Ainv(i));
    energy = double(points.energy_meV(i));
    [~, trace_idx] = min(abs(q_values - q));
    if isempty(trace_idx) || ~isfinite(energy) || energy < min(energy_axis) || ...
            energy > max(energy_axis)
        continue
    end
    y = interp1(energy_axis, normalized(:, trace_idx), energy, 'linear', NaN) + ...
        offsets(trace_idx);
    if ~isfinite(y)
        continue
    end
    marker = 'o';
    if startsWith(local_table_text(points.source_mode, i), 'combined_q_binning')
        marker = 's';
    end
    plot(ax, energy, y, marker, 'MarkerSize', 5.8, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.75);

    rows(end + 1, :) = {branch_label, q, double(points.q_abs_Ainv(i)), ...
        energy, trace_idx, q_values(trace_idx), y, ...
        local_table_text(points.source_mode, i), ...
        double(points.source_q_count(i)), ...
        local_table_text(points.source_q_Ainv, i)}; %#ok<AGROW>
end
end


function [energy_axis, energy_mask] = local_waterfall_energy_window(qe, opts)
full_energy = double(qe.energy_meV(:));
e_min = max(opts.waterfall_start_meV, max(0, min(full_energy)));
e_max = min(max(full_energy), max(opts.energy_window_meV(2), 2000));
if isfield(opts, 'waterfall_end_meV') && isfinite(opts.waterfall_end_meV)
    e_max = min(max(full_energy), opts.waterfall_end_meV);
end
energy_mask = full_energy >= e_min & full_energy <= e_max;
energy_axis = full_energy(energy_mask);
end


function [energy_axis, energy_mask] = local_waterfall_norm_window(qe, opts)
full_energy = double(qe.energy_meV(:));
window = opts.waterfall_area_norm_window_meV;
energy_mask = full_energy >= window(1) & full_energy <= window(2);
if ~any(energy_mask)
    energy_mask = full_energy >= 250;
end
if ~any(energy_mask)
    energy_mask = true(size(full_energy));
end
energy_axis = full_energy(energy_mask);
end


function [q_values, traces, groups] = local_waterfall_trace_set(qe, extract, ...
    energy_mask)
q_axis = double(qe.q_Ainv(:));
groups = local_binning_map_q_groups(extract, q_axis);
traces = zeros(nnz(energy_mask), 0);
q_values = zeros(0, 1);
for i = 1:numel(groups)
    q_idx = groups(i).indices;
    q_idx = q_idx(q_idx >= 1 & q_idx <= numel(q_axis));
    if isempty(q_idx)
        continue
    end
    traces(:, end + 1) = mean(double(qe.intensity(energy_mask, q_idx)), ...
        2, 'omitnan'); %#ok<AGROW>
    q_values(end + 1, 1) = mean(q_axis(q_idx), 'omitnan'); %#ok<AGROW>
end
end


function groups = local_binning_map_q_groups(extract, q_axis)
groups = struct('indices', {}, 'source_mode', {}, 'q_mean', {});
map = extract.binning_map;
for i = 1:height(map)
    q_idx = local_parse_index_list(local_table_text(map.source_q_index, i));
    q_idx = q_idx(q_idx >= 1 & q_idx <= numel(q_axis));
    if isempty(q_idx)
        continue
    end
    item.indices = q_idx;
    item.source_mode = local_table_text(map.source_mode, i);
    item.q_mean = mean(q_axis(q_idx), 'omitnan');
    groups(end + 1) = item; %#ok<AGROW>
end
if isempty(groups)
    return
end
[~, order] = sort([groups.q_mean], 'ascend');
groups = groups(order);
end


function normalized = local_visual_normalize_traces(energy_axis, traces, ...
    norm_energy_axis, norm_traces, opts)
normalized = local_area_normalize_traces(traces, norm_energy_axis, norm_traces);
normalized = local_apply_waterfall_residual_and_gain(energy_axis, ...
    normalized, opts);
end


function normalized = local_area_normalize_traces(traces, norm_energy_axis, ...
    norm_traces)
normalized = zeros(size(traces));
for i = 1:size(traces, 2)
    y_norm = double(norm_traces(:, i));
    valid = isfinite(norm_energy_axis) & isfinite(y_norm);
    area = NaN;
    if nnz(valid) >= 2
        area = trapz(norm_energy_axis(valid), y_norm(valid));
    end
    if ~isfinite(area) || abs(area) <= eps
        area = 1;
    end
    normalized(:, i) = double(traces(:, i)) ./ area;
end
scale_values = abs(normalized(isfinite(normalized)));
if isempty(scale_values)
    return
end
display_scale = prctile(scale_values, 95);
if isfinite(display_scale) && display_scale > eps
    normalized = normalized .* (0.85 ./ display_scale);
end
end


function traces = local_apply_waterfall_residual_and_gain(energy_axis, ...
    traces, opts)
if isfield(opts, 'waterfall_residual') && opts.waterfall_residual
    for i = 1:size(traces, 2)
        y = double(traces(:, i));
        baseline = local_waterfall_asls_baseline(y, 1e6, 0.01, 12);
        traces(:, i) = y - baseline;
    end
    scale_values = abs(traces(isfinite(traces)));
    if ~isempty(scale_values)
        scale = prctile(scale_values, 95);
        if isfinite(scale) && scale > eps
            traces = traces .* (0.85 ./ scale);
        end
    end
end
if isfield(opts, 'waterfall_gain') && isfinite(opts.waterfall_gain) && ...
        opts.waterfall_gain > 0
    traces = traces .* opts.waterfall_gain;
end
end


function baseline = local_waterfall_asls_baseline(y, lambda_value, ...
    asymmetry, iterations)
y = double(y(:));
n = numel(y);
if n < 3
    baseline = y;
    return
end
lambda_value = max(double(lambda_value), 1);
asymmetry = min(max(double(asymmetry), 1e-4), 0.49);
d = diff(speye(n), 2);
w = ones(n, 1);
for i = 1:max(1, round(iterations))
    w_mat = spdiags(w, 0, n, n);
    baseline = (w_mat + lambda_value * (d' * d)) \ (w .* y);
    w = asymmetry * (y > baseline) + (1 - asymmetry) * (y <= baseline);
end
end


function color = local_waterfall_color(index, n_traces)
if n_traces <= 1
    t = 0.5;
else
    t = (index - 1) ./ (n_traces - 1);
end
if t < 0.5
    u = t ./ 0.5;
    color = (1 - u) .* [1.0 0.0 0.0] + u .* [0.05 0.05 0.05];
else
    u = (t - 0.5) ./ 0.5;
    color = (1 - u) .* [0.05 0.05 0.05] + u .* [0.0 0.80 0.10];
end
end


function values = local_parse_index_list(text)
parts = strsplit(char(text), '|');
values = str2double(parts);
values = values(isfinite(values));
values = max(1, round(values(:).'));
end


function text = local_table_text(value, row)
if iscell(value)
    text = char(value{row});
elseif isstring(value)
    text = char(value(row));
else
    text = char(string(value(row)));
end
end


function img = local_read_rgb(path)
img = imread(path);
if ndims(img) == 2
    img = repmat(img, 1, 1, 3);
elseif size(img, 3) > 3
    img = img(:, :, 1:3);
end
end


function canvas = local_concat_horizontal(images, gap_px)
n = numel(images);
heights = cellfun(@(x) size(x, 1), images);
widths = cellfun(@(x) size(x, 2), images);
max_h = max(heights);
total_w = sum(widths) + gap_px * (n - 1);
canvas = uint8(255 * ones(max_h, total_w, 3));
x0 = 1;
for i = 1:n
    img = images{i};
    h = size(img, 1);
    w = size(img, 2);
    y0 = floor((max_h - h) / 2) + 1;
    canvas(y0:(y0+h-1), x0:(x0+w-1), :) = img;
    x0 = x0 + w + gap_px;
end
end
