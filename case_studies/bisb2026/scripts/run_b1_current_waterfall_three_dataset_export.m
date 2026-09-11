function out = run_b1_current_waterfall_three_dataset_export(options)
%RUN_B1_CURRENT_WATERFALL_THREE_DATASET_EXPORT Export current-qcalib B1 waterfall.
%
% Thesis-facing plot. This reproduces the readable signed-q waterfall style
% from the previous diagnostic figure, but uses only the current GUI-history
% outputs and overlays only the accepted B1 apex points.

arguments
    options.outputTag {mustBeTextScalar} = "260519"
    options.maxTracesPerPanel (1, 1) double = 85
    options.energyWindowMeV (1, 2) double = [250 1700]
    options.areaNormWindowMeV (1, 2) double = [50 3800]
    options.waterfallGain (1, 1) double = 2.0
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
results_root = fullfile(project_root, 'paper_results');
out_dir = fullfile(results_root, ...
    sprintf('b1_current_waterfall_three_dataset_%s', char(options.outputTag)));
if ~isfolder(out_dir)
    mkdir(out_dir);
end

sessions = local_sessions();
panel_paths = strings(1, numel(sessions));
panel_results = cell(1, numel(sessions));
summary_rows = {};
for i = 1:numel(sessions)
    [panel_results{i}, rows] = local_plot_session(results_root, ...
        out_dir, sessions(i), options);
    panel_paths(i) = string(panel_results{i}.png);
    summary_rows = [summary_rows; rows]; %#ok<AGROW>
end

images = cell(1, numel(panel_paths));
for i = 1:numel(panel_paths)
    images{i} = local_read_rgb(panel_paths(i));
end
canvas = local_concat_horizontal(images, 55);

out = struct();
out.output_dir = out_dir;
out.png = fullfile(out_dir, 'b1_current_waterfall_three_dataset.png');
out.pdf = fullfile(out_dir, 'b1_current_waterfall_three_dataset.pdf');
out.csv = fullfile(out_dir, 'b1_current_waterfall_trace_summary.csv');
out.report = fullfile(out_dir, 'b1_current_waterfall_report.md');
out.panels = panel_results;

imwrite(canvas, out.png);
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [50 50 min(size(canvas, 2), 5600) min(size(canvas, 1), 2800)]);
ax = axes(fig);
image(ax, canvas);
axis(ax, 'image');
axis(ax, 'off');
exportgraphics(fig, out.pdf, 'ContentType', 'vector');
if ~isfile(out.pdf)
    print(fig, out.pdf, '-dpdf', '-bestfit');
end
close(fig);

if isempty(summary_rows)
    summary = table();
else
    summary = cell2table(summary_rows, 'VariableNames', { ...
        'dataset_label', 'session_folder', 'q_Ainv', 'source_q_Ainv', ...
        'q_column', 'apex_energy_meV', 'plot_trace_index', 'plot_y'});
end
writetable(summary, out.csv);
local_write_report(out, sessions, panel_results, options);

fprintf('B1 current waterfall export written:\n');
fprintf('  %s\n', out.png);
fprintf('  %s\n', out.pdf);
fprintf('  %s\n', out.csv);
fprintf('  %s\n', out.report);
end


function sessions = local_sessions()
sessions = struct( ...
    'label', { ...
        '590 10w defocus 1film', ...
        'n0 10w defocus repeat 1film', ...
        '20w defocus 2film'}, ...
    'folder', { ...
        '590_gui_history_area_260506', ...
        'n0_PL2_10w_gui_history_area_260506', ...
        'no_PL2_20w_2film_gui_history_area_260506_highq_refined'});
end


function [panel, rows] = local_plot_session(results_root, out_dir, session, options)
folder = fullfile(results_root, session.folder);
mat_path = fullfile(folder, 'analysis_results.mat');
points_path = fullfile(folder, 'branch1_points.csv');
if ~isfile(mat_path)
    error('run_b1_current_waterfall_three_dataset_export:MissingMat', ...
        'Missing analysis MAT: %s', mat_path);
end
if ~isfile(points_path)
    error('run_b1_current_waterfall_three_dataset_export:MissingCsv', ...
        'Missing B1 points CSV: %s', points_path);
end

saved = load(mat_path, 'output');
qe = saved.output.qe_pp;
points = readtable(points_path);
points = local_filter_points(points, qe, options.energyWindowMeV);

[energy_axis, energy_mask] = local_energy_window(qe, options.energyWindowMeV);
[norm_energy_axis, norm_energy_mask] = local_energy_window(qe, ...
    options.areaNormWindowMeV);
[q_values, q_columns, traces] = local_trace_set(qe, points, energy_mask, ...
    options.maxTracesPerPanel);
[~, ~, norm_traces] = local_trace_set(qe, points, norm_energy_mask, ...
    options.maxTracesPerPanel);
normalized = local_visual_normalize_traces(energy_axis, traces, ...
    norm_energy_axis, norm_traces, options);

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

rows = local_overlay_b1_points(ax, points, energy_axis, normalized, ...
    offsets, q_values, q_columns, session);

hold(ax, 'off');
xlim(ax, [max(0, min(energy_axis)), min(max(energy_axis), 2000)]);
if isempty(offsets)
    ylim(ax, [-0.2 1.0]);
else
    ylim(ax, [-0.2, offsets(end) + 1.45]);
end
ax.YTick = [];
ax.YColor = [0 0 0];
ax.XColor = [0 0 0];
ax.LineWidth = 1.35;
ax.FontSize = 18;
ax.TickDir = 'out';
ax.XTick = [0 1000 2000];
xlabel(ax, 'Energy loss (meV)', 'FontSize', 24);
title(ax, {'B1 waterfall extraction overlay', ...
    char(session.label), ...
    'signed q | current q-calib | B1 apex only'}, ...
    'FontSize', 12, 'Interpreter', 'none');
box(ax, 'off');

safe_name = regexprep(char(session.folder), '[^A-Za-z0-9_]+', '_');
panel = struct();
panel.label = session.label;
panel.folder = session.folder;
panel.png = fullfile(out_dir, sprintf('%s_b1_current_waterfall_overlay.png', ...
    safe_name));
panel.pdf = fullfile(out_dir, sprintf('%s_b1_current_waterfall_overlay.pdf', ...
    safe_name));
panel.n_total_b1_points = height(points);
panel.n_plotted_traces = size(normalized, 2);
panel.q_min_Ainv = min(q_values);
panel.q_max_Ainv = max(q_values);
panel.energy_min_meV = min(points.energy_meV);
panel.energy_max_meV = max(points.energy_meV);

exportgraphics(fig, panel.png, 'Resolution', 300);
exportgraphics(fig, panel.pdf, 'ContentType', 'vector');
close(fig);
end


function points = local_filter_points(points, qe, energy_window)
required = {'q_Ainv', 'energy_meV'};
missing = setdiff(required, points.Properties.VariableNames);
if ~isempty(missing)
    error('run_b1_current_waterfall_three_dataset_export:BadCsv', ...
        'B1 point table is missing columns: %s', strjoin(missing, ', '));
end
q_axis = double(qe.q_Ainv(:));
valid = isfinite(points.q_Ainv) & isfinite(points.energy_meV) & ...
    points.energy_meV >= energy_window(1) & ...
    points.energy_meV <= energy_window(2) & ...
    points.q_Ainv >= min(q_axis) & points.q_Ainv <= max(q_axis);
points = sortrows(points(valid, :), 'q_Ainv');
if isempty(points)
    error('run_b1_current_waterfall_three_dataset_export:NoPoints', ...
        'No finite B1 points remain after filtering.');
end
end


function [energy_axis, energy_mask] = local_energy_window(qe, window_meV)
full_energy = double(qe.energy_meV(:));
energy_mask = full_energy >= window_meV(1) & full_energy <= window_meV(2);
if ~any(energy_mask)
    energy_mask = full_energy >= max(0, min(full_energy));
end
energy_axis = full_energy(energy_mask);
end


function [q_values, q_columns, traces] = local_trace_set(qe, points, ...
    energy_mask, max_traces)
q_axis = double(qe.q_Ainv(:));
dq = median(abs(diff(unique(q_axis))), 'omitnan');
if ~isfinite(dq) || dq <= 0
    dq = 0;
end
q_min = min(points.q_Ainv) - 0.5 * dq;
q_max = max(points.q_Ainv) + 0.5 * dq;
q_columns = find(q_axis >= q_min & q_axis <= q_max);
if isempty(q_columns)
    [~, q_columns] = min(abs(q_axis - mean(points.q_Ainv, 'omitnan')));
end
max_traces = max(3, round(max_traces));
if numel(q_columns) > max_traces
    keep = unique(round(linspace(1, numel(q_columns), max_traces)));
    q_columns = q_columns(keep);
end
q_values = q_axis(q_columns);
traces = double(qe.intensity(energy_mask, q_columns));
end


function normalized = local_visual_normalize_traces(energy_axis, traces, ...
    norm_energy_axis, norm_traces, options)
normalized = local_area_normalize_traces(traces, norm_energy_axis, norm_traces);
normalized = local_apply_waterfall_residual_and_gain(energy_axis, normalized, ...
    options.waterfallGain);
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
        if nnz(valid) >= 2
            area = trapz(norm_energy_axis(valid), abs(y_norm(valid)));
        end
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
    traces, gain_value)
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
if isfinite(gain_value) && gain_value > 0
    traces = traces .* gain_value;
end
traces(~isfinite(energy_axis), :) = NaN;
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


function rows = local_overlay_b1_points(ax, points, energy_axis, normalized, ...
    offsets, q_values, q_columns, session)
rows = cell(height(points), 8);
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
    plot(ax, energy, y, 'o', 'MarkerSize', 5.8, ...
        'MarkerFaceColor', [0.00 0.40 1.00], ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.75);
    rows(i, :) = {char(session.label), char(session.folder), q, ...
        q_values(trace_idx), q_columns(trace_idx), energy, trace_idx, y};
end
rows = rows(~cellfun(@isempty, rows(:, 1)), :);
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
    canvas(y0:(y0 + h - 1), x0:(x0 + w - 1), :) = img;
    x0 = x0 + w + gap_px;
end
end


function local_write_report(out, sessions, panel_results, options)
fid = fopen(out.report, 'w');
if fid < 0
    error('run_b1_current_waterfall_three_dataset_export:ReportWriteFailed', ...
        'Could not write report: %s', out.report);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# B1 current-qcalib waterfall export\n\n');
fprintf(fid, 'Generated by `run_b1_current_waterfall_three_dataset_export`.\n\n');
fprintf(fid, '## Scope\n\n');
fprintf(fid, '- Uses current GUI-history result folders only.\n');
fprintf(fid, '- Reads `analysis_results.mat` spectra and `branch1_points.csv` for each dataset.\n');
fprintf(fid, '- Does not reuse archived PRE_QCALIB waterfall outputs.\n');
fprintf(fid, '- Uses signed-q traces inside the accepted B1 q range.\n');
fprintf(fid, '- Overlays only one marker series: the accepted B1 Fano apex points.\n');
fprintf(fid, '- Does not plot or use double-peak lower/upper branch markers.\n');
fprintf(fid, '- Area normalization, residual display, vertical offsets, and gain are display operations only.\n\n');
fprintf(fid, '## Display settings\n\n');
fprintf(fid, '- Energy window: %.0f-%.0f meV.\n', ...
    options.energyWindowMeV(1), options.energyWindowMeV(2));
fprintf(fid, '- Area normalization window: %.0f-%.0f meV.\n', ...
    options.areaNormWindowMeV(1), options.areaNormWindowMeV(2));
fprintf(fid, '- Waterfall gain: %.2f.\n', options.waterfallGain);
fprintf(fid, '- Maximum plotted traces per panel: %d.\n\n', ...
    options.maxTracesPerPanel);
fprintf(fid, '## Panels\n\n');
for i = 1:numel(sessions)
    p = panel_results{i};
    fprintf(fid, '- `%s`: `%s`, B1 points %d, plotted traces %d, q %.5f to %.5f 1/A, apex %.1f to %.1f meV.\n', ...
        sessions(i).label, sessions(i).folder, p.n_total_b1_points, ...
        p.n_plotted_traces, p.q_min_Ainv, p.q_max_Ainv, ...
        p.energy_min_meV, p.energy_max_meV);
end
fprintf(fid, '\n## Outputs\n\n');
fprintf(fid, '- PNG: `%s`\n', out.png);
fprintf(fid, '- PDF: `%s`\n', out.pdf);
fprintf(fid, '- CSV: `%s`\n', out.csv);
end
