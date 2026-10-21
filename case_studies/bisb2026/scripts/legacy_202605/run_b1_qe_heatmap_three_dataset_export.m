function out = run_b1_qe_heatmap_three_dataset_export(options)
%RUN_B1_QE_HEATMAP_THREE_DATASET_EXPORT Export B1-window q-E heatmaps only.
%
% Thesis-facing export layer. It reads current GUI-history analysis outputs,
% crops the existing area-normalized physical q-E maps to the B1 discussion
% energy window, and deliberately does not overlay fitted branch points.

arguments
    options.outputTag {mustBeTextScalar} = "260521"
    options.energyWindowMeV (1, 2) double = [200 1700]
    options.useSnapQWindow (1, 1) logical = true
    options.renormalizeDisplayWindow (1, 1) logical = true
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

results_root = fullfile(project_root, 'paper_results');
out_dir = fullfile(results_root, ...
    sprintf('b1_qe_heatmaps_three_dataset_%s', char(options.outputTag)));
if ~isfolder(out_dir)
    mkdir(out_dir);
end

sessions = local_sessions();
panel_results = cell(1, numel(sessions));
summary_rows = cell(numel(sessions), 1);

for i = 1:numel(sessions)
    [panel_results{i}, summary_rows{i}] = local_plot_session( ...
        results_root, out_dir, sessions(i), options);
end

summary = vertcat(summary_rows{:});
summary_csv = fullfile(out_dir, 'b1_qe_heatmap_summary.csv');
writetable(summary, summary_csv);

combined = local_plot_combined(out_dir, sessions, panel_results, options);

out = struct();
out.output_dir = out_dir;
out.summary_csv = summary_csv;
out.combined_png = combined.png;
out.combined_pdf = combined.pdf;
out.panels = panel_results;
out.report = fullfile(out_dir, 'b1_qe_heatmap_report.md');
local_write_report(out, sessions, panel_results, summary, options);

fprintf('B1 q-E heatmap export written:\n');
fprintf('  %s\n', out.combined_png);
fprintf('  %s\n', out.combined_pdf);
fprintf('  %s\n', out.summary_csv);
fprintf('  %s\n', out.report);
end


function sessions = local_sessions()
sessions = struct( ...
    'label', { ...
        '10w 1film', ...
        '10w 1film repeat', ...
        '20w 2film'}, ...
    'folder', { ...
        '590_gui_history_area_260506', ...
        'n0_PL2_10w_gui_history_area_260506', ...
        'no_PL2_20w_2film_gui_history_area_260506_highq_refined'}, ...
    'safe_name', { ...
        '590', ...
        'n0_10w_repeat', ...
        '20w_2film'});
end


function [panel, summary] = local_plot_session(results_root, out_dir, ...
    session, options)
folder = fullfile(results_root, session.folder);
mat_path = fullfile(folder, 'analysis_results.mat');
if ~isfile(mat_path)
    error('run_b1_qe_heatmap_three_dataset_export:MissingMat', ...
        'Missing analysis MAT: %s', mat_path);
end

saved = load(mat_path, 'output');
qe = saved.output.qe_pp;
snap = saved.output.snap;

[q_axis, e_axis, map, q_window, e_window] = local_cropped_map(qe, snap, ...
    options);
clim_vals = local_color_limits(map);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 80 960 700]);
ax = axes(fig);
imagesc(ax, q_axis, e_axis, map);
axis(ax, 'xy');
colormap(ax, turbo);
if all(isfinite(clim_vals)) && clim_vals(1) < clim_vals(2)
    clim(ax, clim_vals);
end
colorbar(ax);
box(ax, 'on');
ax.LineWidth = 1.1;
ax.FontSize = 16;
ax.TickDir = 'out';
xlabel(ax, 'q (1/A)');
ylabel(ax, 'Energy relative to ZLP (meV)');

panel = struct();
panel.label = session.label;
panel.folder = session.folder;
panel.png = fullfile(out_dir, sprintf('b1_qe_heatmap_%s.png', ...
    session.safe_name));
panel.pdf = fullfile(out_dir, sprintf('b1_qe_heatmap_%s.pdf', ...
    session.safe_name));
panel.q_min_Ainv = min(q_axis);
panel.q_max_Ainv = max(q_axis);
panel.energy_min_meV = min(e_axis);
panel.energy_max_meV = max(e_axis);
panel.clim_low = clim_vals(1);
panel.clim_high = clim_vals(2);
panel.source_mat = mat_path;

local_export_figure(fig, panel.png, panel.pdf);
close(fig);

summary = table( ...
    string(session.label), string(session.folder), string(mat_path), ...
    q_window(1), q_window(2), e_window(1), e_window(2), ...
    min(q_axis), max(q_axis), min(e_axis), max(e_axis), ...
    size(map, 2), size(map, 1), options.renormalizeDisplayWindow, ...
    clim_vals(1), clim_vals(2), ...
    string(panel.png), string(panel.pdf), ...
    'VariableNames', { ...
        'dataset_label', 'session_folder', 'source_mat', ...
        'requested_q_min_Ainv', 'requested_q_max_Ainv', ...
        'requested_energy_min_meV', 'requested_energy_max_meV', ...
        'plotted_q_min_Ainv', 'plotted_q_max_Ainv', ...
        'plotted_energy_min_meV', 'plotted_energy_max_meV', ...
        'n_q_columns', 'n_energy_channels', 'display_window_area_norm', ...
        'color_limit_low', 'color_limit_high', ...
        'panel_png', 'panel_pdf'});
end


function combined = local_plot_combined(out_dir, sessions, panel_results, options)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60 80 1500 1080]);
positions = [
    0.075 0.565 0.360 0.350
    0.565 0.565 0.360 0.350
    0.320 0.105 0.360 0.350
    ];
for i = 1:numel(panel_results)
    saved = load(fullfile(fileparts(panel_results{i}.source_mat), ...
        'analysis_results.mat'), 'output');
    [q_axis, e_axis, map] = local_cropped_map(saved.output.qe_pp, ...
        saved.output.snap, options);
    ax = axes(fig, 'Position', positions(i, :)); %#ok<LAXES>
    imagesc(ax, q_axis, e_axis, map);
    axis(ax, 'xy');
    colormap(ax, turbo);
    clim_vals = [panel_results{i}.clim_low, panel_results{i}.clim_high];
    if all(isfinite(clim_vals)) && clim_vals(1) < clim_vals(2)
        clim(ax, clim_vals);
    end
    box(ax, 'on');
    ax.LineWidth = 1.0;
    ax.FontSize = 13;
    ax.TickDir = 'out';
    pbaspect(ax, [1.35 1 1]);
    title(ax, sprintf('(%c) %s', char(96 + i), sessions(i).label), ...
        'Interpreter', 'none');
    xlabel(ax, 'q (1/A)');
    if i == 1
        ylabel(ax, 'Energy relative to ZLP (meV)');
    else
        ylabel(ax, '');
    end
    colorbar(ax);
end

combined = struct();
combined.png = fullfile(out_dir, 'b1_qe_heatmaps_three_dataset.png');
combined.pdf = fullfile(out_dir, 'b1_qe_heatmaps_three_dataset.pdf');
local_export_figure(fig, combined.png, combined.pdf);
close(fig);
end


function [q_axis, e_axis, map, q_window, e_window] = local_cropped_map(qe, ...
    snap, options)
full_q = double(qe.q_Ainv(:)).';
full_e = double(qe.energy_meV(:));
if options.useSnapQWindow
    q_window = sort(double([snap.qStart, snap.qEnd]));
else
    q_window = [min(full_q), max(full_q)];
end
e_window = sort(double(options.energyWindowMeV));

q_mask = full_q >= q_window(1) & full_q <= q_window(2);
e_mask = full_e >= e_window(1) & full_e <= e_window(2);
if ~any(q_mask)
    error('run_b1_qe_heatmap_three_dataset_export:EmptyQWindow', ...
        'No q columns in requested window [%.5g, %.5g].', ...
        q_window(1), q_window(2));
end
if ~any(e_mask)
    error('run_b1_qe_heatmap_three_dataset_export:EmptyEnergyWindow', ...
        'No energy channels in requested window [%.1f, %.1f] meV.', ...
        e_window(1), e_window(2));
end

q_axis = full_q(q_mask);
e_axis = full_e(e_mask);
map = double(qe.intensity(e_mask, q_mask));
if options.renormalizeDisplayWindow
    map = local_area_normalize_columns(e_axis, map);
end
end


function map_norm = local_area_normalize_columns(e_axis, map)
map_norm = zeros(size(map));
for j = 1:size(map, 2)
    y = double(map(:, j));
    valid = isfinite(e_axis) & isfinite(y);
    area = NaN;
    if nnz(valid) >= 2
        area = trapz(e_axis(valid), y(valid));
    end
    if ~isfinite(area) || abs(area) <= eps
        if nnz(valid) >= 2
            area = trapz(e_axis(valid), abs(y(valid)));
        end
    end
    if ~isfinite(area) || abs(area) <= eps
        area = 1;
    end
    map_norm(:, j) = y ./ area;
end
end


function clim_vals = local_color_limits(map)
vals = map(isfinite(map));
if isempty(vals)
    clim_vals = [NaN NaN];
    return
end
vals = sort(vals(:));
lo = local_percentile(vals, 2);
hi = local_percentile(vals, 98);
if lo == hi
    hi = lo + eps;
end
clim_vals = [lo hi];
end


function value = local_percentile(sorted_vals, pct)
n = numel(sorted_vals);
if n == 1
    value = sorted_vals(1);
    return
end
pos = 1 + (n - 1) * pct / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    value = sorted_vals(lo);
else
    frac = pos - lo;
    value = sorted_vals(lo) * (1 - frac) + sorted_vals(hi) * frac;
end
end


function local_export_figure(fig, png_path, pdf_path)
exportgraphics(fig, png_path, 'Resolution', 300);
try
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
catch
    print(fig, pdf_path, '-dpdf', '-bestfit');
end
end


function local_write_report(out, sessions, panel_results, summary, options)
fid = fopen(out.report, 'w');
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# B1 q-E heatmap-only export\n\n');
fprintf(fid, '## Scope\n\n');
fprintf(fid, '- Purpose: thesis-facing B1 discussion heatmaps.\n');
fprintf(fid, '- Overlay policy: no fitted branch points, no candidate points, no model curves.\n');
fprintf(fid, '- Energy display window: `[%.0f, %.0f] meV`.\n', ...
    options.energyWindowMeV(1), options.energyWindowMeV(2));
fprintf(fid, '- q display window: current GUI snap q window from each analysis result.\n');
if options.renormalizeDisplayWindow
    fprintf(fid, '- Area normalization: columns renormalized within the displayed energy window.\n');
else
    fprintf(fid, '- Area normalization: inherited from the current GUI-history output without display-window renormalization.\n');
end
fprintf(fid, '- Input level: current derived `analysis_results.mat` files; no raw-data reprocessing.\n');
fprintf(fid, '- Intensity preprocessing before display-window normalization: inherited from current GUI-history outputs.\n\n');

fprintf(fid, '## Outputs\n\n');
fprintf(fid, '- Combined PNG: `%s`\n', out.combined_png);
fprintf(fid, '- Combined PDF: `%s`\n', out.combined_pdf);
fprintf(fid, '- Summary CSV: `%s`\n\n', out.summary_csv);

fprintf(fid, '## Panels\n\n');
fprintf(fid, '| Panel | Dataset | Source result folder | Plotted q range (1/A) | Plotted energy range (meV) | Panel PNG |\n');
fprintf(fid, '|---:|---|---|---:|---:|---|\n');
for i = 1:numel(panel_results)
    p = panel_results{i};
    fprintf(fid, '| %d | %s | `%s` | %.5g to %.5g | %.0f to %.0f | `%s` |\n', ...
        i, sessions(i).label, sessions(i).folder, ...
        p.q_min_Ainv, p.q_max_Ainv, ...
        p.energy_min_meV, p.energy_max_meV, p.png);
end

fprintf(fid, '\n## Color scaling note\n\n');
fprintf(fid, ['Each panel uses its own 2nd-98th percentile display color limits ', ...
    'to keep the B1-window morphology readable. Therefore the heatmaps are ', ...
    'qualitative spatial/spectral views, not a basis for comparing absolute ', ...
    'spectral weight across datasets.\n\n']);

fprintf(fid, '## Summary table preview\n\n');
fprintf(fid, '- Rows written: %d\n', height(summary));
end
