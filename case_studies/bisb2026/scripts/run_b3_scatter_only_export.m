function output = run_b3_scatter_only_export()
%RUN_B3_SCATTER_ONLY_EXPORT Export a no-fit B3 scatter comparison.

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

results_root = fullfile(project_root, 'paper_results');
output_dir = fullfile(results_root, 'b3_scatter_only_260507');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

datasets = local_dataset_config(results_root);
combined = table();

for i = 1:numel(datasets)
    csv_path = fullfile(datasets(i).input_dir, 'branch3_points.csv');
    if ~isfile(csv_path)
        error('run_b3_scatter_only_export:MissingInput', ...
            'Missing B3 point CSV: %s', csv_path);
    end

    tbl = readtable(csv_path);
    local_require_columns(tbl, {'q_Ainv', 'energy_meV'});

    n = height(tbl);
    tbl.session_key = repmat({datasets(i).session_key}, n, 1);
    tbl.session_label = repmat({datasets(i).session_label}, n, 1);
    tbl.source_csv = repmat({csv_path}, n, 1);
    tbl.energy_eV = tbl.energy_meV ./ 1000;
    tbl = movevars(tbl, {'session_key', 'session_label', 'source_csv'}, ...
        'Before', 1);
    tbl = movevars(tbl, 'energy_eV', 'After', 'energy_meV');
    combined = [combined; tbl]; %#ok<AGROW>
end

combined_csv = fullfile(output_dir, 'b3_scatter_points_three_sessions.csv');
writetable(combined, combined_csv);

png_path = fullfile(output_dir, 'b3_scatter_three_sessions_no_fit.png');
pdf_path = fullfile(output_dir, 'b3_scatter_three_sessions_no_fit.pdf');
local_plot_b3_scatter(combined, datasets, png_path, pdf_path);

averaged = local_qabs_average(combined, datasets);
averaged_csv = fullfile(output_dir, 'b3_qabs_averaged_points.csv');
writetable(averaged, averaged_csv);

qabs_png_path = fullfile(output_dir, ...
    'b3_apex_dispersion_qabs_averaged_scatter_errorbars.png');
qabs_pdf_path = fullfile(output_dir, ...
    'b3_apex_dispersion_qabs_averaged_scatter_errorbars.pdf');
local_plot_b3_qabs_curve(averaged, datasets, qabs_png_path, qabs_pdf_path);

output = struct();
output.output_dir = output_dir;
output.combined_csv = combined_csv;
output.averaged_csv = averaged_csv;
output.png_path = png_path;
output.pdf_path = pdf_path;
output.qabs_png_path = qabs_png_path;
output.qabs_pdf_path = qabs_pdf_path;
output.n_points = height(combined);
output.n_qabs_points = height(averaged);

fprintf('B3 scatter-only export complete.\n');
fprintf('  Output directory: %s\n', output_dir);
fprintf('  Combined points: %d\n', output.n_points);
fprintf('  PNG: %s\n', png_path);
fprintf('  PDF: %s\n', pdf_path);
fprintf('  |q|-averaged points: %d\n', output.n_qabs_points);
fprintf('  |q|-averaged PNG: %s\n', qabs_png_path);
fprintf('  |q|-averaged PDF: %s\n', qabs_pdf_path);
end


function datasets = local_dataset_config(results_root)
datasets = struct( ...
    'session_key', {}, ...
    'session_label', {}, ...
    'input_dir', {}, ...
    'marker', {}, ...
    'color', {});

datasets(end + 1) = struct( ...
    'session_key', '590_PL2_10w', ...
    'session_label', '590 10w', ...
    'input_dir', fullfile(results_root, '590_gui_history_area_260506'), ...
    'marker', 'o', ...
    'color', [0.120, 0.470, 0.900]);

datasets(end + 1) = struct( ...
    'session_key', 'n0_PL2_10w_repeat', ...
    'session_label', 'n0 10w repeat', ...
    'input_dir', fullfile(results_root, 'n0_PL2_10w_gui_history_area_260506'), ...
    'marker', 'o', ...
    'color', [0.160, 0.500, 0.220]);

datasets(end + 1) = struct( ...
    'session_key', 'no_PL2_20w_2film', ...
    'session_label', '20w 2film', ...
    'input_dir', fullfile(results_root, ...
        'no_PL2_20w_2film_gui_history_area_260506_highq_refined'), ...
    'marker', 'o', ...
    'color', [0.930, 0.280, 0.300]);
end


function local_require_columns(tbl, required)
names = tbl.Properties.VariableNames;
for i = 1:numel(required)
    if ~ismember(required{i}, names)
        error('run_b3_scatter_only_export:MissingColumn', ...
            'Branch CSV is missing required column "%s".', required{i});
    end
end
end


function local_plot_b3_scatter(combined, datasets, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 620]);
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(datasets)
    mask = strcmp(combined.session_key, datasets(i).session_key);
    sub = combined(mask, :);
    scatter(ax, sub.q_Ainv, sub.energy_eV, 38, ...
        'Marker', datasets(i).marker, ...
        'MarkerEdgeColor', datasets(i).color, ...
        'MarkerFaceColor', datasets(i).color, ...
        'MarkerFaceAlpha', 0.70, ...
        'MarkerEdgeAlpha', 0.95, ...
        'LineWidth', 0.8, ...
        'DisplayName', datasets(i).session_label);
end

hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
ax.LineWidth = 0.9;
ax.FontName = 'Arial';
ax.FontSize = 11;

xlim(ax, [-0.155, 0.155]);
y_min = min(combined.energy_eV, [], 'omitnan');
y_max = max(combined.energy_eV, [], 'omitnan');
ylim(ax, [floor((y_min - 0.04) * 20) / 20, ceil((y_max + 0.04) * 20) / 20]);

xlabel(ax, 'q (A^{-1})');
ylabel(ax, 'B3 energy (eV)');
title(ax, 'B3 branch points, scatter only');
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off');

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function averaged = local_qabs_average(combined, datasets)
averaged = table();

for i = 1:numel(datasets)
    mask = strcmp(combined.session_key, datasets(i).session_key);
    sub = combined(mask, :);
    if isempty(sub)
        continue;
    end

    q_abs_group = round(abs(sub.q_Ainv), 6);
    q_values = unique(q_abs_group);
    q_values = q_values(q_values > 0);

    rows = table();
    for qi = 1:numel(q_values)
        q_mask = abs(q_abs_group - q_values(qi)) < 1e-12;
        q_sub = sub(q_mask, :);
        energy_eV = double(q_sub.energy_eV);
        finite_energy = energy_eV(isfinite(energy_eV));
        if isempty(finite_energy)
            continue;
        end

        ci_eV = local_ci_half_eV(q_sub);
        finite_ci = ci_eV(isfinite(ci_eV) & ci_eV >= 0);
        if isempty(finite_ci)
            ci_rms_eV = NaN;
        else
            ci_rms_eV = sqrt(mean(finite_ci .^ 2));
        end

        if numel(finite_energy) > 1
            pair_spread_eV = std(finite_energy, 0);
        else
            pair_spread_eV = NaN;
        end

        err_candidates = [ci_rms_eV, pair_spread_eV];
        err_candidates = err_candidates(isfinite(err_candidates));
        if isempty(err_candidates)
            err_eV = 0;
        else
            err_eV = max(err_candidates);
        end

        rows = [rows; table( ...
            {datasets(i).session_key}, ...
            {datasets(i).session_label}, ...
            q_values(qi), ...
            mean(finite_energy, 'omitnan'), ...
            err_eV, ...
            numel(finite_energy), ...
            min(q_sub.q_Ainv), ...
            max(q_sub.q_Ainv), ...
            'VariableNames', {'session_key', 'session_label', ...
            'q_abs_Ainv', 'energy_mean_eV', 'energy_err_eV', ...
            'n_raw_points', 'q_signed_min_Ainv', ...
            'q_signed_max_Ainv'})]; %#ok<AGROW>
    end

    rows = sortrows(rows, 'q_abs_Ainv');
    averaged = [averaged; rows]; %#ok<AGROW>
end
end


function ci_eV = local_ci_half_eV(tbl)
if ismember('E_ci_half_meV', tbl.Properties.VariableNames)
    ci_eV = double(tbl.E_ci_half_meV) ./ 1000;
elseif all(ismember({'E_ci_lo', 'E_ci_hi'}, tbl.Properties.VariableNames))
    ci_eV = abs(double(tbl.E_ci_hi) - double(tbl.E_ci_lo)) ./ 2000;
else
    ci_eV = NaN(height(tbl), 1);
end
end


function local_plot_b3_qabs_curve(averaged, datasets, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1024, 760]);
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(datasets)
    mask = strcmp(averaged.session_key, datasets(i).session_key);
    sub = averaged(mask, :);
    if isempty(sub)
        continue;
    end

    [q_sorted, order] = sort(sub.q_abs_Ainv);
    E_sorted = sub.energy_mean_eV(order);
    err_sorted = sub.energy_err_eV(order);
    col = datasets(i).color;

    errorbar(ax, q_sorted, E_sorted, err_sorted, ...
        'LineStyle', 'none', ...
        'Marker', 'none', ...
        'Color', local_lighten(col, 0.62), ...
        'LineWidth', 0.55, ...
        'CapSize', 0, ...
        'HandleVisibility', 'off');

    scatter(ax, q_sorted, E_sorted, 34, ...
        'Marker', datasets(i).marker, ...
        'MarkerFaceColor', col, ...
        'MarkerEdgeColor', col, ...
        'LineWidth', 0.5, ...
        'DisplayName', datasets(i).session_label);
end

hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.12;
ax.LineWidth = 0.9;
ax.FontName = 'Arial';
ax.FontSize = 10;
ax.XColor = [0.12, 0.18, 0.30];
ax.YColor = [0.12, 0.18, 0.30];

xlim(ax, [0, 0.15]);
y_min = min(averaged.energy_mean_eV - averaged.energy_err_eV, [], 'omitnan');
y_max = max(averaged.energy_mean_eV + averaged.energy_err_eV, [], 'omitnan');
ylim(ax, [floor((y_min - 0.02) * 20) / 20, ceil((y_max + 0.02) * 20) / 20]);

xticks(ax, 0:0.05:0.15);
xlabel(ax, '|q| (A^{-1})', 'FontWeight', 'bold', 'FontSize', 12);
ylabel(ax, 'E (eV)', 'FontWeight', 'bold', 'FontSize', 12);
text(ax, 0.02, 0.96, 'B3 apex dispersion, |q|-averaged', ...
    'Units', 'normalized', ...
    'FontName', 'Arial', ...
    'FontWeight', 'bold', ...
    'FontSize', 11, ...
    'Color', [0.08, 0.10, 0.14], ...
    'VerticalAlignment', 'top');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9);

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function color = local_lighten(color_in, amount)
color = color_in + amount .* (1 - color_in);
color = min(max(color, 0), 1);
end
