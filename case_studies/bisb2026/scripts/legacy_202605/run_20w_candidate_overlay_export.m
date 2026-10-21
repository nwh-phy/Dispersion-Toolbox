function output = run_20w_candidate_overlay_export()
%RUN_20W_CANDIDATE_OVERLAY_EXPORT Export conservative and exploratory B1 views.
%
% This script is a presentation/export layer. It does not alter the GUI
% pipeline, the branch CSVs, or the quasi-2D fits. Candidate-only points are
% loaded from the high-q audit result and explicitly marked as excluded from
% the main fit.

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

main_dir = fullfile(project_root, 'paper_results', ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined');
audit_dir = fullfile(project_root, 'paper_results', ...
    '20w_B1_highq_audit_260506');
output_dir = fullfile(project_root, 'paper_results', ...
    '20w_B1_candidate_overlay_260506');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

paths = local_input_paths(main_dir, audit_dir);
result_data = load(paths.analysis_results, 'output');
candidate_tbl = local_candidate_points_table(readtable(paths.rescue_candidates));
model_tbl = readtable(paths.dispersion_model_summary);
fit_model = local_current_b1_model(model_tbl);

candidate_csv = fullfile(output_dir, 'branch1_candidate_points.csv');
writetable(candidate_tbl, candidate_csv);

fig_paths = struct();
fig_paths.conservative_qe_map = fullfile(output_dir, ...
    '20w_B1_conservative_qe_map.png');
fig_paths.exploratory_qe_map = fullfile(output_dir, ...
    '20w_B1_exploratory_qe_map_candidate_overlay.png');
fig_paths.conservative_dispersion = fullfile(output_dir, ...
    '20w_B1_conservative_dispersion.png');
fig_paths.exploratory_dispersion = fullfile(output_dir, ...
    '20w_B1_exploratory_dispersion_candidate_overlay.png');

qe_pp = result_data.output.qe_pp;
branches = result_data.output.branches;
snap = result_data.output.snap;

local_plot_qe_map_overlay(qe_pp, branches, candidate_tbl, snap, ...
    fig_paths.conservative_qe_map, false);
local_plot_qe_map_overlay(qe_pp, branches, candidate_tbl, snap, ...
    fig_paths.exploratory_qe_map, true);
local_plot_dispersion_overlay(branches, candidate_tbl, fit_model, ...
    fig_paths.conservative_dispersion, false);
local_plot_dispersion_overlay(branches, candidate_tbl, fit_model, ...
    fig_paths.exploratory_dispersion, true);

report_path = fullfile(output_dir, '20w_B1_candidate_overlay_report.md');
local_write_report(report_path, main_dir, audit_dir, candidate_tbl, ...
    fit_model, fig_paths);

mat_path = fullfile(output_dir, '20w_B1_candidate_overlay_results.mat');
save(mat_path, 'candidate_tbl', 'fit_model', 'fig_paths');

output = struct();
output.output_dir = output_dir;
output.candidate_csv = candidate_csv;
output.report_path = report_path;
output.figure_paths = fig_paths;
output.candidate_points = candidate_tbl;

fprintf('20w B1 candidate overlay export complete.\n');
fprintf('  Output: %s\n', output_dir);
fprintf('  Candidate-only points: %d\n', height(candidate_tbl));
end


function paths = local_input_paths(main_dir, audit_dir)
paths = struct();
paths.analysis_results = fullfile(main_dir, 'analysis_results.mat');
paths.dispersion_model_summary = fullfile(main_dir, ...
    'dispersion_model_summary.csv');
paths.rescue_candidates = fullfile(audit_dir, ...
    'branch1_highq_rescue_candidates.csv');

names = fieldnames(paths);
for i = 1:numel(names)
    if ~isfile(paths.(names{i}))
        error('run_20w_candidate_overlay_export:MissingInput', ...
            'Missing required input: %s', paths.(names{i}));
    end
end
end


function candidate_tbl = local_candidate_points_table(rescue_tbl)
if isempty(rescue_tbl)
    candidate_tbl = table();
    return
end

candidate_tbl = table();
candidate_tbl.q_Ainv = rescue_tbl.q_Ainv;
candidate_tbl.energy_meV = rescue_tbl.candidate_energy_meV;
candidate_tbl.E_ci_half_meV = rescue_tbl.old_CI_half_meV;
candidate_tbl.R2 = rescue_tbl.old_R2;
candidate_tbl.gamma_over_E = rescue_tbl.old_gamma_over_E;
candidate_tbl.model_residual_abs_meV = rescue_tbl.model_residual_abs_meV;
candidate_tbl.symmetry_abs_delta_meV = rescue_tbl.symmetry_abs_delta_meV;
candidate_tbl.source = rescue_tbl.candidate_source;
candidate_tbl.evidence_level = rescue_tbl.candidate_class;
candidate_tbl.included_in_main_fit = false(height(candidate_tbl), 1);
candidate_tbl.plot_layer = repmat({'candidate_only_not_for_fit'}, ...
    height(candidate_tbl), 1);
candidate_tbl.note = repmat({ ...
    'candidate-only points are not included in quasi-2D fitting'}, ...
    height(candidate_tbl), 1);
end


function fit_model = local_current_b1_model(model_tbl)
mask = model_tbl.branch == 1 & strcmp(model_tbl.model, 'quasi2d_plasmon') & ...
    model_tbl.success == 1;
if ~any(mask)
    error('run_20w_candidate_overlay_export:MissingModel', ...
        'Could not find successful current B1 quasi2d model.');
end
row = model_tbl(find(mask, 1), :);
fit_model = struct();
fit_model.A_fit = row.param1(1);
fit_model.rho0_A = row.rho0_A(1);
fit_model.q_c_Ainv = row.q_c_Ainv(1);
fit_model.E_flat_meV = row.E_flat_meV(1);
fit_model.R_squared = row.R2(1);
fit_model.epsilon_bg = 1;
end


function local_plot_qe_map_overlay(qe, branches, candidate_tbl, snap, out_path, show_candidates)
q_mask = qe.q_Ainv >= min(snap.qStart, snap.qEnd) & ...
    qe.q_Ainv <= max(snap.qStart, snap.qEnd);
e_mask = qe.energy_meV >= min(snap.energyMin, snap.energyMax) & ...
    qe.energy_meV <= max(snap.energyMin, snap.energyMax);
map = double(qe.intensity(e_mask, q_mask));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 900]);
ax = axes(fig);
imagesc(ax, qe.q_Ainv(q_mask), qe.energy_meV(e_mask), map);
axis(ax, 'xy');
colormap(ax, turbo);
clim_vals = local_color_limits(map);
if all(isfinite(clim_vals)) && clim_vals(1) < clim_vals(2)
    clim(ax, clim_vals);
end
colorbar(ax);
hold(ax, 'on');

for b = 1:numel(branches)
    br = branches{b};
    if isempty(br)
        continue
    end
    col = local_branch_color(b);
    plot(ax, br(:, 1), br(:, 2), 'o', ...
        'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
        'Color', col, 'MarkerSize', 5.5, ...
        'LineStyle', 'none', 'DisplayName', sprintf('Branch %d current', b));
end

if show_candidates && ~isempty(candidate_tbl)
    candidate_col = local_candidate_color();
    errorbar(ax, candidate_tbl.q_Ainv, candidate_tbl.energy_meV, ...
        candidate_tbl.E_ci_half_meV, 's', ...
        'Color', candidate_col, ...
        'MarkerEdgeColor', candidate_col, ...
        'MarkerFaceColor', 'none', ...
        'LineStyle', 'none', 'LineWidth', 1.0, 'MarkerSize', 5.5, ...
        'DisplayName', 'B1 candidate-only');
end

hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, 'q (1/A)');
ylabel(ax, 'Energy relative to ZLP (meV)');
if show_candidates
    title(ax, '20w q-E map with B1 candidate-only overlay');
else
    title(ax, '20w q-E map, conservative current branches');
end
legend(ax, 'Location', 'best', 'FontSize', 7);
local_export_figure(fig, out_path);
end


function local_plot_dispersion_overlay(branches, candidate_tbl, fit_model, out_path, show_candidates)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 900]);
ax = axes(fig);
hold(ax, 'on');

for b = 1:numel(branches)
    br = branches{b};
    if isempty(br)
        continue
    end
    col = local_branch_color(b);
    ci_half = local_branch_ci_half(br);
    errorbar(ax, br(:, 1), br(:, 2), ci_half, 'o', ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'w', ...
        'LineStyle', 'none', 'LineWidth', 1.0, 'MarkerSize', 5.5, ...
        'DisplayName', sprintf('Branch %d current', b));
end

q_fit = linspace(-0.015, 0.015, 301)';
E_fit = local_predict_quasi2d(fit_model, abs(q_fit));
plot(ax, q_fit, E_fit, 'k-', 'LineWidth', 1.6, ...
    'DisplayName', sprintf('Current B1 fit R^2=%.3f', fit_model.R_squared));

if show_candidates && ~isempty(candidate_tbl)
    candidate_col = local_candidate_color();
    errorbar(ax, candidate_tbl.q_Ainv, candidate_tbl.energy_meV, ...
        candidate_tbl.E_ci_half_meV, 's', ...
        'Color', candidate_col, ...
        'MarkerEdgeColor', candidate_col, ...
        'MarkerFaceColor', 'none', ...
        'LineStyle', 'none', 'LineWidth', 1.0, 'MarkerSize', 5.5, ...
        'DisplayName', 'B1 candidate-only, excluded from fit');
end

hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, 'q (1/A)');
ylabel(ax, 'Energy (meV)');
if show_candidates
    title(ax, '20w B1 exploratory candidate-only overlay');
else
    title(ax, '20w B1 conservative current dispersion');
end
legend(ax, 'Location', 'bestoutside', 'FontSize', 8);
local_export_figure(fig, out_path);
end


function ci_half = local_branch_ci_half(br)
ci_half = NaN(size(br, 1), 1);
if size(br, 2) >= 7
    ci_half = 0.5 .* (br(:, 7) - br(:, 6));
end
fallback = max(1, 0.001 .* abs(br(:, 2)));
bad = ~isfinite(ci_half) | ci_half <= 0;
ci_half(bad) = fallback(bad);
end


function E = local_predict_quasi2d(fit_model, q_abs)
E = sqrt(fit_model.A_fit .* q_abs ./ ...
    (fit_model.epsilon_bg + fit_model.rho0_A .* q_abs));
end


function col = local_branch_color(branch_id)
colors = [
    0.08 0.36 0.78
    0.12 0.65 0.30
    0.90 0.58 0.05
    0.55 0.25 0.75
    ];
idx = min(max(branch_id, 1), size(colors, 1));
col = colors(idx, :);
end


function col = local_candidate_color()
col = [0.72 0.10 0.48];
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


function p = local_percentile(sorted_vals, pct)
n = numel(sorted_vals);
if n == 1
    p = sorted_vals(1);
    return
end
pos = 1 + (n - 1) * pct / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    p = sorted_vals(lo);
else
    frac = pos - lo;
    p = sorted_vals(lo) * (1 - frac) + sorted_vals(hi) * frac;
end
end


function local_write_report(report_path, main_dir, audit_dir, candidate_tbl, ...
    fit_model, fig_paths)
fid = fopen(report_path, 'w');
if fid < 0
    error('run_20w_candidate_overlay_export:ReportOpenFailed', ...
        'Unable to write report: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# 20w B1 Candidate-only Overlay Export\n\n');
fprintf(fid, 'Main result directory: `%s`\n\n', main_dir);
fprintf(fid, 'Audit source directory: `%s`\n\n', audit_dir);
fprintf(fid, '## Rule\n\n');
fprintf(fid, 'The conservative figures show only current accepted branch points. ');
fprintf(fid, 'The exploratory figures add B1 high-q candidate-only points. ');
fprintf(fid, 'These candidate-only points are not included in quasi-2D fitting.\n\n');

fprintf(fid, '## Candidate-only layer\n\n');
fprintf(fid, '- Candidate-only points: %d\n', height(candidate_tbl));
fprintf(fid, '- `included_in_main_fit`: false for every candidate-only row\n');
fprintf(fid, '- Current B1 fit: rho0 = %.3g A, qc = %.4g 1/A, Eflat = %.1f meV, R2 = %.4f\n\n', ...
    fit_model.rho0_A, fit_model.q_c_Ainv, fit_model.E_flat_meV, ...
    fit_model.R_squared);

fprintf(fid, '## Outputs\n\n');
fprintf(fid, '- `branch1_candidate_points.csv`\n');
fprintf(fid, '- `%s`\n', local_file_name(fig_paths.conservative_qe_map));
fprintf(fid, '- `%s`\n', local_file_name(fig_paths.exploratory_qe_map));
fprintf(fid, '- `%s`\n', local_file_name(fig_paths.conservative_dispersion));
fprintf(fid, '- `%s`\n\n', local_file_name(fig_paths.exploratory_dispersion));

fprintf(fid, '## Recommended use\n\n');
fprintf(fid, 'Use conservative figures for the main result. ');
fprintf(fid, 'Use exploratory figures only when discussing possible B1 high-q continuation with the advisor.\n');
end


function name = local_file_name(path_value)
[~, name, ext] = fileparts(path_value);
name = [name ext];
end


function local_export_figure(fig, out_path)
set(fig, 'PaperPositionMode', 'auto');
print(fig, out_path, '-dpng', '-r300');
close(fig);
end
