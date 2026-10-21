function output = run_20w_lowq_sqrt_zero_diagnostic()
%RUN_20W_LOWQ_SQRT_ZERO_DIAGNOSTIC Test low-q sqrt zero crossing for 20w B1.
%
% This diagnostic reads the q-averaged 20w B1 points from the existing
% no-background B1 enhancement table. It does not modify the principal
% thickness-constrained quasi-2D plasmon fit.

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

input_csv = fullfile(project_root, 'paper_results', ...
    'b1_physical_fit_enhancements_260507', ...
    'b1_enhancement_points_qabs.csv');
output_dir = fullfile(project_root, 'paper_results', ...
    '20w_lowq_sqrt_zero_260508');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

qmax_values_Ainv = [0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.0475];
primary_qmax_Ainv = 0.040;
sigma_floor_meV = 10;
session_key = "no_PL2_20w_2film";

points = local_load_20w_points(input_csv, session_key, ...
    max(qmax_values_Ainv), sigma_floor_meV);
summary = local_sqrt_zero_sweep(points, qmax_values_Ainv, ...
    sigma_floor_meV);

points_csv = fullfile(output_dir, '20w_lowq_sqrt_zero_points.csv');
summary_csv = fullfile(output_dir, '20w_lowq_sqrt_zero_summary.csv');
writetable(points, points_csv);
writetable(summary, summary_csv);

fig_png = fullfile(output_dir, '20w_lowq_sqrt_zero_diagnostic.png');
fig_pdf = fullfile(output_dir, '20w_lowq_sqrt_zero_diagnostic.pdf');
local_plot_primary_fit(points, summary, primary_qmax_Ainv, ...
    fig_png, fig_pdf);

report_path = fullfile(output_dir, '20w_lowq_sqrt_zero_report.md');
local_write_report(report_path, input_csv, summary, primary_qmax_Ainv, ...
    fig_png);

mat_path = fullfile(output_dir, '20w_lowq_sqrt_zero_results.mat');
save(mat_path, 'points', 'summary', 'qmax_values_Ainv', ...
    'primary_qmax_Ainv', 'sigma_floor_meV', 'session_key');

output = struct();
output.output_dir = output_dir;
output.points_csv = points_csv;
output.summary_csv = summary_csv;
output.fig_png = fig_png;
output.fig_pdf = fig_pdf;
output.report_path = report_path;
output.mat_path = mat_path;
output.points = points;
output.summary = summary;

fprintf('20w low-q sqrt zero diagnostic complete.\n');
fprintf('  Output: %s\n', output_dir);
fprintf('  Points used up to %.4f A^-1: %d\n', ...
    primary_qmax_Ainv, sum(points.q_abs_Ainv <= primary_qmax_Ainv));
end


function points = local_load_20w_points(input_csv, session_key, ...
    qmax_Ainv, sigma_floor_meV)
if ~isfile(input_csv)
    error('run_20w_lowq_sqrt_zero_diagnostic:MissingInput', ...
        'Missing required input CSV: %s', input_csv);
end

tbl = readtable(input_csv);
required = {'session_key', 'q_abs_Ainv', 'energy_mean_meV', ...
    'energy_err_meV', 'R2_mean_all', 'include_for_fit'};
local_require_columns(tbl, required);

mask = strcmp(string(tbl.session_key), session_key) & ...
    tbl.include_for_fit == 1 & ...
    isfinite(tbl.q_abs_Ainv) & tbl.q_abs_Ainv > 0 & ...
    tbl.q_abs_Ainv <= qmax_Ainv & ...
    isfinite(tbl.energy_mean_meV) & tbl.energy_mean_meV > 0;
points = tbl(mask, :);
points = sortrows(points, 'q_abs_Ainv');
points.sqrt_q_Ainv = sqrt(points.q_abs_Ainv);
points.sigma_used_meV = max(points.energy_err_meV, sigma_floor_meV);

if height(points) < 3
    error('run_20w_lowq_sqrt_zero_diagnostic:InsufficientPoints', ...
        'Need at least 3 low-q 20w B1 points.');
end
end


function local_require_columns(tbl, required)
names = tbl.Properties.VariableNames;
for i = 1:numel(required)
    if ~ismember(required{i}, names)
        error('run_20w_lowq_sqrt_zero_diagnostic:MissingColumn', ...
            'Input table is missing required column "%s".', required{i});
    end
end
end


function summary = local_sqrt_zero_sweep(points, qmax_values_Ainv, ...
    sigma_floor_meV)
rows = table();
for i = 1:numel(qmax_values_Ainv)
    qmax = qmax_values_Ainv(i);
    sub = points(points.q_abs_Ainv <= qmax, :);
    if height(sub) < 3
        continue;
    end

    free_fit = local_fit_sqrt_free_intercept(sub, sigma_floor_meV);
    zero_fit = local_fit_sqrt_zero_intercept(sub, sigma_floor_meV);

    row = table( ...
        qmax, height(sub), min(sub.q_abs_Ainv), max(sub.q_abs_Ainv), ...
        free_fit.slope_meV_sqrtA, free_fit.intercept_meV, ...
        free_fit.intercept_ci95_low_meV, ...
        free_fit.intercept_ci95_high_meV, ...
        free_fit.intercept_sigma_meV, ...
        free_fit.intercept_z, ...
        free_fit.RMSE_meV, free_fit.reduced_chi2, free_fit.R_squared, ...
        zero_fit.slope_meV_sqrtA, zero_fit.RMSE_meV, ...
        zero_fit.reduced_chi2, zero_fit.R_squared, ...
        zero_fit.RMSE_meV - free_fit.RMSE_meV, ...
        'VariableNames', {'qmax_Ainv', 'n_points', ...
        'q_abs_min_Ainv', 'q_abs_max_data_Ainv', ...
        'free_slope_meV_sqrtA', 'free_intercept_meV', ...
        'free_intercept_ci95_low_meV', ...
        'free_intercept_ci95_high_meV', ...
        'free_intercept_sigma_meV', 'free_intercept_z', ...
        'free_RMSE_meV', 'free_reduced_chi2', 'free_R_squared', ...
        'zero_slope_meV_sqrtA', 'zero_RMSE_meV', ...
        'zero_reduced_chi2', 'zero_R_squared', 'delta_RMSE_zero_minus_free_meV'});
    rows = [rows; row]; %#ok<AGROW>
end
summary = rows;
end


function fit = local_fit_sqrt_free_intercept(points, sigma_floor_meV)
x = double(points.sqrt_q_Ainv(:));
y = double(points.energy_mean_meV(:));
sigma = max(double(points.energy_err_meV(:)), sigma_floor_meV);
w = 1 ./ (sigma .^ 2);

X = [x, ones(size(x))];
W = diag(w);
beta = (X' * W * X) \ (X' * W * y);
residuals = y - X * beta;

n = numel(y);
dof = max(n - 2, 1);
chi2 = sum((residuals ./ sigma) .^ 2);
reduced_chi2 = chi2 / dof;
cov_beta = pinv(X' * W * X) * reduced_chi2;
intercept_sigma = sqrt(max(cov_beta(2, 2), 0));

fit = struct();
fit.slope_meV_sqrtA = beta(1);
fit.intercept_meV = beta(2);
fit.intercept_sigma_meV = intercept_sigma;
fit.intercept_ci95_low_meV = beta(2) - 1.96 * intercept_sigma;
fit.intercept_ci95_high_meV = beta(2) + 1.96 * intercept_sigma;
fit.intercept_z = beta(2) / max(intercept_sigma, eps);
fit.RMSE_meV = sqrt(mean(residuals .^ 2));
fit.reduced_chi2 = reduced_chi2;
fit.R_squared = local_r_squared(y, residuals);
fit.residuals_meV = residuals;
end


function fit = local_fit_sqrt_zero_intercept(points, sigma_floor_meV)
x = double(points.sqrt_q_Ainv(:));
y = double(points.energy_mean_meV(:));
sigma = max(double(points.energy_err_meV(:)), sigma_floor_meV);
w = 1 ./ (sigma .^ 2);

slope = sum(w .* x .* y) ./ sum(w .* x .^ 2);
residuals = y - slope .* x;
n = numel(y);
dof = max(n - 1, 1);
chi2 = sum((residuals ./ sigma) .^ 2);

fit = struct();
fit.slope_meV_sqrtA = slope;
fit.RMSE_meV = sqrt(mean(residuals .^ 2));
fit.reduced_chi2 = chi2 / dof;
fit.R_squared = local_r_squared(y, residuals);
fit.residuals_meV = residuals;
end


function r2 = local_r_squared(y, residuals)
ss_res = sum(residuals .^ 2);
ss_tot = sum((y - mean(y)) .^ 2);
r2 = 1 - ss_res / max(ss_tot, eps);
end


function local_plot_primary_fit(points, summary, primary_qmax_Ainv, ...
    fig_png, fig_pdf)
row = summary(abs(summary.qmax_Ainv - primary_qmax_Ainv) < 1e-12, :);
if isempty(row)
    row = summary(end, :);
end

mask = points.q_abs_Ainv <= row.qmax_Ainv(1);
sub = points(mask, :);
x_curve = linspace(0, max(points.sqrt_q_Ainv) * 1.05, 200)';
y_free = row.free_slope_meV_sqrtA(1) .* x_curve + ...
    row.free_intercept_meV(1);
y_zero = row.zero_slope_meV_sqrtA(1) .* x_curve;

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 760 520]);
ax = axes(fig);
hold(ax, 'on');

outside = points(~mask, :);
if ~isempty(outside)
    errorbar(ax, outside.sqrt_q_Ainv, outside.energy_mean_meV, ...
        outside.energy_err_meV, 'o', 'Color', [0.65 0.65 0.65], ...
        'MarkerFaceColor', [0.85 0.85 0.85], ...
        'LineWidth', 0.9, 'DisplayName', 'excluded higher q');
end

errorbar(ax, sub.sqrt_q_Ainv, sub.energy_mean_meV, ...
    sub.energy_err_meV, 'o', 'Color', [0.10 0.34 0.76], ...
    'MarkerFaceColor', [0.10 0.34 0.76], 'LineWidth', 1.0, ...
    'DisplayName', sprintf('20w B1, |q| <= %.3f A^{-1}', ...
    row.qmax_Ainv(1)));
plot(ax, x_curve, y_free, '-', 'Color', [0.86 0.22 0.18], ...
    'LineWidth', 1.8, 'DisplayName', ...
    sprintf('free intercept b = %.0f meV', row.free_intercept_meV(1)));
plot(ax, x_curve, y_zero, '--', 'Color', [0.10 0.55 0.24], ...
    'LineWidth', 1.6, 'DisplayName', 'forced through zero');

xlabel(ax, 'sqrt(|q|) (A^{-1/2})');
ylabel(ax, 'B1 energy (meV)');
title(ax, '20w B1 low-q sqrt zero-crossing diagnostic');
grid(ax, 'on');
box(ax, 'on');
legend(ax, 'Location', 'northwest');

exportgraphics(fig, fig_png, 'Resolution', 220);
exportgraphics(fig, fig_pdf, 'ContentType', 'vector');
close(fig);
end


function local_write_report(report_path, input_csv, summary, ...
    primary_qmax_Ainv, fig_png)
row = summary(abs(summary.qmax_Ainv - primary_qmax_Ainv) < 1e-12, :);
if isempty(row)
    row = summary(end, :);
end

fid = fopen(report_path, 'w');
if fid < 0
    error('run_20w_lowq_sqrt_zero_diagnostic:ReportOpenFailed', ...
        'Could not open report for writing: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# 20w B1 Low-q Sqrt Zero Diagnostic\n\n');
fprintf(fid, '- Source table: `%s`\n', input_csv);
fprintf(fid, '- Session: `no_PL2_20w_2film`\n');
fprintf(fid, '- Model comparison: `E = a sqrt(|q|) + b` versus `E = a sqrt(|q|)`.\n');
fprintf(fid, '- Primary low-q window: `|q| <= %.4f A^-1`.\n\n', ...
    row.qmax_Ainv(1));

fprintf(fid, '## Primary Result\n\n');
fprintf(fid, '- Points: %d, q range %.4f-%.4f A^-1.\n', ...
    row.n_points(1), row.q_abs_min_Ainv(1), row.q_abs_max_data_Ainv(1));
fprintf(fid, '- Free intercept: %.2f meV, 95%% CI [%.2f, %.2f] meV.\n', ...
    row.free_intercept_meV(1), ...
    row.free_intercept_ci95_low_meV(1), ...
    row.free_intercept_ci95_high_meV(1));
fprintf(fid, '- RMSE: free %.2f meV, forced-zero %.2f meV.\n', ...
    row.free_RMSE_meV(1), row.zero_RMSE_meV(1));
fprintf(fid, '- Interpretation: in the very-low-q 20w window, the bare sqrt line still prefers a positive intercept; use this as a diagnostic, not as a replacement for the screened quasi-2D model.\n\n');

[~, fig_name, fig_ext] = fileparts(fig_png);
fprintf(fid, '![20w low-q sqrt zero diagnostic](%s%s)\n', ...
    fig_name, fig_ext);
end
