function output = run_20w_lowq_gap_diagnostic()
%RUN_20W_LOWQ_GAP_DIAGNOSTIC Fit finite-q=0 energy diagnostics for 20w B1.
%
% This script tests whether the 20w 2film low-q B1 points require a finite
% q->0 energy scale. It is diagnostic only and does not change the principal
% B1 quasi-2D plasmon fit.

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
    '20w_lowq_gap_model_260508');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

session_key = "no_PL2_20w_2film";
qmax_values_Ainv = [0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.0475];
primary_qmax_Ainv = 0.040;
sigma_floor_meV = 10;
epsilon_bg = 4.5;

points = local_load_20w_points(input_csv, session_key, ...
    max(qmax_values_Ainv), sigma_floor_meV);
summary = local_gap_model_sweep(points, qmax_values_Ainv, ...
    epsilon_bg, sigma_floor_meV);

points_csv = fullfile(output_dir, '20w_lowq_gap_model_points.csv');
summary_csv = fullfile(output_dir, '20w_lowq_gap_model_summary.csv');
writetable(points, points_csv);
writetable(summary, summary_csv);

fig_png = fullfile(output_dir, '20w_lowq_gap_model_diagnostic.png');
fig_pdf = fullfile(output_dir, '20w_lowq_gap_model_diagnostic.pdf');
local_plot_primary_models(points, summary, primary_qmax_Ainv, ...
    fig_png, fig_pdf);

report_path = fullfile(output_dir, '20w_lowq_gap_model_report.md');
local_write_report(report_path, input_csv, summary, primary_qmax_Ainv, ...
    epsilon_bg, fig_png);

mat_path = fullfile(output_dir, '20w_lowq_gap_model_results.mat');
save(mat_path, 'points', 'summary', 'session_key', 'qmax_values_Ainv', ...
    'primary_qmax_Ainv', 'sigma_floor_meV', 'epsilon_bg');

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

fprintf('20w low-q gap diagnostic complete.\n');
fprintf('  Output: %s\n', output_dir);
fprintf('  Primary qmax: %.4f A^-1, epsilon_bg: %.3g\n', ...
    primary_qmax_Ainv, epsilon_bg);
end


function points = local_load_20w_points(input_csv, session_key, ...
    qmax_Ainv, sigma_floor_meV)
if ~isfile(input_csv)
    error('run_20w_lowq_gap_diagnostic:MissingInput', ...
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
points.energy2_meV2 = points.energy_mean_meV .^ 2;
points.sigma_used_meV = max(points.energy_err_meV, sigma_floor_meV);
points.energy2_sigma_meV2 = 2 .* points.energy_mean_meV .* ...
    points.sigma_used_meV;

if height(points) < 3
    error('run_20w_lowq_gap_diagnostic:InsufficientPoints', ...
        'Need at least 3 low-q 20w B1 points.');
end
end


function local_require_columns(tbl, required)
names = tbl.Properties.VariableNames;
for i = 1:numel(required)
    if ~ismember(required{i}, names)
        error('run_20w_lowq_gap_diagnostic:MissingColumn', ...
            'Input table is missing required column "%s".', required{i});
    end
end
end


function summary = local_gap_model_sweep(points, qmax_values_Ainv, ...
    epsilon_bg, sigma_floor_meV)
rows = table();
for i = 1:numel(qmax_values_Ainv)
    qmax = qmax_values_Ainv(i);
    sub = points(points.q_abs_Ainv <= qmax, :);
    if height(sub) < 3
        continue;
    end

    fits = { ...
        local_fit_linear_E2_gap(sub, sigma_floor_meV), ...
        local_fit_linear_E2_zero_gap(sub, sigma_floor_meV), ...
        local_fit_screened_gap(sub, epsilon_bg, sigma_floor_meV), ...
        local_fit_screened_zero_gap(sub, epsilon_bg, sigma_floor_meV)};

    for j = 1:numel(fits)
        rows = [rows; local_fit_summary_row(fits{j}, sub, qmax)]; %#ok<AGROW>
    end
end
summary = rows;
end


function fit = local_fit_linear_E2_gap(points, sigma_floor_meV)
q = double(points.q_abs_Ainv(:));
E = double(points.energy_mean_meV(:));
sigma_E = max(double(points.energy_err_meV(:)), sigma_floor_meV);
y = E .^ 2;
sigma_y = 2 .* E .* sigma_E;
w = 1 ./ max(sigma_y, eps) .^ 2;

X = [ones(size(q)), q];
W = diag(w);
beta = (X' * W * X) \ (X' * W * y);
Delta2 = max(beta(1), 0);
C = max(beta(2), 0);

resid_y = y - X * beta;
dof_y = max(numel(y) - 2, 1);
red_chi2_y = sum((resid_y ./ sigma_y) .^ 2) / dof_y;
cov_beta = pinv(X' * W * X) * red_chi2_y;
Delta2_sigma = sqrt(max(cov_beta(1, 1), 0));

E_pred = sqrt(max(Delta2 + C .* q, 0));
fit = local_build_fit_struct('linear_E2_gap', q, E, sigma_E, E_pred, 2);
fit.Delta2_meV2 = Delta2;
fit.Delta2_sigma_meV2 = Delta2_sigma;
fit.Delta_meV = sqrt(Delta2);
fit.Delta_ci95_low_meV = sqrt(max(Delta2 - 1.96 * Delta2_sigma, 0));
fit.Delta_ci95_high_meV = sqrt(max(Delta2 + 1.96 * Delta2_sigma, 0));
fit.C_meV2_A = C;
fit.A_meV2_A = NaN;
fit.rho_A = NaN;
fit.epsilon_bg = NaN;
end


function fit = local_fit_linear_E2_zero_gap(points, sigma_floor_meV)
q = double(points.q_abs_Ainv(:));
E = double(points.energy_mean_meV(:));
sigma_E = max(double(points.energy_err_meV(:)), sigma_floor_meV);
y = E .^ 2;
sigma_y = 2 .* E .* sigma_E;
w = 1 ./ max(sigma_y, eps) .^ 2;

C = sum(w .* q .* y) ./ sum(w .* q .^ 2);
C = max(C, 0);
E_pred = sqrt(max(C .* q, 0));

fit = local_build_fit_struct('linear_E2_zero_gap', q, E, sigma_E, E_pred, 1);
fit.Delta2_meV2 = 0;
fit.Delta2_sigma_meV2 = NaN;
fit.Delta_meV = 0;
fit.Delta_ci95_low_meV = 0;
fit.Delta_ci95_high_meV = 0;
fit.C_meV2_A = C;
fit.A_meV2_A = NaN;
fit.rho_A = NaN;
fit.epsilon_bg = NaN;
end


function fit = local_fit_screened_gap(points, epsilon_bg, sigma_floor_meV)
q = double(points.q_abs_Ainv(:));
E = double(points.energy_mean_meV(:));
sigma_E = max(double(points.energy_err_meV(:)), sigma_floor_meV);

lin = local_fit_linear_E2_gap(points, sigma_floor_meV);
Delta0 = min(max(lin.Delta_meV, 0), 0.95 * min(E));
A0 = max(lin.C_meV2_A * epsilon_bg, 1);
rho0 = 30;

model = @(p, x) sqrt(max(p(1) .^ 2 + ...
    p(2) .* x ./ (epsilon_bg + p(3) .* x), 0));
p = local_lsq_fit(@(p, x) model(p, x), [Delta0, A0, rho0], ...
    [0, 0, 0.05], [max(E), Inf, 5000], q, E, sigma_E);
E_pred = model(p, q);

fit = local_build_fit_struct('screened_gap', q, E, sigma_E, E_pred, 3);
fit.Delta_meV = p(1);
fit.Delta2_meV2 = p(1) .^ 2;
fit.Delta2_sigma_meV2 = NaN;
fit.Delta_ci95_low_meV = NaN;
fit.Delta_ci95_high_meV = NaN;
fit.C_meV2_A = NaN;
fit.A_meV2_A = p(2);
fit.rho_A = p(3);
fit.epsilon_bg = epsilon_bg;
end


function fit = local_fit_screened_zero_gap(points, epsilon_bg, sigma_floor_meV)
q = double(points.q_abs_Ainv(:));
E = double(points.energy_mean_meV(:));
sigma_E = max(double(points.energy_err_meV(:)), sigma_floor_meV);

lin = local_fit_linear_E2_zero_gap(points, sigma_floor_meV);
A0 = max(lin.C_meV2_A * epsilon_bg, 1);
rho0 = 30;

model = @(p, x) sqrt(max(p(1) .* x ./ ...
    (epsilon_bg + p(2) .* x), 0));
p = local_lsq_fit(@(p, x) model(p, x), [A0, rho0], ...
    [0, 0.05], [Inf, 5000], q, E, sigma_E);
E_pred = model(p, q);

fit = local_build_fit_struct('screened_zero_gap', q, E, sigma_E, E_pred, 2);
fit.Delta_meV = 0;
fit.Delta2_meV2 = 0;
fit.Delta2_sigma_meV2 = NaN;
fit.Delta_ci95_low_meV = 0;
fit.Delta_ci95_high_meV = 0;
fit.C_meV2_A = NaN;
fit.A_meV2_A = p(1);
fit.rho_A = p(2);
fit.epsilon_bg = epsilon_bg;
end


function p = local_lsq_fit(model, p0, lb, ub, q, E, sigma_E)
try
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxFunctionEvaluations', 20000, 'MaxIterations', 3000, ...
        'FunctionTolerance', 1e-12, 'StepTolerance', 1e-12);
    weighted_model = @(p_in, q_in) model(p_in, q_in) ./ sigma_E;
    p = lsqcurvefit(weighted_model, p0, q, E ./ sigma_E, lb, ub, opts);
catch
    cost = @(p_in) sum(((model(local_apply_bounds(p_in, lb, ub), q) - E) ...
        ./ sigma_E) .^ 2);
    fmin_opts = optimset('Display', 'off', 'MaxFunEvals', 20000, ...
        'MaxIter', 3000, 'TolFun', 1e-12, 'TolX', 1e-12);
    p = fminsearch(cost, p0, fmin_opts);
    p = local_apply_bounds(p, lb, ub);
end
end


function p = local_apply_bounds(p, lb, ub)
p = max(p, lb);
finite_ub = isfinite(ub);
p(finite_ub) = min(p(finite_ub), ub(finite_ub));
end


function fit = local_build_fit_struct(model_key, q, E, sigma_E, E_pred, n_params)
residuals = E - E_pred;
n = numel(E);
dof = max(n - n_params, 1);
chi2 = sum((residuals ./ sigma_E) .^ 2);

fit = struct();
fit.model_key = model_key;
fit.n_points = n;
fit.n_params = n_params;
fit.dof = dof;
fit.chi2 = chi2;
fit.reduced_chi2 = chi2 / dof;
fit.RMSE_meV = sqrt(mean(residuals .^ 2));
fit.R_squared = local_r_squared(E, residuals);
fit.AIC = chi2 + 2 * n_params;
fit.BIC = chi2 + n_params * log(n);
fit.residuals_meV = residuals;
fit.E_pred_meV = E_pred;
fit.q_abs_min_Ainv = min(q);
fit.q_abs_max_data_Ainv = max(q);
end


function r2 = local_r_squared(y, residuals)
ss_res = sum(residuals .^ 2);
ss_tot = sum((y - mean(y)) .^ 2);
r2 = 1 - ss_res / max(ss_tot, eps);
end


function row = local_fit_summary_row(fit, points, qmax_Ainv)
row = table( ...
    {fit.model_key}, qmax_Ainv, fit.n_points, fit.n_params, fit.dof, ...
    fit.q_abs_min_Ainv, fit.q_abs_max_data_Ainv, ...
    fit.Delta_meV, fit.Delta_ci95_low_meV, fit.Delta_ci95_high_meV, ...
    fit.Delta2_meV2, fit.Delta2_sigma_meV2, fit.C_meV2_A, ...
    fit.A_meV2_A, fit.rho_A, fit.epsilon_bg, ...
    fit.RMSE_meV, fit.reduced_chi2, fit.R_squared, fit.AIC, fit.BIC, ...
    'VariableNames', {'model_key', 'qmax_Ainv', 'n_points', ...
    'n_params', 'dof', 'q_abs_min_Ainv', 'q_abs_max_data_Ainv', ...
    'Delta_meV', 'Delta_ci95_low_meV', 'Delta_ci95_high_meV', ...
    'Delta2_meV2', 'Delta2_sigma_meV2', 'C_meV2_A', ...
    'A_meV2_A', 'rho_A', 'epsilon_bg', 'RMSE_meV', ...
    'reduced_chi2', 'R_squared', 'AIC', 'BIC'});

if height(points) ~= fit.n_points
    error('run_20w_lowq_gap_diagnostic:InternalSizeMismatch', ...
        'Summary row point count does not match fit point count.');
end
end


function local_plot_primary_models(points, summary, primary_qmax_Ainv, ...
    fig_png, fig_pdf)
mask = points.q_abs_Ainv <= primary_qmax_Ainv;
sub = points(mask, :);
plot_rows = summary(abs(summary.qmax_Ainv - primary_qmax_Ainv) < 1e-12, :);

q_curve = linspace(0, max(points.q_abs_Ainv) * 1.05, 220)';

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1050 480]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile;
hold(ax1, 'on');
outside = points(~mask, :);
if ~isempty(outside)
    errorbar(ax1, outside.q_abs_Ainv, outside.energy_mean_meV, ...
        outside.energy_err_meV, 'o', 'Color', [0.65 0.65 0.65], ...
        'MarkerFaceColor', [0.85 0.85 0.85], 'LineWidth', 0.9, ...
        'DisplayName', 'excluded higher q');
end
errorbar(ax1, sub.q_abs_Ainv, sub.energy_mean_meV, sub.energy_err_meV, ...
    'o', 'Color', [0.10 0.34 0.76], 'MarkerFaceColor', [0.10 0.34 0.76], ...
    'LineWidth', 1.0, 'DisplayName', '20w B1 low q');
local_plot_model_curve(ax1, q_curve, plot_rows, 'linear_E2_gap', ...
    [0.86 0.22 0.18], '-', 'E2 gap');
local_plot_model_curve(ax1, q_curve, plot_rows, 'linear_E2_zero_gap', ...
    [0.10 0.55 0.24], '--', 'E2 zero gap');
local_plot_model_curve(ax1, q_curve, plot_rows, 'screened_gap', ...
    [0.48 0.18 0.72], '-', 'screened gap');
local_plot_model_curve(ax1, q_curve, plot_rows, 'screened_zero_gap', ...
    [0.05 0.55 0.75], '--', 'screened zero gap');
xlabel(ax1, '|q| (A^{-1})');
ylabel(ax1, 'B1 energy (meV)');
title(ax1, sprintf('20w 2film low-q models, |q| <= %.3f A^{-1}', ...
    primary_qmax_Ainv));
grid(ax1, 'on');
box(ax1, 'on');
legend(ax1, 'Location', 'northwest');

ax2 = nexttile;
hold(ax2, 'on');
errorbar(ax2, sub.q_abs_Ainv, sub.energy2_meV2 ./ 1e6, ...
    sub.energy2_sigma_meV2 ./ 1e6, 'o', ...
    'Color', [0.10 0.34 0.76], 'MarkerFaceColor', [0.10 0.34 0.76], ...
    'LineWidth', 1.0, 'DisplayName', '20w B1 low q');
local_plot_E2_curve(ax2, q_curve, plot_rows, 'linear_E2_gap', ...
    [0.86 0.22 0.18], '-', 'E2 gap');
local_plot_E2_curve(ax2, q_curve, plot_rows, 'linear_E2_zero_gap', ...
    [0.10 0.55 0.24], '--', 'E2 zero gap');
xlabel(ax2, '|q| (A^{-1})');
ylabel(ax2, 'B1 energy^2 (10^6 meV^2)');
title(ax2, 'Linearized gap diagnostic');
grid(ax2, 'on');
box(ax2, 'on');
legend(ax2, 'Location', 'northwest');

exportgraphics(fig, fig_png, 'Resolution', 220);
exportgraphics(fig, fig_pdf, 'ContentType', 'vector');
close(fig);
end


function local_plot_model_curve(ax, q_curve, rows, model_key, color, style, label)
row = rows(strcmp(rows.model_key, model_key), :);
if isempty(row)
    return;
end
E = local_predict_from_row(row, q_curve);
plot(ax, q_curve, E, style, 'Color', color, 'LineWidth', 1.7, ...
    'DisplayName', label);
end


function local_plot_E2_curve(ax, q_curve, rows, model_key, color, style, label)
row = rows(strcmp(rows.model_key, model_key), :);
if isempty(row)
    return;
end
E = local_predict_from_row(row, q_curve);
plot(ax, q_curve, (E .^ 2) ./ 1e6, style, 'Color', color, ...
    'LineWidth', 1.7, 'DisplayName', label);
end


function E = local_predict_from_row(row, q)
key = string(row.model_key{1});
switch key
    case "linear_E2_gap"
        E = sqrt(max(row.Delta2_meV2(1) + row.C_meV2_A(1) .* q, 0));
    case "linear_E2_zero_gap"
        E = sqrt(max(row.C_meV2_A(1) .* q, 0));
    case "screened_gap"
        E = sqrt(max(row.Delta_meV(1) .^ 2 + row.A_meV2_A(1) .* q ./ ...
            (row.epsilon_bg(1) + row.rho_A(1) .* q), 0));
    case "screened_zero_gap"
        E = sqrt(max(row.A_meV2_A(1) .* q ./ ...
            (row.epsilon_bg(1) + row.rho_A(1) .* q), 0));
    otherwise
        E = NaN(size(q));
end
end


function local_write_report(report_path, input_csv, summary, ...
    primary_qmax_Ainv, epsilon_bg, fig_png)
rows = summary(abs(summary.qmax_Ainv - primary_qmax_Ainv) < 1e-12, :);
lin_gap = rows(strcmp(rows.model_key, 'linear_E2_gap'), :);
lin_zero = rows(strcmp(rows.model_key, 'linear_E2_zero_gap'), :);
scr_gap = rows(strcmp(rows.model_key, 'screened_gap'), :);
scr_zero = rows(strcmp(rows.model_key, 'screened_zero_gap'), :);

fid = fopen(report_path, 'w');
if fid < 0
    error('run_20w_lowq_gap_diagnostic:ReportOpenFailed', ...
        'Could not open report for writing: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# 20w B1 Low-q Gap Diagnostic\n\n');
fprintf(fid, '- Source table: `%s`\n', input_csv);
fprintf(fid, '- Session: `no_PL2_20w_2film`\n');
fprintf(fid, '- Primary low-q window: `|q| <= %.4f A^-1`.\n', ...
    primary_qmax_Ainv);
fprintf(fid, '- Screened-model epsilon_bg: %.3g.\n\n', epsilon_bg);

fprintf(fid, '## Primary Result\n\n');
fprintf(fid, '- Linear `E^2 = Delta^2 + C q`: Delta = %.1f meV, 95%% CI [%.1f, %.1f] meV, RMSE = %.2f meV.\n', ...
    lin_gap.Delta_meV(1), lin_gap.Delta_ci95_low_meV(1), ...
    lin_gap.Delta_ci95_high_meV(1), lin_gap.RMSE_meV(1));
fprintf(fid, '- Linear zero-gap `E^2 = C q`: RMSE = %.2f meV.\n', ...
    lin_zero.RMSE_meV(1));
fprintf(fid, '- Screened gap: Delta = %.1f meV, rho = %.2f A, RMSE = %.2f meV, AIC = %.2f.\n', ...
    scr_gap.Delta_meV(1), scr_gap.rho_A(1), ...
    scr_gap.RMSE_meV(1), scr_gap.AIC(1));
fprintf(fid, '- Screened zero-gap: rho = %.2f A, RMSE = %.2f meV, AIC = %.2f.\n\n', ...
    scr_zero.rho_A(1), scr_zero.RMSE_meV(1), scr_zero.AIC(1));

fprintf(fid, 'Interpretation: the bare zero-gap sqrt model is disfavored for the 20w 2film low-q points. A finite-energy linearized model fits the low-q curvature very well, but a screened zero-gap quasi-2D form also absorbs much of the same curvature. Treat Delta as an effective low-q scale, not yet as a confirmed true gap.\n\n');

[~, fig_name, fig_ext] = fileparts(fig_png);
fprintf(fid, '![20w low-q gap diagnostic](%s%s)\n', fig_name, fig_ext);
end
