function output = run_lowq_powerlaw_deviation_export()
%RUN_LOWQ_POWERLAW_DEVIATION_EXPORT Export B1 sqrt-law deviation figures.
%
% The diagnostic plots E versus sqrt(|q|). A pure low-q 2D plasmon power
% law should be approximately linear and pass through zero in this view.

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
    'lowq_powerlaw_deviation_260508');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

sessions = local_session_config();
primary_qmax_Ainv = 0.040;
display_qmax_Ainv = 0.050;
sigma_floor_meV = 10;

points = local_load_points(input_csv, sessions, display_qmax_Ainv, ...
    sigma_floor_meV);
summary = local_fit_summary(points, sessions, primary_qmax_Ainv, ...
    sigma_floor_meV);
residuals = local_residual_table(points, summary);

points_csv = fullfile(output_dir, 'powerlaw_deviation_points.csv');
summary_csv = fullfile(output_dir, 'powerlaw_deviation_summary.csv');
residual_csv = fullfile(output_dir, 'powerlaw_deviation_residuals.csv');
writetable(points, points_csv);
writetable(summary, summary_csv);
writetable(residuals, residual_csv);

sqrt_png = fullfile(output_dir, 'powerlaw_deviation_sqrt_comparison.png');
sqrt_pdf = fullfile(output_dir, 'powerlaw_deviation_sqrt_comparison.pdf');
resid_png = fullfile(output_dir, 'powerlaw_deviation_residuals.png');
resid_pdf = fullfile(output_dir, 'powerlaw_deviation_residuals.pdf');

local_plot_sqrt_comparison(points, summary, sessions, primary_qmax_Ainv, ...
    sqrt_png, sqrt_pdf);
local_plot_residuals(residuals, sessions, primary_qmax_Ainv, ...
    resid_png, resid_pdf);

report_path = fullfile(output_dir, 'powerlaw_deviation_report.md');
local_write_report(report_path, input_csv, summary, ...
    primary_qmax_Ainv, sqrt_png, resid_png);

mat_path = fullfile(output_dir, 'powerlaw_deviation_results.mat');
save(mat_path, 'points', 'summary', 'residuals', 'sessions', ...
    'primary_qmax_Ainv', 'display_qmax_Ainv', 'sigma_floor_meV');

output = struct();
output.output_dir = output_dir;
output.points_csv = points_csv;
output.summary_csv = summary_csv;
output.residual_csv = residual_csv;
output.sqrt_png = sqrt_png;
output.sqrt_pdf = sqrt_pdf;
output.resid_png = resid_png;
output.resid_pdf = resid_pdf;
output.report_path = report_path;
output.mat_path = mat_path;
output.points = points;
output.summary = summary;
output.residuals = residuals;

fprintf('Low-q power-law deviation export complete.\n');
fprintf('  Output: %s\n', output_dir);
fprintf('  Primary qmax: %.4f A^-1\n', primary_qmax_Ainv);
end


function sessions = local_session_config()
sessions = struct( ...
    'session_key', {}, ...
    'session_label', {}, ...
    'thickness_label', {}, ...
    'color', {});

sessions(end + 1) = struct( ...
    'session_key', "590_PL2_10w", ...
    'session_label', "590 10w defocus", ...
    'thickness_label', "1film", ...
    'color', [0.120, 0.470, 0.900]);

sessions(end + 1) = struct( ...
    'session_key', "n0_PL2_10w_repeat", ...
    'session_label', "n0 10w repeat", ...
    'thickness_label', "1film", ...
    'color', [0.160, 0.500, 0.220]);

sessions(end + 1) = struct( ...
    'session_key', "no_PL2_20w_2film", ...
    'session_label', "20w", ...
    'thickness_label', "2film", ...
    'color', [0.930, 0.280, 0.300]);
end


function points = local_load_points(input_csv, sessions, qmax_Ainv, ...
    sigma_floor_meV)
if ~isfile(input_csv)
    error('run_lowq_powerlaw_deviation_export:MissingInput', ...
        'Missing required input CSV: %s', input_csv);
end

tbl = readtable(input_csv);
required = {'session_key', 'q_abs_Ainv', 'energy_mean_meV', ...
    'energy_err_meV', 'R2_mean_all', 'include_for_fit'};
local_require_columns(tbl, required);

points = table();
for i = 1:numel(sessions)
    key = sessions(i).session_key;
    mask = strcmp(string(tbl.session_key), key) & ...
        tbl.include_for_fit == 1 & ...
        isfinite(tbl.q_abs_Ainv) & tbl.q_abs_Ainv > 0 & ...
        tbl.q_abs_Ainv <= qmax_Ainv & ...
        isfinite(tbl.energy_mean_meV) & tbl.energy_mean_meV > 0;
    sub = tbl(mask, :);
    if height(sub) < 3
        error('run_lowq_powerlaw_deviation_export:InsufficientPoints', ...
            'Need at least 3 low-q B1 points for %s.', key);
    end

    sub.session_order = repmat(i, height(sub), 1);
    sub.thickness_label = repmat(sessions(i).thickness_label, height(sub), 1);
    sub.sqrt_q_Ainv = sqrt(sub.q_abs_Ainv);
    sub.sigma_used_meV = max(sub.energy_err_meV, sigma_floor_meV);
    points = [points; sub]; %#ok<AGROW>
end

points = sortrows(points, {'session_order', 'q_abs_Ainv'});
end


function local_require_columns(tbl, required)
names = tbl.Properties.VariableNames;
for i = 1:numel(required)
    if ~ismember(required{i}, names)
        error('run_lowq_powerlaw_deviation_export:MissingColumn', ...
            'Input table is missing required column "%s".', required{i});
    end
end
end


function summary = local_fit_summary(points, sessions, primary_qmax_Ainv, ...
    sigma_floor_meV)
rows = table();
for i = 1:numel(sessions)
    sub = points(points.session_order == i & ...
        points.q_abs_Ainv <= primary_qmax_Ainv, :);
    free_fit = local_fit_sqrt_free_intercept(sub, sigma_floor_meV);
    zero_fit = local_fit_sqrt_zero_intercept(sub, sigma_floor_meV);

    row = table( ...
        i, sessions(i).session_key, sessions(i).session_label, ...
        sessions(i).thickness_label, primary_qmax_Ainv, height(sub), ...
        min(sub.q_abs_Ainv), max(sub.q_abs_Ainv), ...
        free_fit.slope_meV_sqrtA, free_fit.free_intercept_meV, ...
        free_fit.free_intercept_ci95_low_meV, ...
        free_fit.free_intercept_ci95_high_meV, ...
        free_fit.free_intercept_sigma_meV, free_fit.free_RMSE_meV, ...
        free_fit.free_R_squared, zero_fit.zero_slope_meV_sqrtA, ...
        zero_fit.zero_RMSE_meV, zero_fit.zero_R_squared, ...
        zero_fit.zero_RMSE_meV - free_fit.free_RMSE_meV, ...
        'VariableNames', {'session_order', 'session_key', ...
        'session_label', 'thickness_label', 'qmax_Ainv', 'n_points', ...
        'q_abs_min_Ainv', 'q_abs_max_data_Ainv', ...
        'free_slope_meV_sqrtA', 'free_intercept_meV', ...
        'free_intercept_ci95_low_meV', ...
        'free_intercept_ci95_high_meV', ...
        'free_intercept_sigma_meV', 'free_RMSE_meV', ...
        'free_R_squared', 'zero_slope_meV_sqrtA', ...
        'zero_RMSE_meV', 'zero_R_squared', ...
        'delta_RMSE_zero_minus_free_meV'});
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
dof = max(numel(y) - 2, 1);
chi2 = sum((residuals ./ sigma) .^ 2);
cov_beta = pinv(X' * W * X) * (chi2 / dof);
intercept_sigma = sqrt(max(cov_beta(2, 2), 0));

fit = struct();
fit.slope_meV_sqrtA = beta(1);
fit.free_intercept_meV = beta(2);
fit.free_intercept_sigma_meV = intercept_sigma;
fit.free_intercept_ci95_low_meV = beta(2) - 1.96 * intercept_sigma;
fit.free_intercept_ci95_high_meV = beta(2) + 1.96 * intercept_sigma;
fit.free_RMSE_meV = sqrt(mean(residuals .^ 2));
fit.free_R_squared = local_r_squared(y, residuals);
end


function fit = local_fit_sqrt_zero_intercept(points, sigma_floor_meV)
x = double(points.sqrt_q_Ainv(:));
y = double(points.energy_mean_meV(:));
sigma = max(double(points.energy_err_meV(:)), sigma_floor_meV);
w = 1 ./ (sigma .^ 2);

slope = sum(w .* x .* y) ./ sum(w .* x .^ 2);
residuals = y - slope .* x;

fit = struct();
fit.zero_slope_meV_sqrtA = slope;
fit.zero_RMSE_meV = sqrt(mean(residuals .^ 2));
fit.zero_R_squared = local_r_squared(y, residuals);
end


function r2 = local_r_squared(y, residuals)
ss_res = sum(residuals .^ 2);
ss_tot = sum((y - mean(y)) .^ 2);
r2 = 1 - ss_res / max(ss_tot, eps);
end


function residuals = local_residual_table(points, summary)
residuals = table();
for i = 1:height(summary)
    sub = points(points.session_order == summary.session_order(i), :);
    zero_pred = summary.zero_slope_meV_sqrtA(i) .* sub.sqrt_q_Ainv;
    free_pred = summary.free_slope_meV_sqrtA(i) .* sub.sqrt_q_Ainv + ...
        summary.free_intercept_meV(i);
    part = table( ...
        sub.session_order, sub.session_key, sub.thickness_label, ...
        sub.q_abs_Ainv, sub.sqrt_q_Ainv, sub.energy_mean_meV, ...
        sub.energy_err_meV, zero_pred, sub.energy_mean_meV - zero_pred, ...
        free_pred, sub.energy_mean_meV - free_pred, ...
        'VariableNames', {'session_order', 'session_key', ...
        'thickness_label', 'q_abs_Ainv', 'sqrt_q_Ainv', ...
        'energy_mean_meV', 'energy_err_meV', 'zero_pred_meV', ...
        'zero_residual_meV', 'free_pred_meV', 'free_residual_meV'});
    residuals = [residuals; part]; %#ok<AGROW>
end
end


function local_plot_sqrt_comparison(points, summary, sessions, ...
    primary_qmax_Ainv, fig_png, fig_pdf)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1280 520]);
tiledlayout(fig, 1, numel(sessions), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for i = 1:numel(sessions)
    ax = nexttile;
    hold(ax, 'on');
    row = summary(summary.session_order == i, :);
    sub_all = points(points.session_order == i, :);
    in_fit = sub_all.q_abs_Ainv <= primary_qmax_Ainv;
    sub = sub_all(in_fit, :);
    outside = sub_all(~in_fit, :);
    x_curve = linspace(0, max(sub_all.sqrt_q_Ainv) * 1.06, 220)';
    y_free = row.free_slope_meV_sqrtA(1) .* x_curve + ...
        row.free_intercept_meV(1);
    y_zero = row.zero_slope_meV_sqrtA(1) .* x_curve;

    if ~isempty(outside)
        errorbar(ax, outside.sqrt_q_Ainv, outside.energy_mean_meV, ...
            outside.energy_err_meV, 'o', 'Color', [0.65 0.65 0.65], ...
            'MarkerFaceColor', [0.85 0.85 0.85], ...
            'LineWidth', 0.8, 'DisplayName', 'excluded higher q');
    end
    errorbar(ax, sub.sqrt_q_Ainv, sub.energy_mean_meV, ...
        sub.energy_err_meV, 'o', 'Color', sessions(i).color, ...
        'MarkerFaceColor', sessions(i).color, 'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s B1', sessions(i).thickness_label));
    plot(ax, x_curve, y_free, '-', 'Color', [0.86 0.22 0.18], ...
        'LineWidth', 1.8, 'DisplayName', ...
        sprintf('free b = %.0f meV', row.free_intercept_meV(1)));
    plot(ax, x_curve, y_zero, '--', 'Color', [0.10 0.55 0.24], ...
        'LineWidth', 1.6, 'DisplayName', 'forced b = 0');

    title(ax, sprintf('%s %s', sessions(i).session_label, ...
        sessions(i).thickness_label));
    xlabel(ax, 'sqrt(|q|) (A^{-1/2})');
    if i == 1
        ylabel(ax, 'B1 energy (meV)');
    end
    grid(ax, 'on');
    box(ax, 'on');
    legend(ax, 'Location', 'northwest');
end

sgtitle(fig, sprintf('B1 deviation from E proportional to sqrt(q), |q| <= %.3f A^{-1}', ...
    primary_qmax_Ainv));
exportgraphics(fig, fig_png, 'Resolution', 220);
exportgraphics(fig, fig_pdf, 'ContentType', 'vector');
close(fig);
end


function local_plot_residuals(residuals, sessions, primary_qmax_Ainv, ...
    fig_png, fig_pdf)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100 100 1280 430]);
tiledlayout(fig, 1, numel(sessions), 'TileSpacing', 'compact', ...
    'Padding', 'compact');

fit_residuals = residuals(residuals.q_abs_Ainv <= primary_qmax_Ainv, :);
resid_extent = max(abs(fit_residuals.zero_residual_meV) + ...
    fit_residuals.energy_err_meV, [], 'omitnan');
ylim_abs = max(40, ceil(resid_extent / 10) * 10);

for i = 1:numel(sessions)
    ax = nexttile;
    hold(ax, 'on');
    sub = residuals(residuals.session_order == i, :);
    in_fit = sub.q_abs_Ainv <= primary_qmax_Ainv;
    plot(ax, [0, max(sub.sqrt_q_Ainv) * 1.05], [0, 0], '-', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 0.9, ...
        'DisplayName', 'zero residual');
    errorbar(ax, sub.sqrt_q_Ainv(in_fit), sub.zero_residual_meV(in_fit), ...
        sub.energy_err_meV(in_fit), 'o', 'Color', sessions(i).color, ...
        'MarkerFaceColor', sessions(i).color, 'LineWidth', 1.0, ...
        'DisplayName', 'residual to b=0');

    title(ax, sprintf('%s %s', sessions(i).session_label, ...
        sessions(i).thickness_label));
    xlabel(ax, 'sqrt(|q|) (A^{-1/2})');
    if i == 1
        ylabel(ax, 'E - a sqrt(|q|) (meV)');
    end
    grid(ax, 'on');
    box(ax, 'on');
    ylim(ax, [-ylim_abs, ylim_abs]);
    legend(ax, 'Location', 'best');
end

sgtitle(fig, 'Residual from forced-zero sqrt power law');
exportgraphics(fig, fig_png, 'Resolution', 220);
exportgraphics(fig, fig_pdf, 'ContentType', 'vector');
close(fig);
end


function local_write_report(report_path, input_csv, summary, ...
    primary_qmax_Ainv, sqrt_png, resid_png)
fid = fopen(report_path, 'w');
if fid < 0
    error('run_lowq_powerlaw_deviation_export:ReportOpenFailed', ...
        'Could not open report for writing: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# B1 Low-q Power-law Deviation\n\n');
fprintf(fid, '- Source table: `%s`\n', input_csv);
fprintf(fid, '- Primary fit window: `|q| <= %.4f A^-1`.\n', ...
    primary_qmax_Ainv);
fprintf(fid, '- Power-law check: `E = a sqrt(|q|)` versus `E = a sqrt(|q|) + b`.\n\n');

fprintf(fid, '| session | thickness | n | b free (meV) | RMSE b=0 | RMSE free |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|\n');
for i = 1:height(summary)
    fprintf(fid, '| %s | %s | %d | %.1f [%.1f, %.1f] | %.2f | %.2f |\n', ...
        summary.session_label(i), summary.thickness_label(i), ...
        summary.n_points(i), summary.free_intercept_meV(i), ...
        summary.free_intercept_ci95_low_meV(i), ...
        summary.free_intercept_ci95_high_meV(i), ...
        summary.zero_RMSE_meV(i), summary.free_RMSE_meV(i));
end

fprintf(fid, '\nThe 1film data stay close to the forced-zero sqrt power law. The 2film data prefer a positive free intercept, making the finite-thickness deviation visually explicit.\n\n');

[~, sqrt_name, sqrt_ext] = fileparts(sqrt_png);
[~, resid_name, resid_ext] = fileparts(resid_png);
fprintf(fid, '![sqrt comparison](%s%s)\n\n', sqrt_name, sqrt_ext);
fprintf(fid, '![sqrt residuals](%s%s)\n', resid_name, resid_ext);
end
