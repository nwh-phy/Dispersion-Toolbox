function output = run_epsilon_bg_sensitivity()
%RUN_EPSILON_BG_SENSITIVITY Refit B1/B3 with different dielectric backgrounds.
%
% This workflow uses exported branch point CSV files from the current
% Area-normalized Fano analysis. It directly fits
%
%   E(q) = sqrt(A |q| / (epsilon_bg + rho0 |q|))
%
% so epsilon_bg can represent an effective encapsulation environment.

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

results_root = fullfile(project_root, 'paper_results');
output_dir = fullfile(results_root, 'epsilon_bg_sensitivity_260506');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

epsilon_bg_values = [1, 4.5, 10, 15];
datasets = local_dataset_config(results_root);
branches = [1, 3];

rows = local_empty_rows();

for di = 1:numel(datasets)
    ds = datasets(di);
    for bi = 1:numel(branches)
        branch_id = branches(bi);
        csv_path = fullfile(ds.input_dir, sprintf('branch%d_points.csv', branch_id));
        if ~isfile(csv_path)
            warning('run_epsilon_bg_sensitivity:MissingBranchCsv', ...
                'Missing branch CSV: %s', csv_path);
            continue;
        end

        tbl = readtable(csv_path);
        [q_Ainv, energy_meV, confidence] = local_extract_branch_columns(tbl);

        for ei = 1:numel(epsilon_bg_values)
            epsilon_bg = epsilon_bg_values(ei);
            fit = local_fit_quasi2d_epsilon_bg(q_Ainv, energy_meV, ...
                confidence, epsilon_bg);

            rows = local_append_row(rows, ds, branch_id, epsilon_bg, fit);
        end
    end
end

summary = struct2table(rows);
summary = sortrows(summary, {'branch', 'session_key', 'epsilon_bg'});

summary_csv = fullfile(output_dir, 'epsilon_bg_sensitivity_summary.csv');
writetable(summary, summary_csv);

local_plot_metric(summary, output_dir, 'rho0_A', ...
    'rho0 (A)', 'rho0_vs_epsilon_bg.png');
local_plot_metric(summary, output_dir, 'q_c_Ainv', ...
    'q_c (1/A)', 'qc_vs_epsilon_bg.png');
local_plot_metric(summary, output_dir, 'E_flat_meV', ...
    'E_flat (meV)', 'eflat_vs_epsilon_bg.png');

report_path = fullfile(output_dir, 'epsilon_bg_sensitivity_report.md');
local_write_report(summary, report_path);

mat_path = fullfile(output_dir, 'epsilon_bg_sensitivity_results.mat');
save(mat_path, 'summary', 'epsilon_bg_values', 'datasets');

output = struct();
output.output_dir = output_dir;
output.summary_csv = summary_csv;
output.report_path = report_path;
output.summary = summary;

fprintf('Epsilon background sensitivity complete.\n');
fprintf('  Output: %s\n', output_dir);
fprintf('  Summary rows: %d\n', height(summary));
end


function datasets = local_dataset_config(results_root)
datasets = struct( ...
    'session_key', {}, ...
    'session_label', {}, ...
    'input_dir', {});

datasets(end + 1) = struct( ...
    'session_key', '590_PL2_10w', ...
    'session_label', '590 PL2 10w', ...
    'input_dir', fullfile(results_root, '590_gui_history_area_260506'));

datasets(end + 1) = struct( ...
    'session_key', 'n0_PL2_10w_repeat', ...
    'session_label', 'n0 PL2 10w repeat', ...
    'input_dir', fullfile(results_root, 'n0_PL2_10w_gui_history_area_260506'));

datasets(end + 1) = struct( ...
    'session_key', 'no_PL2_20w_2film', ...
    'session_label', 'no PL2 20w 2film', ...
    'input_dir', fullfile(results_root, ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined'));
end


function rows = local_empty_rows()
rows = struct( ...
    'session_key', {}, ...
    'session_label', {}, ...
    'branch', {}, ...
    'epsilon_bg', {}, ...
    'epsilon_s_equiv', {}, ...
    'n_points', {}, ...
    'q_abs_min_Ainv', {}, ...
    'q_abs_max_Ainv', {}, ...
    'A_fit', {}, ...
    'rho0_A', {}, ...
    'q_c_Ainv', {}, ...
    'E_flat_meV', {}, ...
    'R_squared', {}, ...
    'RMSE_meV', {}, ...
    'rho0_at_upper_bound', {});
end


function [q_Ainv, energy_meV, confidence] = local_extract_branch_columns(tbl)
required = {'q_Ainv', 'energy_meV'};
for i = 1:numel(required)
    if ~ismember(required{i}, tbl.Properties.VariableNames)
        error('run_epsilon_bg_sensitivity:MissingColumn', ...
            'Branch CSV is missing required column "%s".', required{i});
    end
end

q_Ainv = tbl.q_Ainv;
energy_meV = tbl.energy_meV;
if ismember('R2', tbl.Properties.VariableNames)
    confidence = tbl.R2;
else
    confidence = ones(size(q_Ainv));
end
end


function fit = local_fit_quasi2d_epsilon_bg(q_Ainv, energy_meV, confidence, epsilon_bg)
q_raw = double(q_Ainv(:));
q = abs(q_raw);
E = double(energy_meV(:));
w = double(confidence(:));
if numel(w) ~= numel(q)
    w = ones(size(q));
end

valid = isfinite(q) & isfinite(E) & q > 0 & E > 0 & isfinite(w);
q_raw = q_raw(valid);
q = q(valid);
E = E(valid);
w = w(valid);

positive_w = w(isfinite(w) & w > 0);
if isempty(positive_w)
    w(:) = 1;
else
    floor_w = max(median(positive_w) * 0.05, eps);
    w(~isfinite(w) | w <= 0) = floor_w;
    w = w ./ median(w);
end

if numel(q) < 3
    error('run_epsilon_bg_sensitivity:InsufficientData', ...
        'Need at least 3 valid branch points for epsilon_bg fitting.');
end

rho0_max_A = 5000;
rho0_min_A = 0.1;
model_fn = @(p, q_in) sqrt(abs(p(1)) .* abs(q_in) ./ ...
    (epsilon_bg + abs(p(2)) .* abs(q_in)));

p0 = local_initial_guess(q, E, w, epsilon_bg, rho0_min_A, rho0_max_A);
lb = [0, rho0_min_A];
ub = [Inf, rho0_max_A];

try
    fit_opts = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxFunctionEvaluations', 10000, ...
        'MaxIterations', 2000, ...
        'FunctionTolerance', 1e-12, ...
        'StepTolerance', 1e-12);
    weighted_model = @(p, q_in) sqrt(w) .* model_fn(p, q_in);
    weighted_data = sqrt(w) .* E;
    p_fit = lsqcurvefit(weighted_model, p0, q, weighted_data, ...
        lb, ub, fit_opts);
catch
    cost = @(x) local_log_param_cost(x, model_fn, q, E, w, ...
        rho0_min_A, rho0_max_A);
    fmin_opts = optimset('Display', 'off', 'MaxFunEvals', 10000, ...
        'MaxIter', 2000, 'TolFun', 1e-12, 'TolX', 1e-12);
    x0 = log(max(p0, [eps, rho0_min_A]));
    x_fit = fminsearch(cost, x0, fmin_opts);
    p_fit = exp(x_fit);
    p_fit(2) = min(max(p_fit(2), rho0_min_A), rho0_max_A);
end

p_fit = abs(p_fit);
A_fit = p_fit(1);
rho0_fit = p_fit(2);
E_pred = model_fn(p_fit, q);
residuals = E - E_pred;
SS_res = sum(w .* residuals .^ 2);
SS_tot = sum(w .* (E - mean(E)) .^ 2);
R_squared = 1 - SS_res / max(SS_tot, eps);
RMSE_meV = sqrt(mean(residuals .^ 2));
E_flat_meV = sqrt(A_fit / rho0_fit);
q_c_Ainv = epsilon_bg / rho0_fit;

fit = struct();
fit.A_fit = A_fit;
fit.rho0_A = rho0_fit;
fit.q_c_Ainv = q_c_Ainv;
fit.E_flat_meV = E_flat_meV;
fit.R_squared = R_squared;
fit.RMSE_meV = RMSE_meV;
fit.n_points = numel(q);
fit.q_abs_min_Ainv = min(q);
fit.q_abs_max_Ainv = max(q);
fit.rho0_at_upper_bound = abs(rho0_fit - rho0_max_A) < 1e-6;
fit.q_data = q_raw;
fit.E_data = E;
fit.E_pred = E_pred;
end


function p0 = local_initial_guess(q, E, w, epsilon_bg, rho0_min_A, rho0_max_A)
rho_candidates = [0.5, 1, 2, 5, 10, 25, 50, 100, ...
    250, 500, 1000, 2000, 3500];
rho_candidates = rho_candidates .* max(epsilon_bg, 1);
rho_candidates = unique(min(max(rho_candidates, rho0_min_A), rho0_max_A));

best_cost = Inf;
p0 = [max(E) ^ 2 * epsilon_bg / max(q), max(10 * epsilon_bg, rho0_min_A)];

for i = 1:numel(rho_candidates)
    rho = rho_candidates(i);
    x = q ./ (epsilon_bg + rho .* q);
    A = max(sum(w .* x .* E .^ 2) / max(sum(w .* x .^ 2), eps), eps);
    pred = sqrt(A .* q ./ (epsilon_bg + rho .* q));
    cost = sum(w .* (pred - E) .^ 2);
    if cost < best_cost
        best_cost = cost;
        p0 = [A, rho];
    end
end
end


function cost = local_log_param_cost(x, model_fn, q, E, w, rho0_min_A, rho0_max_A)
p = exp(x(:))';
p(2) = min(max(p(2), rho0_min_A), rho0_max_A);
residuals = model_fn(p, q) - E;
cost = sum(w .* residuals .^ 2);
if p(2) <= rho0_min_A || p(2) >= rho0_max_A
    cost = cost + 1e6;
end
end


function rows = local_append_row(rows, ds, branch_id, epsilon_bg, fit)
idx = numel(rows) + 1;
rows(idx).session_key = ds.session_key;
rows(idx).session_label = ds.session_label;
rows(idx).branch = branch_id;
rows(idx).epsilon_bg = epsilon_bg;
rows(idx).epsilon_s_equiv = 2 * epsilon_bg - 1;
rows(idx).n_points = fit.n_points;
rows(idx).q_abs_min_Ainv = fit.q_abs_min_Ainv;
rows(idx).q_abs_max_Ainv = fit.q_abs_max_Ainv;
rows(idx).A_fit = fit.A_fit;
rows(idx).rho0_A = fit.rho0_A;
rows(idx).q_c_Ainv = fit.q_c_Ainv;
rows(idx).E_flat_meV = fit.E_flat_meV;
rows(idx).R_squared = fit.R_squared;
rows(idx).RMSE_meV = fit.RMSE_meV;
rows(idx).rho0_at_upper_bound = fit.rho0_at_upper_bound;
end


function local_plot_metric(summary, output_dir, metric_name, y_label, file_name)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 900, 520]);
hold on;
grid on;
box on;

keys = cell(height(summary), 1);
for i = 1:height(summary)
    keys{i} = sprintf('%s B%d', summary.session_key{i}, summary.branch(i));
end
unique_keys = unique(keys, 'stable');
colors = lines(numel(unique_keys));

for ki = 1:numel(unique_keys)
    mask = strcmp(keys, unique_keys{ki});
    sub = summary(mask, :);
    [x, order] = sort(sub.epsilon_bg);
    y = sub.(metric_name);
    y = y(order);
    plot(x, y, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
        'DisplayName', strrep(unique_keys{ki}, '_', '\_'), ...
        'Color', colors(ki, :));
end

xlabel('epsilon\_bg');
ylabel(y_label);
title(strrep(file_name(1:end-4), '_', ' '));
legend('Location', 'bestoutside');

out_path = fullfile(output_dir, file_name);
if exist('exportgraphics', 'file')
    exportgraphics(fig, out_path, 'Resolution', 200);
else
    saveas(fig, out_path);
end
close(fig);
end


function local_write_report(summary, report_path)
fid = fopen(report_path, 'w');
if fid < 0
    error('run_epsilon_bg_sensitivity:ReportOpenFailed', ...
        'Unable to open report for writing: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Epsilon Background Sensitivity Report\n\n');
fprintf(fid, 'This report refits B1 and B3 branch points from the current ');
fprintf(fid, 'Area-normalized Fano-apex pipeline with a direct quasi-2D model:\n\n');
fprintf(fid, '`E(q) = sqrt(A |q| / (epsilon_bg + rho0 |q|))`.\n\n');
fprintf(fid, 'B2 is not refit here because the previous model comparison selected ');
fprintf(fid, '`optical_constant` for B2 in all three data sets.\n\n');
fprintf(fid, '## Sweep setup\n\n');
fprintf(fid, '- epsilon_bg values: 1, 4.5, 10, 15\n');
fprintf(fid, '- equivalent code epsilon_s: `epsilon_s = 2 * epsilon_bg - 1`\n');
fprintf(fid, '- rho0 upper bound: 5000 A\n');
fprintf(fid, '- branch weights: exported per-point R2, normalized by median R2\n');
fprintf(fid, '- source data: branch1_points.csv and branch3_points.csv\n\n');

fprintf(fid, '## Key numeric summary\n\n');
fprintf(fid, '| Session | Branch | rho0 at eps=1 (A) | rho0 at eps=10 (A) | ');
fprintf(fid, 'qc at eps=1 (1/A) | qc at eps=10 (1/A) | Eflat at eps=10 (meV) | R2 at eps=10 |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|\n');

groups = local_summary_groups(summary);
for gi = 1:numel(groups)
    mask = strcmp(summary.session_key, groups(gi).session_key) & ...
        summary.branch == groups(gi).branch;
    sub = summary(mask, :);
    eps1 = sub(abs(sub.epsilon_bg - 1) < 1e-9, :);
    eps10 = sub(abs(sub.epsilon_bg - 10) < 1e-9, :);
    if isempty(eps1) || isempty(eps10)
        continue;
    end
    fprintf(fid, '| %s | B%d | %.3g | %.3g | %.4g | %.4g | %.1f | %.4f |\n', ...
        eps1.session_label{1}, eps1.branch(1), eps1.rho0_A(1), ...
        eps10.rho0_A(1), eps1.q_c_Ainv(1), eps10.q_c_Ainv(1), ...
        eps10.E_flat_meV(1), eps10.R_squared(1));
end

fprintf(fid, '\n## Interpretation\n\n');
fprintf(fid, '1. Increasing epsilon_bg mainly rescales rho0 upward. ');
fprintf(fid, 'For a fixed branch, qc and Eflat remain nearly unchanged when the ');
fprintf(fid, 'fit is not bound-limited, because qc = epsilon_bg / rho0.\n');
fprintf(fid, '2. The default epsilon_bg = 1 result should be treated as a ');
fprintf(fid, 'suspended/vacuum-like baseline. For MoS2 encapsulation, a larger ');
fprintf(fid, 'effective epsilon_bg gives a larger absolute rho0.\n');
fprintf(fid, '3. Relative comparisons between the three data sets are more robust ');
fprintf(fid, 'than the absolute rho0 value, provided the same encapsulation ');
fprintf(fid, 'environment is assumed for all data sets.\n');
fprintf(fid, '4. This sensitivity sweep does not rescue low-confidence 20w B1 ');
fprintf(fid, 'high-q points. It only changes the dielectric-background ');
fprintf(fid, 'interpretation of fitted parameters.\n\n');

fprintf(fid, '## Full table\n\n');
fprintf(fid, '| Session | Branch | epsilon_bg | epsilon_s_equiv | N | rho0 (A) | qc (1/A) | Eflat (meV) | R2 | RMSE (meV) |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:height(summary)
    fprintf(fid, '| %s | B%d | %.3g | %.3g | %d | %.4g | %.4g | %.1f | %.4f | %.1f |\n', ...
        summary.session_label{i}, summary.branch(i), summary.epsilon_bg(i), ...
        summary.epsilon_s_equiv(i), summary.n_points(i), summary.rho0_A(i), ...
        summary.q_c_Ainv(i), summary.E_flat_meV(i), summary.R_squared(i), ...
        summary.RMSE_meV(i));
end
end


function groups = local_summary_groups(summary)
groups = struct('session_key', {}, 'branch', {});
for i = 1:height(summary)
    exists = false;
    for gi = 1:numel(groups)
        if strcmp(groups(gi).session_key, summary.session_key{i}) && ...
                groups(gi).branch == summary.branch(i)
            exists = true;
            break;
        end
    end
    if ~exists
        groups(end + 1).session_key = summary.session_key{i}; %#ok<AGROW>
        groups(end).branch = summary.branch(i);
    end
end
end
