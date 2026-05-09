function output = run_b1_physical_fit_analysis(options)
%RUN_B1_PHYSICAL_FIT_ANALYSIS Thickness-constrained B1 quasi-2D fit.

arguments
    options.datasets = []
    options.branchFileName {mustBeTextScalar} = "branch1_points.csv"
    options.outputTag {mustBeTextScalar} = "b1_physical_fit_260507"
    options.filePrefix {mustBeTextScalar} = "b1"
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

results_root = fullfile(project_root, 'paper_results');
branch_file_name = char(string(options.branchFileName));
file_prefix = char(string(options.filePrefix));
output_dir = fullfile(results_root, char(string(options.outputTag)));
if ~isfolder(output_dir)
    mkdir(output_dir);
end

primary_epsilon_bg = 4.5;
epsilon_bg_values = [1, 4.5, 10, 15];
q_abs_max_Ainv = 0.15;
r2_min_for_fit = 0.7;
sigma_floor_meV = 10;

datasets = local_dataset_config(results_root);
if ~isempty(options.datasets)
    datasets = options.datasets;
end
raw_points = local_load_b1_points(datasets, q_abs_max_Ainv, ...
    branch_file_name);
qabs_points = local_qabs_average(raw_points, datasets, ...
    r2_min_for_fit, sigma_floor_meV);

points_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_physical_fit_points_qabs.csv', file_prefix));
writetable(qabs_points, points_csv);

fit_points = qabs_points(qabs_points.include_for_fit, :);
if height(fit_points) < 8
    error('run_b1_physical_fit_analysis:InsufficientFitPoints', ...
        'Need at least 8 included |q|-averaged B1 points.');
end

main_fit = local_fit_model(fit_points, datasets, ...
    'thickness_free_ratio', primary_epsilon_bg);
comparison = local_model_comparison(fit_points, datasets, primary_epsilon_bg);
sensitivity = local_epsilon_sensitivity(fit_points, datasets, ...
    epsilon_bg_values);
summary = local_main_summary(main_fit, fit_points, datasets);

summary_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_physical_fit_summary.csv', file_prefix));
comparison_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_model_comparison.csv', file_prefix));
sensitivity_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_epsilon_bg_sensitivity.csv', file_prefix));
writetable(summary, summary_csv);
writetable(comparison, comparison_csv);
writetable(sensitivity, sensitivity_csv);

panels_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_physical_fit_panels.png', file_prefix));
panels_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_physical_fit_panels.pdf', file_prefix));
residuals_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_residuals.png', file_prefix));
residuals_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_residuals.pdf', file_prefix));
ratio_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_rho0_thickness_ratio.png', file_prefix));
ratio_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_rho0_thickness_ratio.pdf', file_prefix));

local_plot_fit_panels(qabs_points, fit_points, datasets, main_fit, ...
    comparison, panels_png, panels_pdf);
local_plot_residuals(qabs_points, datasets, main_fit, ...
    residuals_png, residuals_pdf);
local_plot_rho0_ratio(main_fit, ratio_png, ratio_pdf);

mat_path = fullfile(output_dir, local_prefixed_file( ...
    'b1_physical_fit_results.mat', file_prefix));
save(mat_path, 'raw_points', 'qabs_points', 'fit_points', 'main_fit', ...
    'comparison', 'sensitivity', 'summary', 'datasets', ...
    'primary_epsilon_bg', 'epsilon_bg_values');

output = struct();
output.output_dir = output_dir;
output.points_csv = points_csv;
output.summary_csv = summary_csv;
output.comparison_csv = comparison_csv;
output.sensitivity_csv = sensitivity_csv;
output.panels_png = panels_png;
output.panels_pdf = panels_pdf;
output.residuals_png = residuals_png;
output.residuals_pdf = residuals_pdf;
output.ratio_png = ratio_png;
output.ratio_pdf = ratio_pdf;
output.mat_path = mat_path;
output.summary = summary;
output.comparison = comparison;
output.sensitivity = sensitivity;
output.main_fit = main_fit;

fprintf('B1 physical fit analysis complete.\n');
fprintf('  Output directory: %s\n', output_dir);
fprintf('  Included |q|-averaged points: %d\n', height(fit_points));
fprintf('  Main epsilon_bg: %.3g\n', primary_epsilon_bg);
fprintf('  rho0_1film = %.4g A\n', main_fit.rho0_1film_A);
fprintf('  rho0_2film = %.4g A\n', main_fit.rho0_2film_A);
fprintf('  ratio r = %.4g\n', main_fit.thickness_ratio_r);
fprintf('  R2 = %.5f, RMSE = %.2f meV\n', ...
    main_fit.R_squared, main_fit.RMSE_meV);
fprintf('  Summary: %s\n', summary_csv);
end


function datasets = local_dataset_config(results_root)
datasets = struct( ...
    'session_key', {}, ...
    'session_label', {}, ...
    'input_dir', {}, ...
    'thickness_class', {}, ...
    'thickness_factor', {}, ...
    'color', {}, ...
    'marker', {});

datasets(end + 1) = struct( ...
    'session_key', '590_PL2_10w', ...
    'session_label', '590 10w defocus 1film', ...
    'input_dir', fullfile(results_root, '590_gui_history_area_260506'), ...
    'thickness_class', '1film', ...
    'thickness_factor', 1, ...
    'color', [0.120, 0.470, 0.900], ...
    'marker', 'o');

datasets(end + 1) = struct( ...
    'session_key', 'n0_PL2_10w_repeat', ...
    'session_label', 'n0 10w defocus repeat 1film', ...
    'input_dir', fullfile(results_root, 'n0_PL2_10w_gui_history_area_260506'), ...
    'thickness_class', '1film', ...
    'thickness_factor', 1, ...
    'color', [0.160, 0.500, 0.220], ...
    'marker', 'o');

datasets(end + 1) = struct( ...
    'session_key', 'no_PL2_20w_2film', ...
    'session_label', '20w defocus 2film', ...
    'input_dir', fullfile(results_root, ...
        'no_PL2_20w_2film_gui_history_area_260506_highq_refined'), ...
    'thickness_class', '2film', ...
    'thickness_factor', 2, ...
    'color', [0.930, 0.280, 0.300], ...
    'marker', 'o');
end


function raw = local_load_b1_points(datasets, q_abs_max_Ainv, ...
    branch_file_name)
raw = table();

for i = 1:numel(datasets)
    csv_path = fullfile(datasets(i).input_dir, branch_file_name);
    if ~isfile(csv_path)
        error('run_b1_physical_fit_analysis:MissingInput', ...
            'Missing B1 point CSV: %s', csv_path);
    end

    tbl = readtable(csv_path);
    local_require_columns(tbl, {'q_Ainv', 'energy_meV'});

    n = height(tbl);
    q_signed = double(tbl.q_Ainv);
    q_abs = abs(q_signed);
    energy_meV = double(tbl.energy_meV);
    R2 = local_numeric_column(tbl, 'R2', NaN(n, 1));
    E_ci_half_meV = local_energy_ci_half(tbl);

    valid = isfinite(q_abs) & q_abs > 0 & q_abs <= q_abs_max_Ainv & ...
        isfinite(energy_meV) & energy_meV > 0;

    part = table( ...
        repmat(i, n, 1), ...
        repmat({datasets(i).session_key}, n, 1), ...
        repmat({datasets(i).session_label}, n, 1), ...
        repmat({datasets(i).thickness_class}, n, 1), ...
        repmat(datasets(i).thickness_factor, n, 1), ...
        repmat({csv_path}, n, 1), ...
        q_signed, ...
        q_abs, ...
        energy_meV, ...
        energy_meV ./ 1000, ...
        R2, ...
        E_ci_half_meV, ...
        valid, ...
        'VariableNames', {'session_index', 'session_key', ...
        'session_label', 'thickness_class', 'thickness_factor', ...
        'source_csv', 'q_signed_Ainv', 'q_abs_Ainv', ...
        'energy_meV', 'energy_eV', 'R2', 'E_ci_half_meV', ...
        'in_q_window'});

    raw = [raw; part]; %#ok<AGROW>
end
end


function local_require_columns(tbl, required)
names = tbl.Properties.VariableNames;
for i = 1:numel(required)
    if ~ismember(required{i}, names)
        error('run_b1_physical_fit_analysis:MissingColumn', ...
            'B1 CSV is missing required column "%s".', required{i});
    end
end
end


function col = local_numeric_column(tbl, name, fallback)
if ismember(name, tbl.Properties.VariableNames)
    col = double(tbl.(name));
else
    col = fallback;
end
end


function ci = local_energy_ci_half(tbl)
n = height(tbl);
if ismember('E_ci_half_meV', tbl.Properties.VariableNames)
    ci = double(tbl.E_ci_half_meV);
elseif all(ismember({'E_ci_lo', 'E_ci_hi'}, tbl.Properties.VariableNames))
    ci = abs(double(tbl.E_ci_hi) - double(tbl.E_ci_lo)) ./ 2;
else
    ci = NaN(n, 1);
end
end


function averaged = local_qabs_average(raw, datasets, r2_min_for_fit, ...
    sigma_floor_meV)
averaged = table();

for i = 1:numel(datasets)
    sub = raw(raw.session_index == i & raw.in_q_window, :);
    if isempty(sub)
        continue;
    end

    q_group = round(sub.q_abs_Ainv, 6);
    q_values = unique(q_group);
    q_values = q_values(q_values > 0);

    for qi = 1:numel(q_values)
        q_mask = abs(q_group - q_values(qi)) < 1e-12;
        group_all = sub(q_mask, :);
        fit_raw_mask = isfinite(group_all.R2) & ...
            group_all.R2 >= r2_min_for_fit;
        group_fit = group_all(fit_raw_mask, :);

        include_for_fit = height(group_fit) > 0;
        if include_for_fit
            source = group_fit;
        else
            source = group_all;
        end

        E = double(source.energy_meV);
        finite_E = E(isfinite(E));
        if isempty(finite_E)
            continue;
        end

        ci_half = double(source.E_ci_half_meV);
        ci_half = ci_half(isfinite(ci_half) & ci_half >= 0);
        if isempty(ci_half)
            ci_half_rms = NaN;
            sigma_ci = NaN;
        else
            ci_half_rms = sqrt(mean(ci_half .^ 2));
            sigma_ci = sqrt(mean((ci_half ./ 1.96) .^ 2));
        end

        if numel(finite_E) > 1
            pair_spread = std(finite_E, 0);
        else
            pair_spread = 0;
        end

        if ~isfinite(sigma_ci)
            sigma_ci = max(pair_spread, sigma_floor_meV);
        end
        sigma_fit = sqrt(sigma_ci .^ 2 + pair_spread .^ 2);
        sigma_fit = max(sigma_fit, sigma_floor_meV);
        energy_err = max([ci_half_rms, pair_spread, sigma_floor_meV], ...
            [], 'omitnan');

        R2_vals = double(group_all.R2);
        R2_fit_vals = double(group_fit.R2);
        row = table( ...
            i, ...
            {datasets(i).session_key}, ...
            {datasets(i).session_label}, ...
            {datasets(i).thickness_class}, ...
            datasets(i).thickness_factor, ...
            q_values(qi), ...
            mean(finite_E, 'omitnan'), ...
            mean(finite_E, 'omitnan') ./ 1000, ...
            energy_err, ...
            energy_err ./ 1000, ...
            sigma_fit, ...
            height(group_all), ...
            height(group_fit), ...
            min(group_all.q_signed_Ainv), ...
            max(group_all.q_signed_Ainv), ...
            mean(R2_vals, 'omitnan'), ...
            min(R2_vals, [], 'omitnan'), ...
            mean(R2_fit_vals, 'omitnan'), ...
            ci_half_rms, ...
            pair_spread, ...
            include_for_fit, ...
            'VariableNames', {'session_index', 'session_key', ...
            'session_label', 'thickness_class', 'thickness_factor', ...
            'q_abs_Ainv', 'energy_mean_meV', 'energy_mean_eV', ...
            'energy_err_meV', 'energy_err_eV', 'sigma_fit_meV', ...
            'n_raw_points', 'n_fit_raw_points', 'q_signed_min_Ainv', ...
            'q_signed_max_Ainv', 'R2_mean_all', 'R2_min_all', ...
            'R2_mean_fit_raw', 'E_ci_half_rms_meV', ...
            'pair_spread_meV', 'include_for_fit'});

        averaged = [averaged; row]; %#ok<AGROW>
    end
end

averaged = sortrows(averaged, {'session_index', 'q_abs_Ainv'});
end


function comparison = local_model_comparison(points, datasets, epsilon_bg)
model_keys = {'thickness_free_ratio', 'independent_rho0', ...
    'shared_rho0', 'fixed_double_ratio'};
comparison = table();

for i = 1:numel(model_keys)
    fit = local_fit_model(points, datasets, model_keys{i}, epsilon_bg);
    comparison = [comparison; local_comparison_row(fit, datasets)]; %#ok<AGROW>
end
end


function sensitivity = local_epsilon_sensitivity(points, datasets, ...
    epsilon_bg_values)
sensitivity = table();

for i = 1:numel(epsilon_bg_values)
    fit = local_fit_model(points, datasets, 'thickness_free_ratio', ...
        epsilon_bg_values(i));
    row = local_comparison_row(fit, datasets);
    sensitivity = [sensitivity; row]; %#ok<AGROW>
end
end


function fit = local_fit_model(points, datasets, model_key, epsilon_bg)
init = local_initial_from_single_sessions(points, datasets, epsilon_bg);
x0 = local_initial_x(init, datasets, model_key);

cost_fn = @(x) local_fit_cost(x, points, datasets, model_key, epsilon_bg);
opts = optimset('Display', 'off', 'MaxFunEvals', 50000, ...
    'MaxIter', 10000, 'TolFun', 1e-12, 'TolX', 1e-12);
x_fit = fminsearch(cost_fn, x0, opts);

params = local_params_from_x(x_fit, datasets, model_key);
[E_pred, residuals] = local_predict(points, params, epsilon_bg);
sigma = max(double(points.sigma_fit_meV), eps);
w = 1 ./ (sigma .^ 2);
chi2 = sum((residuals ./ sigma) .^ 2);
n_points = height(points);
n_params = numel(x_fit);
dof = max(n_points - n_params, 1);
weighted_mean = sum(w .* points.energy_mean_meV) / sum(w);
SS_tot = sum(w .* (points.energy_mean_meV - weighted_mean) .^ 2);
SS_res = sum(w .* residuals .^ 2);
R_squared = 1 - SS_res / max(SS_tot, eps);
RMSE_meV = sqrt(mean(residuals .^ 2));
AIC = n_points * log(max(chi2 / n_points, eps)) + 2 * n_params;
BIC = n_points * log(max(chi2 / n_points, eps)) + ...
    n_params * log(n_points);

fit = struct();
fit.model_key = model_key;
fit.model_label = local_model_label(model_key);
fit.epsilon_bg = epsilon_bg;
fit.x_fit = x_fit;
fit.n_points = n_points;
fit.n_params = n_params;
fit.dof = dof;
fit.chi2 = chi2;
fit.reduced_chi2 = chi2 / dof;
fit.R_squared = R_squared;
fit.RMSE_meV = RMSE_meV;
fit.AIC = AIC;
fit.BIC = BIC;
fit.A_by_session = params.A_by_session;
fit.rho_by_session_A = params.rho_by_session_A;
fit.rho0_1film_A = params.rho0_1film_A;
fit.rho0_2film_A = params.rho0_2film_A;
fit.thickness_ratio_r = params.thickness_ratio_r;
fit.E_pred_meV = E_pred;
fit.residuals_meV = residuals;
end


function init = local_initial_from_single_sessions(points, datasets, ...
    epsilon_bg)
n_sessions = numel(datasets);
A = zeros(n_sessions, 1);
rho = zeros(n_sessions, 1);

for i = 1:n_sessions
    sub = points(points.session_index == i, :);
    [A(i), rho(i)] = local_single_session_fit(sub, epsilon_bg);
end

rho_1film = mean(rho([datasets.thickness_factor] == 1));
rho_2film = mean(rho([datasets.thickness_factor] == 2));
if ~isfinite(rho_1film) || rho_1film <= 0
    rho_1film = max(mean(rho), 1);
end
if ~isfinite(rho_2film) || rho_2film <= 0
    rho_2film = max(2 * rho_1film, 1);
end

init = struct();
init.A_by_session = max(A, eps);
init.rho_by_session_A = max(rho, 0.01);
init.rho0_1film_A = max(rho_1film, 0.01);
init.rho0_2film_A = max(rho_2film, init.rho0_1film_A * 1.25);
init.thickness_ratio_r = max(init.rho0_2film_A / ...
    init.rho0_1film_A, 1.25);
end


function [A_best, rho_best] = local_single_session_fit(points, epsilon_bg)
q = double(points.q_abs_Ainv);
E = double(points.energy_mean_meV);
sigma = max(double(points.sigma_fit_meV), eps);
w = 1 ./ (sigma .^ 2);

rho_grid = logspace(log10(0.05), log10(5000), 400);
best_cost = Inf;
A_best = max(E) ^ 2 * 10;
rho_best = 10;

for i = 1:numel(rho_grid)
    rho = rho_grid(i);
    c = sqrt(q ./ (epsilon_bg + rho .* q));
    alpha = sum(w .* c .* E) / max(sum(w .* c .^ 2), eps);
    A = max(alpha .^ 2, eps);
    residuals = E - alpha .* c;
    cost = sum(w .* residuals .^ 2);
    if cost < best_cost
        best_cost = cost;
        A_best = A;
        rho_best = rho;
    end
end
end


function x0 = local_initial_x(init, datasets, model_key)
switch model_key
    case 'thickness_free_ratio'
        x0 = [ ...
            log(init.rho0_1film_A); ...
            log(init.thickness_ratio_r - 1); ...
            log(init.A_by_session(:))];
    case 'independent_rho0'
        x0 = [log(init.rho_by_session_A(:)); ...
            log(init.A_by_session(:))];
    case 'shared_rho0'
        shared_rho = mean(init.rho_by_session_A);
        x0 = [log(shared_rho); log(init.A_by_session(:))];
    case 'fixed_double_ratio'
        x0 = [log(init.rho0_1film_A); log(init.A_by_session(:))];
    otherwise
        error('run_b1_physical_fit_analysis:UnknownModel', ...
            'Unknown model key: %s', model_key);
end

if any(~isfinite(x0))
    error('run_b1_physical_fit_analysis:InvalidInitialGuess', ...
        'Invalid initial guess for %s.', model_key);
end

if numel(datasets) ~= 3
    error('run_b1_physical_fit_analysis:DatasetCountChanged', ...
        'This analysis expects the current three B1 datasets.');
end
end


function cost = local_fit_cost(x, points, datasets, model_key, epsilon_bg)
if any(~isfinite(x)) || any(abs(x) > 80)
    cost = realmax('double') / 1e6;
    return;
end

params = local_params_from_x(x, datasets, model_key);
[~, residuals] = local_predict(points, params, epsilon_bg);
sigma = max(double(points.sigma_fit_meV), eps);
cost = sum((residuals ./ sigma) .^ 2);

if ~isfinite(cost)
    cost = realmax('double') / 1e6;
end
end


function params = local_params_from_x(x, datasets, model_key)
n_sessions = numel(datasets);
thickness_factor = [datasets.thickness_factor]';

switch model_key
    case 'thickness_free_ratio'
        rho0_1film = exp(x(1));
        r = 1 + exp(x(2));
        rho0_2film = r * rho0_1film;
        rho_by_session = zeros(n_sessions, 1);
        rho_by_session(thickness_factor == 1) = rho0_1film;
        rho_by_session(thickness_factor == 2) = rho0_2film;
        A_by_session = exp(x(3:(2 + n_sessions)));

    case 'independent_rho0'
        rho_by_session = exp(x(1:n_sessions));
        A_by_session = exp(x((n_sessions + 1):(2 * n_sessions)));
        rho0_1film = mean(rho_by_session(thickness_factor == 1));
        rho0_2film = mean(rho_by_session(thickness_factor == 2));
        r = rho0_2film / rho0_1film;

    case 'shared_rho0'
        rho0_shared = exp(x(1));
        rho_by_session = repmat(rho0_shared, n_sessions, 1);
        A_by_session = exp(x(2:(1 + n_sessions)));
        rho0_1film = rho0_shared;
        rho0_2film = rho0_shared;
        r = 1;

    case 'fixed_double_ratio'
        rho0_1film = exp(x(1));
        r = 2;
        rho0_2film = r * rho0_1film;
        rho_by_session = zeros(n_sessions, 1);
        rho_by_session(thickness_factor == 1) = rho0_1film;
        rho_by_session(thickness_factor == 2) = rho0_2film;
        A_by_session = exp(x(2:(1 + n_sessions)));

    otherwise
        error('run_b1_physical_fit_analysis:UnknownModel', ...
            'Unknown model key: %s', model_key);
end

params = struct();
params.A_by_session = A_by_session(:);
params.rho_by_session_A = rho_by_session(:);
params.rho0_1film_A = rho0_1film;
params.rho0_2film_A = rho0_2film;
params.thickness_ratio_r = r;
end


function [E_pred, residuals] = local_predict(points, params, epsilon_bg)
q = double(points.q_abs_Ainv);
session_index = double(points.session_index);
A = params.A_by_session(session_index);
rho = params.rho_by_session_A(session_index);
E_pred = sqrt(A .* q ./ (epsilon_bg + rho .* q));
residuals = double(points.energy_mean_meV) - E_pred;
end


function label = local_model_label(model_key)
switch model_key
    case 'thickness_free_ratio'
        label = '1film shared rho0, 2film free ratio';
    case 'independent_rho0'
        label = 'Independent rho0 per dataset';
    case 'shared_rho0'
        label = 'Single shared rho0';
    case 'fixed_double_ratio'
        label = 'rho0_2film fixed at 2x';
    otherwise
        label = model_key;
end
end


function row = local_comparison_row(fit, datasets)
session_keys = {datasets.session_key};
rho = fit.rho_by_session_A(:);
A = fit.A_by_session(:);
E_flat = sqrt(A ./ rho);
q_c = fit.epsilon_bg ./ rho;

row = table( ...
    {fit.model_key}, ...
    {fit.model_label}, ...
    fit.epsilon_bg, ...
    fit.n_points, ...
    fit.n_params, ...
    fit.dof, ...
    fit.chi2, ...
    fit.reduced_chi2, ...
    fit.R_squared, ...
    fit.RMSE_meV, ...
    fit.AIC, ...
    fit.BIC, ...
    fit.rho0_1film_A, ...
    fit.rho0_2film_A, ...
    fit.thickness_ratio_r, ...
    rho(1), rho(2), rho(3), ...
    A(1), A(2), A(3), ...
    E_flat(1), E_flat(2), E_flat(3), ...
    q_c(1), q_c(2), q_c(3), ...
    'VariableNames', {'model_key', 'model_label', 'epsilon_bg', ...
    'n_points', 'n_params', 'dof', 'chi2', 'reduced_chi2', ...
    'R_squared', 'RMSE_meV', 'AIC', 'BIC', 'rho0_1film_A', ...
    'rho0_2film_A', 'thickness_ratio_r', ...
    ['rho0_' session_keys{1} '_A'], ...
    ['rho0_' session_keys{2} '_A'], ...
    ['rho0_' session_keys{3} '_A'], ...
    ['A_' session_keys{1}], ...
    ['A_' session_keys{2}], ...
    ['A_' session_keys{3}], ...
    ['E_flat_' session_keys{1} '_meV'], ...
    ['E_flat_' session_keys{2} '_meV'], ...
    ['E_flat_' session_keys{3} '_meV'], ...
    ['q_c_' session_keys{1} '_Ainv'], ...
    ['q_c_' session_keys{2} '_Ainv'], ...
    ['q_c_' session_keys{3} '_Ainv']});
end


function summary = local_main_summary(fit, points, datasets)
summary = table();
params = struct();
params.A_by_session = fit.A_by_session;
params.rho_by_session_A = fit.rho_by_session_A;
[E_pred, residuals] = local_predict(points, params, fit.epsilon_bg);

for i = 1:numel(datasets)
    mask = points.session_index == i;
    sub = points(mask, :);
    r = residuals(mask);
    pred = E_pred(mask);
    E = double(sub.energy_mean_meV);
    w = 1 ./ max(double(sub.sigma_fit_meV), eps) .^ 2;
    wmean = sum(w .* E) / sum(w);
    R2 = 1 - sum(w .* r .^ 2) / max(sum(w .* (E - wmean) .^ 2), eps);
    RMSE = sqrt(mean(r .^ 2));
    rho0 = fit.rho_by_session_A(i);
    A_fit = fit.A_by_session(i);

    row = table( ...
        {fit.model_key}, ...
        fit.epsilon_bg, ...
        i, ...
        {datasets(i).session_key}, ...
        {datasets(i).session_label}, ...
        {datasets(i).thickness_class}, ...
        datasets(i).thickness_factor, ...
        height(sub), ...
        min(sub.q_abs_Ainv), ...
        max(sub.q_abs_Ainv), ...
        A_fit, ...
        rho0, ...
        sqrt(A_fit ./ rho0), ...
        fit.epsilon_bg ./ rho0, ...
        fit.rho0_1film_A, ...
        fit.rho0_2film_A, ...
        fit.thickness_ratio_r, ...
        R2, ...
        RMSE, ...
        mean(r, 'omitnan'), ...
        median(r, 'omitnan'), ...
        fit.R_squared, ...
        fit.RMSE_meV, ...
        fit.reduced_chi2, ...
        mean(pred, 'omitnan'), ...
        'VariableNames', {'model_key', 'epsilon_bg', ...
        'session_index', 'session_key', 'session_label', ...
        'thickness_class', 'thickness_factor', 'n_fit_points', ...
        'q_abs_min_Ainv', 'q_abs_max_Ainv', 'A_fit', ...
        'rho0_A', 'E_flat_meV', 'q_c_Ainv', ...
        'rho0_1film_A', 'rho0_2film_A', 'thickness_ratio_r', ...
        'dataset_R_squared', 'dataset_RMSE_meV', ...
        'residual_mean_meV', 'residual_median_meV', ...
        'global_R_squared', 'global_RMSE_meV', ...
        'global_reduced_chi2', 'E_pred_mean_meV'});
    summary = [summary; row]; %#ok<AGROW>
end
end


function local_plot_fit_panels(all_points, fit_points, datasets, fit, ...
    comparison, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1180, 640]);
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');

ax = nexttile(t, 1);
hold(ax, 'on');
for i = 1:numel(datasets)
    sub_all = all_points(all_points.session_index == i, :);
    sub_fit = sub_all(sub_all.include_for_fit, :);
    sub_drop = sub_all(~sub_all.include_for_fit, :);
    col = datasets(i).color;

    if ~isempty(sub_fit)
        errorbar(ax, sub_fit.q_abs_Ainv, sub_fit.energy_mean_eV, ...
            sub_fit.energy_err_eV, 'LineStyle', 'none', ...
            'Marker', 'none', 'Color', local_lighten(col, 0.65), ...
            'LineWidth', 0.65, 'CapSize', 0, 'HandleVisibility', 'off');
        scatter(ax, sub_fit.q_abs_Ainv, sub_fit.energy_mean_eV, 38, ...
            'Marker', datasets(i).marker, 'MarkerFaceColor', col, ...
            'MarkerEdgeColor', col, 'DisplayName', datasets(i).session_label);
    end

    if ~isempty(sub_drop)
        scatter(ax, sub_drop.q_abs_Ainv, sub_drop.energy_mean_eV, 26, ...
            'Marker', datasets(i).marker, 'MarkerFaceColor', 'none', ...
            'MarkerEdgeColor', local_lighten(col, 0.35), ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end

    q_fit = linspace(max(0.001, min(fit_points.q_abs_Ainv)), 0.15, 250)';
    A = fit.A_by_session(i);
    rho = fit.rho_by_session_A(i);
    E_fit = sqrt(A .* q_fit ./ (fit.epsilon_bg + rho .* q_fit)) ./ 1000;
    plot(ax, q_fit, E_fit, '-', 'Color', col, 'LineWidth', 1.6, ...
        'HandleVisibility', 'off');
end
hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
ax.FontName = 'Arial';
ax.FontSize = 9.5;
xlim(ax, [0, 0.155]);
ylim(ax, [0.45, 1.45]);
xlabel(ax, '|q| (A^{-1})');
ylabel(ax, 'B1 energy (eV)');
title(ax, sprintf('B1 fit, epsilon_{bg}=%.1f', ...
    fit.epsilon_bg));
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off');

ax2 = nexttile(t, 2);
bar(ax2, comparison.RMSE_meV, 'FaceColor', [0.35, 0.38, 0.42]);
box(ax2, 'on');
grid(ax2, 'on');
ax2.GridAlpha = 0.16;
ax2.FontName = 'Arial';
ax2.FontSize = 9.5;
ax2.XTick = 1:height(comparison);
ax2.XTickLabel = {'free r', 'indep.', 'shared', 'r=2'};
ax2.XTickLabelRotation = 0;
ylabel(ax2, 'RMSE (meV)');
title(ax2, 'Model comparison');

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_residuals(all_points, datasets, fit, png_path, pdf_path)
params = struct();
params.A_by_session = fit.A_by_session;
params.rho_by_session_A = fit.rho_by_session_A;
[E_pred, residuals] = local_predict(all_points, params, fit.epsilon_bg);

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 900, 560]);
ax = axes(fig);
hold(ax, 'on');
plot(ax, [0, 0.155], [0, 0], '-', 'Color', [0.2, 0.2, 0.2], ...
    'LineWidth', 0.9, 'HandleVisibility', 'off');

for i = 1:numel(datasets)
    mask = all_points.session_index == i;
    sub = all_points(mask, :);
    r = residuals(mask);
    col = datasets(i).color;
    fit_mask = sub.include_for_fit;

    scatter(ax, sub.q_abs_Ainv(fit_mask), r(fit_mask), 36, ...
        'Marker', datasets(i).marker, 'MarkerFaceColor', col, ...
        'MarkerEdgeColor', col, 'DisplayName', datasets(i).session_label);
    scatter(ax, sub.q_abs_Ainv(~fit_mask), r(~fit_mask), 25, ...
        'Marker', datasets(i).marker, 'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', local_lighten(col, 0.35), ...
        'HandleVisibility', 'off');
end

hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
ax.FontName = 'Arial';
ax.FontSize = 10.5;
xlim(ax, [0, 0.155]);
xlabel(ax, '|q| (A^{-1})');
ylabel(ax, 'Residual E_{obs}-E_{fit} (meV)');
title(ax, 'B1 residuals, thickness-constrained model');
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off');

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_rho0_ratio(fit, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 880, 470]);
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(t, 1);
bar(ax1, [fit.rho0_1film_A, fit.rho0_2film_A], ...
    'FaceColor', [0.33, 0.50, 0.65]);
box(ax1, 'on');
grid(ax1, 'on');
ax1.GridAlpha = 0.16;
ax1.XTick = 1:2;
ax1.XTickLabel = {'1film', '2film'};
ylabel(ax1, 'rho0 (A)');
title(ax1, 'Thickness-dependent rho0');

ax2 = nexttile(t, 2);
bar(ax2, fit.thickness_ratio_r, 'FaceColor', [0.50, 0.42, 0.30]);
hold(ax2, 'on');
plot(ax2, [0.5, 1.5], [2, 2], '--', 'Color', [0.2, 0.2, 0.2], ...
    'LineWidth', 1.0);
hold(ax2, 'off');
box(ax2, 'on');
grid(ax2, 'on');
ax2.GridAlpha = 0.16;
ax2.XTick = 1;
ax2.XTickLabel = {'fit'};
xlim(ax2, [0.5, 1.5]);
ylabel(ax2, 'rho0_{2film}/rho0_{1film}');
title(ax2, 'Thickness ratio');

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function c = local_lighten(color, amount)
c = color + (1 - color) .* amount;
c = min(max(c, 0), 1);
end


function name = local_prefixed_file(default_name, file_prefix)
if strcmp(file_prefix, 'b1')
    name = default_name;
elseif startsWith(default_name, 'b1')
    name = [file_prefix default_name(3:end)];
else
    name = [file_prefix '_' default_name];
end
end
