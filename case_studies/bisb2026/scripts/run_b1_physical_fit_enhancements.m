function output = run_b1_physical_fit_enhancements(options)
%RUN_B1_PHYSICAL_FIT_ENHANCEMENTS Robustness checks for B1 physical fit.

arguments
    options.datasets = []
    options.branchFileName {mustBeTextScalar} = "branch1_points.csv"
    options.outputTag {mustBeTextScalar} = "b1_physical_fit_enhancements_260507"
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
primary_qmax_Ainv = 0.15;
qmax_values_Ainv = [0.10, 0.115, 0.13, 0.15, 0.18, 0.20];
qjitter_sigma_values_Ainv = [0, 0.001, 0.0025, 0.005, 0.01];
n_bootstrap = 300;
n_qjitter = 120;
r2_min_for_fit = 0.7;
sigma_floor_meV = 10;
rng_seed = 260507;

rng(rng_seed, 'twister');

datasets = local_dataset_config(results_root);
if ~isempty(options.datasets)
    datasets = options.datasets;
end
raw_points = local_load_b1_points(datasets, max(qmax_values_Ainv), ...
    branch_file_name);
qabs_points = local_qabs_average(raw_points, datasets, ...
    r2_min_for_fit, sigma_floor_meV);

points_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_enhancement_points_qabs.csv', file_prefix));
writetable(qabs_points, points_csv);

base_points = qabs_points(qabs_points.include_for_fit & ...
    qabs_points.q_abs_Ainv <= primary_qmax_Ainv, :);
base_fit = local_fit_thickness_model(base_points, datasets, ...
    primary_epsilon_bg, 'lsq');

bootstrap_ci = local_bootstrap_ci(base_points, datasets, ...
    primary_epsilon_bg, n_bootstrap);
qrange = local_qrange_stability(qabs_points, datasets, ...
    primary_epsilon_bg, qmax_values_Ainv);
qjitter = local_qjitter_sensitivity(base_points, datasets, ...
    primary_epsilon_bg, qjitter_sigma_values_Ainv, n_qjitter);
huber = local_huber_fit_diagnostic(base_points, datasets, ...
    primary_epsilon_bg);
[residual_structure, residual_points] = local_residual_diagnostics( ...
    qabs_points, datasets, base_fit, primary_epsilon_bg);

bootstrap_ci_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_bootstrap_ci.csv', file_prefix));
qrange_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_qrange_stability.csv', file_prefix));
qjitter_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_defocus_qjitter_sensitivity.csv', file_prefix));
huber_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_huber_robust_fit.csv', file_prefix));
residual_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_residual_structure.csv', file_prefix));
residual_points_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_residual_points.csv', file_prefix));
writetable(bootstrap_ci.summary, bootstrap_ci_csv);
writetable(qrange, qrange_csv);
writetable(qjitter, qjitter_csv);
writetable(huber, huber_csv);
writetable(residual_structure, residual_csv);
writetable(residual_points, residual_points_csv);

bootstrap_samples_csv = fullfile(output_dir, local_prefixed_file( ...
    'b1_bootstrap_samples.csv', file_prefix));
writetable(bootstrap_ci.samples, bootstrap_samples_csv);

bootstrap_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_bootstrap_ratio_hist.png', file_prefix));
bootstrap_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_bootstrap_ratio_hist.pdf', file_prefix));
qmax_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_qmax_stability.png', file_prefix));
qmax_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_qmax_stability.pdf', file_prefix));
qjitter_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_defocus_sensitivity.png', file_prefix));
qjitter_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_defocus_sensitivity.pdf', file_prefix));
huber_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_huber_comparison.png', file_prefix));
huber_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_huber_comparison.pdf', file_prefix));
residual_png = fullfile(output_dir, local_prefixed_file( ...
    'b1_residual_structure.png', file_prefix));
residual_pdf = fullfile(output_dir, local_prefixed_file( ...
    'b1_residual_structure.pdf', file_prefix));

local_plot_bootstrap_ratio(bootstrap_ci.samples, bootstrap_ci.summary, ...
    bootstrap_png, bootstrap_pdf);
local_plot_qrange_stability(qrange, qmax_png, qmax_pdf);
local_plot_qjitter_sensitivity(qjitter, qjitter_png, qjitter_pdf);
local_plot_huber_comparison(huber, huber_png, huber_pdf);
local_plot_residual_structure(residual_points, datasets, ...
    residual_png, residual_pdf);

mat_path = fullfile(output_dir, local_prefixed_file( ...
    'b1_physical_fit_enhancements.mat', file_prefix));
save(mat_path, 'datasets', 'raw_points', 'qabs_points', 'base_points', ...
    'base_fit', 'bootstrap_ci', 'qrange', 'qjitter', 'huber', ...
    'residual_structure', 'residual_points', 'primary_epsilon_bg', ...
    'primary_qmax_Ainv', 'qmax_values_Ainv', ...
    'qjitter_sigma_values_Ainv', 'n_bootstrap', 'n_qjitter', ...
    'rng_seed');

output = struct();
output.output_dir = output_dir;
output.points_csv = points_csv;
output.bootstrap_ci_csv = bootstrap_ci_csv;
output.qrange_csv = qrange_csv;
output.qjitter_csv = qjitter_csv;
output.huber_csv = huber_csv;
output.residual_csv = residual_csv;
output.residual_points_csv = residual_points_csv;
output.bootstrap_samples_csv = bootstrap_samples_csv;
output.bootstrap_png = bootstrap_png;
output.qmax_png = qmax_png;
output.qjitter_png = qjitter_png;
output.huber_png = huber_png;
output.residual_png = residual_png;
output.mat_path = mat_path;
output.base_fit = base_fit;
output.bootstrap_ci = bootstrap_ci.summary;
output.qrange = qrange;
output.qjitter = qjitter;
output.huber = huber;
output.residual_structure = residual_structure;

fprintf('B1 physical fit enhancements complete.\n');
fprintf('  Output directory: %s\n', output_dir);
fprintf('  Base qmax: %.3f A^-1, epsilon_bg: %.3g\n', ...
    primary_qmax_Ainv, primary_epsilon_bg);
fprintf('  Base ratio r = %.4g\n', base_fit.thickness_ratio_r);
fprintf('  Bootstrap samples: %d\n', height(bootstrap_ci.samples));
fprintf('  qmax values: %s\n', mat2str(qmax_values_Ainv));
fprintf('  Summary CI: %s\n', bootstrap_ci_csv);
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
    'input_dir', fullfile(results_root, ...
    '590_gui_history_area_260506_wideq030'), ...
    'thickness_class', '1film', ...
    'thickness_factor', 1, ...
    'color', [0.120, 0.470, 0.900], ...
    'marker', 'o');

datasets(end + 1) = struct( ...
    'session_key', 'n0_PL2_10w_repeat', ...
    'session_label', 'n0 10w defocus repeat 1film', ...
    'input_dir', fullfile(results_root, ...
    'n0_PL2_10w_gui_history_area_260506_wideq030'), ...
    'thickness_class', '1film', ...
    'thickness_factor', 1, ...
    'color', [0.160, 0.500, 0.220], ...
    'marker', 'o');

datasets(end + 1) = struct( ...
    'session_key', 'no_PL2_20w_2film', ...
    'session_label', '20w defocus 2film', ...
    'input_dir', fullfile(results_root, ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined_wideq030'), ...
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
        error('run_b1_physical_fit_enhancements:MissingInput', ...
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
        q_signed, q_abs, energy_meV, energy_meV ./ 1000, R2, ...
        E_ci_half_meV, valid, ...
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
        error('run_b1_physical_fit_enhancements:MissingColumn', ...
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
            i, {datasets(i).session_key}, {datasets(i).session_label}, ...
            {datasets(i).thickness_class}, datasets(i).thickness_factor, ...
            q_values(qi), mean(finite_E, 'omitnan'), ...
            mean(finite_E, 'omitnan') ./ 1000, ...
            energy_err, energy_err ./ 1000, sigma_fit, ...
            height(group_all), height(group_fit), ...
            min(group_all.q_signed_Ainv), max(group_all.q_signed_Ainv), ...
            mean(R2_vals, 'omitnan'), min(R2_vals, [], 'omitnan'), ...
            mean(R2_fit_vals, 'omitnan'), ci_half_rms, pair_spread, ...
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


function ci = local_bootstrap_ci(points, datasets, epsilon_bg, n_bootstrap)
base_fit = local_fit_thickness_model(points, datasets, epsilon_bg, 'lsq');
rows = table();

for bi = 1:n_bootstrap
    sample = table();
    for si = 1:numel(datasets)
        sub = points(points.session_index == si, :);
        idx = randi(height(sub), height(sub), 1);
        sample = [sample; sub(idx, :)]; %#ok<AGROW>
    end

    try
        fit = local_fit_thickness_model(sample, datasets, epsilon_bg, ...
            'lsq', base_fit.x_fit);
        rows = [rows; local_fit_row(fit, sprintf('bootstrap_%03d', bi), ...
            'bootstrap', NaN)]; %#ok<AGROW>
    catch ME
        warning('run_b1_physical_fit_enhancements:BootstrapFailed', ...
            'Bootstrap replicate %d failed: %s', bi, ME.message);
    end
end

summary = local_ci_summary(rows, base_fit);
ci = struct();
ci.samples = rows;
ci.summary = summary;
end


function summary = local_ci_summary(samples, base_fit)
metrics = {'rho0_1film_A', 'rho0_2film_A', 'thickness_ratio_r', ...
    'R_squared', 'RMSE_meV', 'reduced_chi2'};
summary = table();

for i = 1:numel(metrics)
    values = double(samples.(metrics{i}));
    values = values(isfinite(values));
    if isempty(values)
        continue;
    end

    row = table( ...
        {metrics{i}}, ...
        local_metric_value(base_fit, metrics{i}), ...
        mean(values), ...
        median(values), ...
        local_percentile(values, 2.5), ...
        local_percentile(values, 16), ...
        local_percentile(values, 84), ...
        local_percentile(values, 97.5), ...
        numel(values), ...
        'VariableNames', {'metric', 'base_value', 'mean_value', ...
        'median_value', 'ci95_low', 'ci68_low', 'ci68_high', ...
        'ci95_high', 'n_samples'});
    summary = [summary; row]; %#ok<AGROW>
end
end


function value = local_metric_value(fit, metric)
switch metric
    case 'rho0_1film_A'
        value = fit.rho0_1film_A;
    case 'rho0_2film_A'
        value = fit.rho0_2film_A;
    case 'thickness_ratio_r'
        value = fit.thickness_ratio_r;
    case 'R_squared'
        value = fit.R_squared;
    case 'RMSE_meV'
        value = fit.RMSE_meV;
    case 'reduced_chi2'
        value = fit.reduced_chi2;
    otherwise
        value = NaN;
end
end


function qrange = local_qrange_stability(qabs_points, datasets, ...
    epsilon_bg, qmax_values_Ainv)
qrange = table();

for i = 1:numel(qmax_values_Ainv)
    qmax = qmax_values_Ainv(i);
    points = qabs_points(qabs_points.include_for_fit & ...
        qabs_points.q_abs_Ainv <= qmax, :);
    if height(points) < 10
        continue;
    end

    fit = local_fit_thickness_model(points, datasets, epsilon_bg, 'lsq');
    row = local_fit_row(fit, sprintf('qmax_%.3f', qmax), ...
        'qrange', qmax);
    row.n_points = height(points);
    row.q_abs_max_Ainv = qmax;
    row.q_abs_data_max_Ainv = max(points.q_abs_Ainv);
    qrange = [qrange; row]; %#ok<AGROW>
end
end


function qjitter = local_qjitter_sensitivity(points, datasets, ...
    epsilon_bg, sigma_values, n_replicates)
base_fit = local_fit_thickness_model(points, datasets, epsilon_bg, 'lsq');
all_rows = table();

for si = 1:numel(sigma_values)
    sigma_q = sigma_values(si);
    for ri = 1:n_replicates
        jittered = points;
        mask_2film = strcmp(jittered.thickness_class, '2film');
        q_jit = jittered.q_abs_Ainv;
        if sigma_q > 0
            q_jit(mask_2film) = abs(q_jit(mask_2film) + ...
                sigma_q .* randn(sum(mask_2film), 1));
            q_jit(q_jit < 1e-4) = 1e-4;
        end
        jittered.q_abs_Ainv = q_jit;

        try
            fit = local_fit_thickness_model(jittered, datasets, ...
                epsilon_bg, 'lsq', base_fit.x_fit);
            row = local_fit_row(fit, sprintf('qjitter_%g_%03d', ...
                sigma_q, ri), 'qjitter', sigma_q);
            row.replicate = ri;
            row.q_jitter_sigma_Ainv = sigma_q;
            all_rows = [all_rows; row]; %#ok<AGROW>
        catch ME
            warning('run_b1_physical_fit_enhancements:QJitterFailed', ...
                'q-jitter replicate failed: %s', ME.message);
        end
    end
end

qjitter = local_qjitter_aggregate(all_rows, sigma_values);
end


function summary = local_qjitter_aggregate(rows, sigma_values)
summary = table();
for i = 1:numel(sigma_values)
    sigma_q = sigma_values(i);
    sub = rows(abs(rows.q_jitter_sigma_Ainv - sigma_q) < 1e-12, :);
    if isempty(sub)
        continue;
    end

    ratio = double(sub.thickness_ratio_r);
    rho1 = double(sub.rho0_1film_A);
    rho2 = double(sub.rho0_2film_A);
    row = table( ...
        sigma_q, height(sub), ...
        median(ratio, 'omitnan'), ...
        local_percentile(ratio, 2.5), ...
        local_percentile(ratio, 97.5), ...
        median(rho1, 'omitnan'), ...
        local_percentile(rho1, 2.5), ...
        local_percentile(rho1, 97.5), ...
        median(rho2, 'omitnan'), ...
        local_percentile(rho2, 2.5), ...
        local_percentile(rho2, 97.5), ...
        median(double(sub.RMSE_meV), 'omitnan'), ...
        'VariableNames', {'q_jitter_sigma_Ainv', 'n_replicates', ...
        'ratio_median', 'ratio_ci95_low', 'ratio_ci95_high', ...
        'rho0_1film_median_A', 'rho0_1film_ci95_low_A', ...
        'rho0_1film_ci95_high_A', 'rho0_2film_median_A', ...
        'rho0_2film_ci95_low_A', 'rho0_2film_ci95_high_A', ...
        'RMSE_median_meV'});
    summary = [summary; row]; %#ok<AGROW>
end
end


function huber = local_huber_fit_diagnostic(points, datasets, epsilon_bg)
fit_lsq = local_fit_thickness_model(points, datasets, epsilon_bg, 'lsq');
fit_huber = local_fit_thickness_model(points, datasets, epsilon_bg, ...
    'huber', fit_lsq.x_fit);

huber = [ ...
    local_fit_row(fit_lsq, 'least_squares', 'huber_comparison', NaN); ...
    local_fit_row(fit_huber, 'huber_delta_1p345', ...
    'huber_comparison', NaN)];
end


function [structure, residual_points] = local_residual_diagnostics( ...
    qabs_points, datasets, fit, epsilon_bg)
points = qabs_points(qabs_points.include_for_fit, :);
[E_pred, residuals] = local_predict(points, fit, epsilon_bg);
sigma = max(double(points.sigma_fit_meV), eps);

residual_points = points;
residual_points.E_pred_meV = E_pred;
residual_points.residual_meV = residuals;
residual_points.normalized_residual = residuals ./ sigma;

structure = table();
all_row = local_residual_row(residual_points, 'all_sessions');
structure = [structure; all_row]; %#ok<AGROW>

for i = 1:numel(datasets)
    sub = residual_points(residual_points.session_index == i, :);
    structure = [structure; local_residual_row(sub, ...
        datasets(i).session_key)]; %#ok<AGROW>
end
end


function row = local_residual_row(points, label)
q = double(points.q_abs_Ainv);
r = double(points.residual_meV);
sigma = double(points.sigma_fit_meV);
R2 = double(points.R2_mean_all);

valid = isfinite(q) & isfinite(r);
qv = q(valid);
rv = r(valid);
if numel(qv) >= 2
    p = polyfit(qv, rv, 1);
    slope = p(1);
    corr_q = local_corr(qv, rv);
else
    slope = NaN;
    corr_q = NaN;
end

corr_sigma = local_corr(abs(rv), sigma(valid));
corr_R2 = local_corr(abs(rv), R2(valid));

q1 = local_percentile(qv, 33.3);
q2 = local_percentile(qv, 66.7);
low = rv(qv <= q1);
mid = rv(qv > q1 & qv <= q2);
high = rv(qv > q2);

row = table( ...
    {label}, numel(rv), min(qv), max(qv), ...
    mean(rv, 'omitnan'), median(rv, 'omitnan'), ...
    sqrt(mean(rv .^ 2, 'omitnan')), slope, corr_q, ...
    corr_sigma, corr_R2, mean(low, 'omitnan'), ...
    mean(mid, 'omitnan'), mean(high, 'omitnan'), ...
    mean(high, 'omitnan') - mean(low, 'omitnan'), ...
    'VariableNames', {'diagnostic_group', 'n_points', ...
    'q_abs_min_Ainv', 'q_abs_max_Ainv', 'residual_mean_meV', ...
    'residual_median_meV', 'residual_RMSE_meV', ...
    'residual_vs_q_slope_meV_A', 'corr_residual_q', ...
    'corr_abs_residual_sigma', 'corr_abs_residual_R2', ...
    'low_q_residual_mean_meV', 'mid_q_residual_mean_meV', ...
    'high_q_residual_mean_meV', 'high_minus_low_residual_meV'});
end


function fit = local_fit_thickness_model(points, datasets, epsilon_bg, ...
    loss_type, x0_override)
if nargin < 5
    x0_override = [];
end

if isempty(x0_override)
    init = local_initial_from_single_sessions(points, datasets, epsilon_bg);
    x0 = [log(init.rho0_1film_A); ...
        log(init.thickness_ratio_r - 1); ...
        log(init.A_by_session(:))];
else
    x0 = x0_override(:);
end

cost_fn = @(x) local_fit_cost(x, points, datasets, epsilon_bg, loss_type);
opts = optimset('Display', 'off', 'MaxFunEvals', 20000, ...
    'MaxIter', 5000, 'TolFun', 1e-11, 'TolX', 1e-11);
x_fit = fminsearch(cost_fn, x0, opts);

params = local_params_from_x(x_fit, datasets);
[E_pred, residuals] = local_predict(points, params, epsilon_bg);
sigma = max(double(points.sigma_fit_meV), eps);
w = 1 ./ (sigma .^ 2);
z = residuals ./ sigma;

if strcmp(loss_type, 'huber')
    objective_value = sum(local_huber_rho(z, 1.345));
else
    objective_value = sum(z .^ 2);
end

n_points = height(points);
n_params = numel(x_fit);
dof = max(n_points - n_params, 1);
weighted_mean = sum(w .* points.energy_mean_meV) / sum(w);
SS_tot = sum(w .* (points.energy_mean_meV - weighted_mean) .^ 2);
SS_res = sum(w .* residuals .^ 2);

fit = struct();
fit.model_key = 'thickness_free_ratio';
fit.loss_type = loss_type;
fit.epsilon_bg = epsilon_bg;
fit.x_fit = x_fit(:);
fit.n_points = n_points;
fit.n_params = n_params;
fit.dof = dof;
fit.objective_value = objective_value;
fit.reduced_chi2 = sum(z .^ 2) / dof;
fit.R_squared = 1 - SS_res / max(SS_tot, eps);
fit.RMSE_meV = sqrt(mean(residuals .^ 2));
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

rho_grid = logspace(log10(0.05), log10(5000), 320);
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


function cost = local_fit_cost(x, points, datasets, epsilon_bg, loss_type)
if any(~isfinite(x)) || any(abs(x) > 80)
    cost = realmax('double') / 1e6;
    return;
end

params = local_params_from_x(x, datasets);
[~, residuals] = local_predict(points, params, epsilon_bg);
sigma = max(double(points.sigma_fit_meV), eps);
z = residuals ./ sigma;

if strcmp(loss_type, 'huber')
    cost = sum(local_huber_rho(z, 1.345));
else
    cost = sum(z .^ 2);
end

if ~isfinite(cost)
    cost = realmax('double') / 1e6;
end
end


function rho = local_huber_rho(z, delta)
az = abs(z);
rho = 0.5 .* z .^ 2;
mask = az > delta;
rho(mask) = delta .* (az(mask) - 0.5 .* delta);
end


function params = local_params_from_x(x, datasets)
n_sessions = numel(datasets);
thickness_factor = [datasets.thickness_factor]';

rho0_1film = exp(x(1));
r = 1 + exp(x(2));
rho0_2film = r * rho0_1film;
rho_by_session = zeros(n_sessions, 1);
rho_by_session(thickness_factor == 1) = rho0_1film;
rho_by_session(thickness_factor == 2) = rho0_2film;
A_by_session = exp(x(3:(2 + n_sessions)));

params = struct();
params.A_by_session = A_by_session(:);
params.rho_by_session_A = rho_by_session(:);
params.rho0_1film_A = rho0_1film;
params.rho0_2film_A = rho0_2film;
params.thickness_ratio_r = r;
end


function [E_pred, residuals] = local_predict(points, fit_or_params, ...
    epsilon_bg)
q = double(points.q_abs_Ainv);
session_index = double(points.session_index);
A = fit_or_params.A_by_session(session_index);
rho = fit_or_params.rho_by_session_A(session_index);
E_pred = sqrt(A .* q ./ (epsilon_bg + rho .* q));
residuals = double(points.energy_mean_meV) - E_pred;
end


function row = local_fit_row(fit, label, analysis_type, qmax_or_sigma)
row = table( ...
    {label}, {analysis_type}, {fit.loss_type}, fit.epsilon_bg, ...
    fit.n_points, fit.n_params, fit.dof, fit.objective_value, ...
    fit.reduced_chi2, fit.R_squared, fit.RMSE_meV, ...
    fit.rho0_1film_A, fit.rho0_2film_A, fit.thickness_ratio_r, ...
    fit.A_by_session(1), fit.A_by_session(2), fit.A_by_session(3), ...
    sqrt(fit.A_by_session(1) ./ fit.rho_by_session_A(1)), ...
    sqrt(fit.A_by_session(2) ./ fit.rho_by_session_A(2)), ...
    sqrt(fit.A_by_session(3) ./ fit.rho_by_session_A(3)), ...
    qmax_or_sigma, NaN, NaN, NaN, ...
    'VariableNames', {'fit_label', 'analysis_type', 'loss_type', ...
    'epsilon_bg', 'n_points', 'n_params', 'dof', 'objective_value', ...
    'reduced_chi2', 'R_squared', 'RMSE_meV', 'rho0_1film_A', ...
    'rho0_2film_A', 'thickness_ratio_r', 'A_590_PL2_10w', ...
    'A_n0_PL2_10w_repeat', 'A_no_PL2_20w_2film', ...
    'E_flat_590_PL2_10w_meV', 'E_flat_n0_PL2_10w_repeat_meV', ...
    'E_flat_no_PL2_20w_2film_meV', 'q_abs_max_Ainv', ...
    'q_abs_data_max_Ainv', 'replicate', 'q_jitter_sigma_Ainv'});
end


function p = local_percentile(values, pct)
values = sort(values(isfinite(values)));
if isempty(values)
    p = NaN;
    return;
end
if numel(values) == 1
    p = values(1);
    return;
end
pos = 1 + (pct / 100) * (numel(values) - 1);
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    p = values(lo);
else
    frac = pos - lo;
    p = values(lo) .* (1 - frac) + values(hi) .* frac;
end
end


function c = local_corr(x, y)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);
if numel(x) < 3 || std(x) == 0 || std(y) == 0
    c = NaN;
    return;
end
x = x - mean(x);
y = y - mean(y);
c = sum(x .* y) ./ sqrt(sum(x .^ 2) .* sum(y .^ 2));
end


function local_plot_bootstrap_ratio(samples, summary, png_path, pdf_path)
ratio = double(samples.thickness_ratio_r);
ratio = ratio(isfinite(ratio));
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 820, 520]);
ax = axes(fig);
histogram(ax, ratio, 28, 'FaceColor', [0.34, 0.45, 0.60], ...
    'EdgeColor', 'none');
hold(ax, 'on');
ratio_summary = summary(strcmp(summary.metric, 'thickness_ratio_r'), :);
if ~isempty(ratio_summary)
    xline(ax, ratio_summary.base_value, '-', 'Base', ...
        'Color', [0.10, 0.10, 0.10], 'LineWidth', 1.3);
    xline(ax, ratio_summary.ci95_low, '--', '95% CI', ...
        'Color', [0.45, 0.12, 0.12], 'LineWidth', 1.0);
    xline(ax, ratio_summary.ci95_high, '--', ...
        'Color', [0.45, 0.12, 0.12], 'LineWidth', 1.0);
end
hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
xlabel(ax, 'rho0_{2film}/rho0_{1film}');
ylabel(ax, 'Bootstrap count');
title(ax, 'B1 thickness ratio bootstrap');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
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


function local_plot_qrange_stability(qrange, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 920, 640]);
t = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

ax1 = nexttile(t, 1);
plot(ax1, qrange.q_abs_max_Ainv, qrange.thickness_ratio_r, '-o', ...
    'LineWidth', 1.5, 'MarkerFaceColor', [0.20, 0.45, 0.65]);
hold(ax1, 'on');
plot(ax1, [min(qrange.q_abs_max_Ainv), max(qrange.q_abs_max_Ainv)], ...
    [2, 2], '--', 'Color', [0.25, 0.25, 0.25]);
hold(ax1, 'off');
box(ax1, 'on');
grid(ax1, 'on');
ax1.GridAlpha = 0.16;
ylabel(ax1, 'rho0 ratio r');
title(ax1, 'qmax stability');

ax2 = nexttile(t, 2);
plot(ax2, qrange.q_abs_max_Ainv, qrange.RMSE_meV, '-s', ...
    'LineWidth', 1.5, 'MarkerFaceColor', [0.55, 0.40, 0.25]);
box(ax2, 'on');
grid(ax2, 'on');
ax2.GridAlpha = 0.16;
xlabel(ax2, 'qmax (A^{-1})');
ylabel(ax2, 'RMSE (meV)');

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_qjitter_sensitivity(qjitter, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 860, 560]);
ax = axes(fig);
x = qjitter.q_jitter_sigma_Ainv;
y = qjitter.ratio_median;
lo = y - qjitter.ratio_ci95_low;
hi = qjitter.ratio_ci95_high - y;
errorbar(ax, x, y, lo, hi, '-o', 'LineWidth', 1.4, ...
    'MarkerFaceColor', [0.55, 0.32, 0.24], ...
    'Color', [0.55, 0.32, 0.24], 'CapSize', 5);
hold(ax, 'on');
plot(ax, [min(x), max(x)], [2, 2], '--', 'Color', [0.25, 0.25, 0.25]);
hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
xlabel(ax, '2film q-jitter sigma (A^{-1})');
ylabel(ax, 'rho0_{2film}/rho0_{1film}');
title(ax, 'Defocus-related q-jitter sensitivity');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_huber_comparison(huber, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 880, 500]);
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');

labels = huber.fit_label;
ax1 = nexttile(t, 1);
bar(ax1, huber.thickness_ratio_r, 'FaceColor', [0.30, 0.48, 0.62]);
box(ax1, 'on');
grid(ax1, 'on');
ax1.GridAlpha = 0.16;
ax1.XTick = 1:height(huber);
ax1.XTickLabel = labels;
ax1.XTickLabelRotation = 15;
ylabel(ax1, 'rho0 ratio r');
title(ax1, 'Robustness of ratio');

ax2 = nexttile(t, 2);
bar(ax2, huber.RMSE_meV, 'FaceColor', [0.45, 0.45, 0.45]);
box(ax2, 'on');
grid(ax2, 'on');
ax2.GridAlpha = 0.16;
ax2.XTick = 1:height(huber);
ax2.XTickLabel = labels;
ax2.XTickLabelRotation = 15;
ylabel(ax2, 'RMSE (meV)');
title(ax2, 'LSQ vs Huber');

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_residual_structure(residual_points, datasets, ...
    png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1040, 560]);
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'loose');

ax1 = nexttile(t, 1);
hold(ax1, 'on');
plot(ax1, [0, max(residual_points.q_abs_Ainv)], [0, 0], '-', ...
    'Color', [0.25, 0.25, 0.25], 'HandleVisibility', 'off');
for i = 1:numel(datasets)
    mask = residual_points.session_index == i;
    scatter(ax1, residual_points.q_abs_Ainv(mask), ...
        residual_points.residual_meV(mask), 34, ...
        'MarkerFaceColor', datasets(i).color, ...
        'MarkerEdgeColor', datasets(i).color, ...
        'DisplayName', datasets(i).session_label);
end
hold(ax1, 'off');
box(ax1, 'on');
grid(ax1, 'on');
ax1.GridAlpha = 0.16;
xlabel(ax1, '|q| (A^{-1})');
ylabel(ax1, 'Residual (meV)');
title(ax1, 'Residual vs q');
legend(ax1, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off');

ax2 = nexttile(t, 2);
scatter(ax2, residual_points.sigma_fit_meV, ...
    abs(residual_points.normalized_residual), 34, ...
    residual_points.q_abs_Ainv, 'filled');
box(ax2, 'on');
grid(ax2, 'on');
ax2.GridAlpha = 0.16;
xlabel(ax2, 'sigma_E (meV)');
ylabel(ax2, '|normalized residual|');
title(ax2, 'Residual vs energy uncertainty');
cb = colorbar(ax2);
cb.Label.String = '|q| (A^{-1})';

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end
