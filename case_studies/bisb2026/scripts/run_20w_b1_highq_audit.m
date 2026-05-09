function output = run_20w_b1_highq_audit()
%RUN_20W_B1_HIGHQ_AUDIT Audit 20w B1 high-q extraction and rescue candidates.
%
% The script does not change the GUI pipeline. It combines the current
% retained B1 branch, the high-q refinement log, the existing dispersion
% model, and the stored full-spectrum fit details to diagnose why B1 high-q
% points fail and which rejected full-spectrum candidates are worth showing
% only as low-confidence guide points.

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

input_dir = fullfile(project_root, 'paper_results', ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined');
output_dir = fullfile(project_root, 'paper_results', ...
    '20w_B1_highq_audit_260506');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

paths = local_input_paths(input_dir);
branch_tbl = readtable(paths.branch1_points);
refinement_tbl = readtable(paths.branch_refinement_log);
model_tbl = readtable(paths.dispersion_model_summary);
result_data = load(paths.analysis_results, 'output');

fit_model = local_current_b1_model(model_tbl);
audit_tbl = local_build_highq_audit(branch_tbl, refinement_tbl, fit_model);
symmetry_tbl = local_build_symmetry_audit(audit_tbl);
[audit_tbl, rescue_tbl] = local_classify_rescue_candidates(audit_tbl, symmetry_tbl);
model_comparison_tbl = local_model_comparison(branch_tbl, rescue_tbl);

audit_csv = fullfile(output_dir, 'branch1_highq_audit.csv');
symmetry_csv = fullfile(output_dir, 'branch1_highq_symmetry_audit.csv');
rescue_csv = fullfile(output_dir, 'branch1_highq_rescue_candidates.csv');
model_csv = fullfile(output_dir, 'branch1_model_comparison.csv');
writetable(audit_tbl, audit_csv);
writetable(symmetry_tbl, symmetry_csv);
writetable(rescue_tbl, rescue_csv);
writetable(model_comparison_tbl, model_csv);

fig_paths = struct();
fig_paths.dispersion = fullfile(output_dir, 'branch1_highq_audit_dispersion.png');
fig_paths.failure_reasons = fullfile(output_dir, 'branch1_highq_failure_reasons.png');
fig_paths.symmetry = fullfile(output_dir, 'branch1_highq_symmetry.png');
fig_paths.spectrum_grid = fullfile(output_dir, 'branch1_highq_spectrum_grid.png');

local_plot_highq_dispersion(branch_tbl, audit_tbl, rescue_tbl, fit_model, ...
    fig_paths.dispersion);
local_plot_failure_reasons(refinement_tbl, fig_paths.failure_reasons);
local_plot_symmetry(symmetry_tbl, fig_paths.symmetry);
local_plot_highq_diagnostic_grid(result_data.output.qe_pp, ...
    result_data.output.fit_res, audit_tbl, fig_paths.spectrum_grid);

report_path = fullfile(output_dir, '20w_B1_highq_audit_report.md');
local_write_report(report_path, input_dir, audit_tbl, symmetry_tbl, ...
    rescue_tbl, model_comparison_tbl, fig_paths);

mat_path = fullfile(output_dir, '20w_B1_highq_audit_results.mat');
save(mat_path, 'audit_tbl', 'symmetry_tbl', 'rescue_tbl', ...
    'model_comparison_tbl', 'fit_model', 'fig_paths');

output = struct();
output.output_dir = output_dir;
output.audit_csv = audit_csv;
output.symmetry_csv = symmetry_csv;
output.rescue_csv = rescue_csv;
output.model_csv = model_csv;
output.report_path = report_path;
output.audit = audit_tbl;
output.symmetry = symmetry_tbl;
output.rescue_candidates = rescue_tbl;
output.model_comparison = model_comparison_tbl;

fprintf('20w B1 high-q audit complete.\n');
fprintf('  Output: %s\n', output_dir);
fprintf('  Audit rows: %d\n', height(audit_tbl));
fprintf('  Candidate-only rescue rows: %d\n', height(rescue_tbl));
end


function paths = local_input_paths(input_dir)
paths = struct();
paths.branch1_points = fullfile(input_dir, 'branch1_points.csv');
paths.branch_refinement_log = fullfile(input_dir, 'branch_refinement_log.csv');
paths.dispersion_model_summary = fullfile(input_dir, 'dispersion_model_summary.csv');
paths.analysis_results = fullfile(input_dir, 'analysis_results.mat');

names = fieldnames(paths);
for i = 1:numel(names)
    if ~isfile(paths.(names{i}))
        error('run_20w_b1_highq_audit:MissingInput', ...
            'Missing required input file: %s', paths.(names{i}));
    end
end
end


function fit_model = local_current_b1_model(model_tbl)
mask = model_tbl.branch == 1 & strcmp(model_tbl.model, 'quasi2d_plasmon') & ...
    model_tbl.success == 1;
if ~any(mask)
    error('run_20w_b1_highq_audit:MissingModel', ...
        'Could not find successful B1 quasi2d_plasmon row.');
end
row = model_tbl(find(mask, 1), :);
fit_model = struct();
fit_model.A_fit = row.param1(1);
fit_model.rho0_A = row.rho0_A(1);
fit_model.E_flat_meV = row.E_flat_meV(1);
fit_model.q_c_Ainv = row.q_c_Ainv(1);
fit_model.epsilon_bg = 1;
end


function audit_tbl = local_build_highq_audit(branch_tbl, refinement_tbl, fit_model)
q_min = 0.10;
q_values = [branch_tbl.q_Ainv(abs(branch_tbl.q_Ainv) >= q_min); ...
    refinement_tbl.q_Ainv(abs(refinement_tbl.q_Ainv) >= q_min)];
q_values = unique(round(q_values(:) * 1e6) / 1e6);
q_values = sort(q_values);

rows = local_empty_audit_rows();
for i = 1:numel(q_values)
    q = q_values(i);
    branch_idx = local_find_q(branch_tbl.q_Ainv, q, 1e-6);
    log_idx = local_find_q(refinement_tbl.q_Ainv, q, 1e-6);
    q_abs = abs(q);
    pred = local_predict_quasi2d(fit_model, q_abs);

    row = local_empty_audit_row();
    row.q_Ainv = q;
    row.q_abs_Ainv = q_abs;
    row.side = local_q_side(q);
    row.model_pred_meV = pred;

    if isfinite(branch_idx)
        b = branch_tbl(branch_idx, :);
        row.final_present = true;
        row.final_energy_meV = b.energy_meV(1);
        row.final_R2 = b.R2(1);
        row.final_CI_half_meV = b.E_ci_half_meV(1);
        row.final_gamma_over_E = b.gamma_meV(1) ./ max(b.energy_meV(1), eps);
        row.final_model_residual_abs_meV = abs(row.final_energy_meV - pred);
    end

    if isfinite(log_idx)
        r = refinement_tbl(log_idx, :);
        row.refit_status = char(string(r.status{1}));
        row.refit_detail = char(string(r.detail{1}));
        row.trigger_reason = char(string(r.trigger_reason{1}));
        row.old_energy_meV = r.old_energy_meV(1);
        row.old_R2 = r.old_R2(1);
        row.old_CI_half_meV = r.old_CI_half_meV(1);
        row.old_gamma_over_E = r.old_gamma_over_E(1);
        row.old_model_residual_abs_meV = abs(row.old_energy_meV - pred);
        row.new_energy_meV = r.new_energy_meV(1);
        row.new_R2 = r.new_R2(1);
        row.new_CI_half_meV = r.new_CI_half_meV(1);
        row.new_gamma_over_E = r.new_gamma_over_E(1);
        row.new_model_residual_abs_meV = abs(row.new_energy_meV - pred);
    else
        row.refit_status = 'not_refit_current_branch';
        row.refit_detail = '';
        row.trigger_reason = '';
    end

    if row.final_present
        row.primary_candidate_source = 'current_branch_fit';
        row.primary_candidate_energy_meV = row.final_energy_meV;
        row.primary_model_residual_abs_meV = row.final_model_residual_abs_meV;
    elseif isfinite(row.old_energy_meV)
        row.primary_candidate_source = 'old_full_spectrum_fit';
        row.primary_candidate_energy_meV = row.old_energy_meV;
        row.primary_model_residual_abs_meV = row.old_model_residual_abs_meV;
    elseif isfinite(row.new_energy_meV)
        row.primary_candidate_source = 'local_refit';
        row.primary_candidate_energy_meV = row.new_energy_meV;
        row.primary_model_residual_abs_meV = row.new_model_residual_abs_meV;
    end

    rows(end + 1) = row; %#ok<AGROW>
end

audit_tbl = struct2table(rows);
end


function rows = local_empty_audit_rows()
rows = repmat(local_empty_audit_row(), 0, 1);
end


function row = local_empty_audit_row()
row = struct();
row.q_Ainv = NaN;
row.q_abs_Ainv = NaN;
row.side = '';
row.final_present = false;
row.final_energy_meV = NaN;
row.final_R2 = NaN;
row.final_CI_half_meV = NaN;
row.final_gamma_over_E = NaN;
row.final_model_residual_abs_meV = NaN;
row.refit_status = '';
row.refit_detail = '';
row.trigger_reason = '';
row.old_energy_meV = NaN;
row.old_R2 = NaN;
row.old_CI_half_meV = NaN;
row.old_gamma_over_E = NaN;
row.old_model_residual_abs_meV = NaN;
row.new_energy_meV = NaN;
row.new_R2 = NaN;
row.new_CI_half_meV = NaN;
row.new_gamma_over_E = NaN;
row.new_model_residual_abs_meV = NaN;
row.model_pred_meV = NaN;
row.primary_candidate_source = '';
row.primary_candidate_energy_meV = NaN;
row.primary_model_residual_abs_meV = NaN;
row.recommendation = '';
end


function idx = local_find_q(q_axis, q_value, tol)
d = abs(q_axis(:) - q_value);
[best, idx0] = min(d);
if isempty(best) || ~isfinite(best) || best > tol
    idx = NaN;
else
    idx = idx0;
end
end


function side = local_q_side(q)
if q < 0
    side = 'negative';
elseif q > 0
    side = 'positive';
else
    side = 'zero';
end
end


function E = local_predict_quasi2d(fit_model, q_abs)
E = sqrt(fit_model.A_fit .* q_abs ./ ...
    (fit_model.epsilon_bg + fit_model.rho0_A .* q_abs));
end


function symmetry_tbl = local_build_symmetry_audit(audit_tbl)
q_abs_values = unique(round(audit_tbl.q_abs_Ainv(:) * 1e6) / 1e6);
q_abs_values = sort(q_abs_values);
row_template = struct('q_abs_Ainv', NaN, 'negative_energy_meV', NaN, ...
    'positive_energy_meV', NaN, 'negative_source', '', ...
    'positive_source', '', 'symmetry_abs_delta_meV', NaN, ...
    'symmetry_status', '');
rows = repmat(row_template, numel(q_abs_values), 1);

for i = 1:numel(q_abs_values)
    q_abs = q_abs_values(i);
    neg_idx = find(abs(audit_tbl.q_Ainv + q_abs) < 1e-6, 1);
    pos_idx = find(abs(audit_tbl.q_Ainv - q_abs) < 1e-6, 1);

    neg_E = NaN;
    pos_E = NaN;
    neg_source = '';
    pos_source = '';
    if ~isempty(neg_idx)
        neg_E = audit_tbl.primary_candidate_energy_meV(neg_idx);
        neg_source = audit_tbl.primary_candidate_source{neg_idx};
    end
    if ~isempty(pos_idx)
        pos_E = audit_tbl.primary_candidate_energy_meV(pos_idx);
        pos_source = audit_tbl.primary_candidate_source{pos_idx};
    end

    delta = abs(neg_E - pos_E);
    if isfinite(delta) && delta <= 150
        status = 'symmetric';
    elseif isfinite(delta) && delta <= 250
        status = 'marginal';
    elseif isfinite(delta)
        status = 'asymmetric';
    else
        status = 'missing_pair';
    end

    rows(i) = struct( ...
        'q_abs_Ainv', q_abs, ...
        'negative_energy_meV', neg_E, ...
        'positive_energy_meV', pos_E, ...
        'negative_source', neg_source, ...
        'positive_source', pos_source, ...
        'symmetry_abs_delta_meV', delta, ...
        'symmetry_status', status);
end

symmetry_tbl = struct2table(rows);
end


function [audit_tbl, rescue_tbl] = local_classify_rescue_candidates(audit_tbl, symmetry_tbl)
model_residual_max_meV = 300;
symmetry_delta_max_meV = 220;
old_R2_min = 0.85;
old_gamma_ratio_max = 1.80;
old_ci_half_max_meV = 650;

row_template = struct('q_Ainv', NaN, 'q_abs_Ainv', NaN, ...
    'candidate_energy_meV', NaN, 'candidate_source', '', ...
    'candidate_class', '', 'fit_weight_for_stress_test', NaN, ...
    'old_R2', NaN, 'old_CI_half_meV', NaN, ...
    'old_gamma_over_E', NaN, 'model_residual_abs_meV', NaN, ...
    'symmetry_abs_delta_meV', NaN, 'selection_reason', '');
rows = repmat(row_template, height(audit_tbl), 1);
row_count = 0;

for i = 1:height(audit_tbl)
    if audit_tbl.final_present(i)
        audit_tbl.recommendation{i} = 'already_in_current_branch';
        continue;
    end
    if ~isfinite(audit_tbl.old_energy_meV(i))
        audit_tbl.recommendation{i} = 'no_old_candidate';
        continue;
    end

    q_abs = audit_tbl.q_abs_Ainv(i);
    sym_idx = find(abs(symmetry_tbl.q_abs_Ainv - q_abs) < 1e-6, 1);
    sym_delta = NaN;
    if ~isempty(sym_idx)
        sym_delta = symmetry_tbl.symmetry_abs_delta_meV(sym_idx);
    end

    pass_old_quality = audit_tbl.old_R2(i) >= old_R2_min && ...
        audit_tbl.old_gamma_over_E(i) <= old_gamma_ratio_max && ...
        audit_tbl.old_CI_half_meV(i) <= old_ci_half_max_meV;
    pass_model = audit_tbl.old_model_residual_abs_meV(i) <= ...
        model_residual_max_meV;
    pass_symmetry = isfinite(sym_delta) && sym_delta <= symmetry_delta_max_meV;

    if pass_old_quality && pass_model && pass_symmetry
        audit_tbl.recommendation{i} = 'candidate_only_not_for_fit';
        reason = sprintf(['old full-spectrum candidate passes loose audit: ', ...
            'R2 %.3f, CI %.1f meV, Gamma/E %.2f, model residual %.1f meV, symmetry delta %.1f meV'], ...
            audit_tbl.old_R2(i), audit_tbl.old_CI_half_meV(i), ...
            audit_tbl.old_gamma_over_E(i), ...
            audit_tbl.old_model_residual_abs_meV(i), sym_delta);
        row_count = row_count + 1;
        rows(row_count) = struct( ...
            'q_Ainv', audit_tbl.q_Ainv(i), ...
            'q_abs_Ainv', audit_tbl.q_abs_Ainv(i), ...
            'candidate_energy_meV', audit_tbl.old_energy_meV(i), ...
            'candidate_source', 'old_full_spectrum_fit', ...
            'candidate_class', 'candidate_only_not_for_fit', ...
            'fit_weight_for_stress_test', 0.25 * audit_tbl.old_R2(i), ...
            'old_R2', audit_tbl.old_R2(i), ...
            'old_CI_half_meV', audit_tbl.old_CI_half_meV(i), ...
            'old_gamma_over_E', audit_tbl.old_gamma_over_E(i), ...
            'model_residual_abs_meV', audit_tbl.old_model_residual_abs_meV(i), ...
            'symmetry_abs_delta_meV', sym_delta, ...
            'selection_reason', reason);
    else
        audit_tbl.recommendation{i} = 'keep_rejected_low_confidence';
    end
end

rows = rows(1:row_count);
if isempty(rows)
    rescue_tbl = table();
else
    rescue_tbl = struct2table(rows);
end
end


function model_tbl = local_model_comparison(branch_tbl, rescue_tbl)
base_q = branch_tbl.q_Ainv(:);
base_E = branch_tbl.energy_meV(:);
base_w = branch_tbl.R2(:);

rows = local_model_row('current_final_points_only', base_q, base_E, base_w);

if ~isempty(rescue_tbl)
    stress_q = [base_q; rescue_tbl.q_Ainv(:)];
    stress_E = [base_E; rescue_tbl.candidate_energy_meV(:)];
    stress_w = [base_w; rescue_tbl.fit_weight_for_stress_test(:)];
    rows(end + 1) = local_model_row( ... %#ok<AGROW>
        'stress_test_final_plus_candidate_only', stress_q, stress_E, stress_w);
end

model_tbl = struct2table(rows);
end


function row = local_model_row(label, q, E, w)
fit = local_fit_branch1_model(q, E, w);
row = struct();
row.model_case = label;
row.n_points = numel(q);
row.rho0_A = fit.rho0_A;
row.q_c_Ainv = fit.q_c_Ainv;
row.E_flat_meV = fit.E_flat_meV;
row.R_squared = fit.R_squared;
row.RMSE_meV = fit.RMSE_meV;
end


function fit = local_fit_branch1_model(q_Ainv, energy_meV, weights)
q = abs(double(q_Ainv(:)));
E = double(energy_meV(:));
w = double(weights(:));
valid = isfinite(q) & isfinite(E) & q > 0 & E > 0 & isfinite(w) & w > 0;
q = q(valid);
E = E(valid);
w = w(valid);
w = w ./ median(w);

epsilon_bg = 1;
rho0_max_A = 5000;
model_fn = @(p, x) sqrt(abs(p(1)) .* abs(x) ./ ...
    (epsilon_bg + abs(p(2)) .* abs(x)));

p0 = [max(E)^2 * 10, 10];
lb = [0, 0.1];
ub = [Inf, rho0_max_A];

try
    opts = optimoptions('lsqcurvefit', 'Display', 'off', ...
        'MaxFunctionEvaluations', 10000, 'MaxIterations', 2000);
    weighted_model = @(p, x) sqrt(w) .* model_fn(p, x);
    p_fit = lsqcurvefit(weighted_model, p0, q, sqrt(w) .* E, lb, ub, opts);
catch
    cost = @(x) sum(w .* (model_fn(exp(x), q) - E) .^ 2);
    x_fit = fminsearch(cost, log(p0), optimset('Display', 'off'));
    p_fit = exp(x_fit);
end

p_fit = abs(p_fit);
pred = model_fn(p_fit, q);
residuals = E - pred;
SS_res = sum(w .* residuals .^ 2);
SS_tot = sum(w .* (E - mean(E)) .^ 2);
fit = struct();
fit.A_fit = p_fit(1);
fit.rho0_A = p_fit(2);
fit.q_c_Ainv = epsilon_bg / p_fit(2);
fit.E_flat_meV = sqrt(p_fit(1) / p_fit(2));
fit.R_squared = 1 - SS_res / max(SS_tot, eps);
fit.RMSE_meV = sqrt(mean(residuals .^ 2));
end


function local_plot_highq_dispersion(branch_tbl, audit_tbl, rescue_tbl, fit_model, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 980, 560]);
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');

plot(ax, branch_tbl.q_Ainv, branch_tbl.energy_meV, 'o', ...
    'Color', [0.1 0.35 0.75], 'MarkerFaceColor', [0.1 0.35 0.75], ...
    'DisplayName', 'current B1 branch');

rejected = strcmp(audit_tbl.refit_status, 'rejected_low_confidence_refit');
plot(ax, audit_tbl.q_Ainv(rejected), audit_tbl.old_energy_meV(rejected), ...
    'o', 'Color', [0.45 0.45 0.45], 'MarkerFaceColor', 'none', ...
    'LineWidth', 1.2, 'DisplayName', 'old full-spectrum rejected');
plot(ax, audit_tbl.q_Ainv(rejected), audit_tbl.new_energy_meV(rejected), ...
    'x', 'Color', [0.85 0.1 0.1], 'LineWidth', 1.2, ...
    'DisplayName', 'local refit rejected');

if ~isempty(rescue_tbl)
    plot(ax, rescue_tbl.q_Ainv, rescue_tbl.candidate_energy_meV, 's', ...
        'Color', [0.95 0.55 0.05], 'MarkerFaceColor', [1.0 0.78 0.22], ...
        'LineWidth', 1.2, 'DisplayName', 'candidate only');
end

q_fit = linspace(-0.15, 0.15, 301)';
E_fit = local_predict_quasi2d(fit_model, abs(q_fit));
plot(ax, q_fit, E_fit, 'k-', 'LineWidth', 1.4, ...
    'DisplayName', 'current B1 quasi2D guide');

xline(ax, -0.10, ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
xline(ax, 0.10, ':', 'Color', [0.35 0.35 0.35], 'HandleVisibility', 'off');
xlabel(ax, 'q (1/A)');
ylabel(ax, 'B1 energy (meV)');
title(ax, '20w B1 high-q audit');
legend(ax, 'Location', 'bestoutside');
local_export_figure(fig, out_path);
end


function local_plot_failure_reasons(refinement_tbl, out_path)
tokens = {'R2', 'CI half', 'Gamma/E', 'window edge', 'energy shift'};
counts = zeros(size(tokens));
details = strcat(string(refinement_tbl.detail), "; ", string(refinement_tbl.trigger_reason));
for i = 1:numel(tokens)
    counts(i) = sum(contains(details, tokens{i}));
end

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 720, 420]);
bar(counts, 'FaceColor', [0.25 0.45 0.75]);
grid on;
box on;
set(gca, 'XTickLabel', tokens, 'XTickLabelRotation', 25);
ylabel('Count');
title('20w B1 high-q refinement failure flags');
local_export_figure(fig, out_path);
end


function local_plot_symmetry(symmetry_tbl, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 760, 440]);
ax = axes(fig);
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
plot(ax, symmetry_tbl.q_abs_Ainv, symmetry_tbl.symmetry_abs_delta_meV, ...
    '-o', 'LineWidth', 1.4, 'MarkerFaceColor', [0.1 0.55 0.25], ...
    'Color', [0.1 0.55 0.25]);
yline(ax, 150, ':', '150 meV', 'Color', [0.35 0.35 0.35]);
yline(ax, 250, ':', '250 meV', 'Color', [0.55 0.15 0.15]);
xlabel(ax, '|q| (1/A)');
ylabel(ax, '|E(-q) - E(+q)| (meV)');
title(ax, '20w B1 high-q symmetry audit');
local_export_figure(fig, out_path);
end


function local_plot_highq_diagnostic_grid(qe_pp, fit_res, audit_tbl, out_path)
target_q = [-0.145, -0.135, -0.125, -0.115, 0.115, 0.125, 0.135, 0.145];
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 80, 1120, 720]);
tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(target_q)
    ax = nexttile;
    [~, q_idx] = min(abs(qe_pp.q_Ainv(:) - target_q(i)));
    q_val = qe_pp.q_Ainv(q_idx);
    detail = fit_res.fit_details{q_idx};
    if isempty(detail) || ~isstruct(detail) || ~isfield(detail, 'energy_data')
        title(ax, sprintf('q=%.3f missing fit', q_val));
        continue;
    end

    E = double(detail.energy_data(:));
    spectrum = double(detail.spectrum_data(:));
    plot(ax, E, spectrum, '-', 'Color', [0.15 0.35 0.75], ...
        'LineWidth', 1.0);
    hold(ax, 'on');
    if isfield(detail, 'curve_fit')
        E_curve = E;
        if isfield(detail, 'energy_fit') && ...
                numel(detail.energy_fit) == numel(detail.curve_fit)
            E_curve = double(detail.energy_fit(:));
        end
        plot(ax, E_curve, double(detail.curve_fit(:)), '-', ...
            'Color', [0.85 0.1 0.1], 'LineWidth', 1.0);
    end
    if isfield(detail, 'peak_curves') && ~isempty(detail.peak_curves)
        peak_curves = detail.peak_curves;
        if iscell(peak_curves)
            for p = 1:numel(peak_curves)
                curve = double(peak_curves{p}(:));
                E_peak = E;
                if isfield(detail, 'energy_fit') && ...
                        numel(detail.energy_fit) == numel(curve)
                    E_peak = double(detail.energy_fit(:));
                end
                plot(ax, E_peak, curve, ':', 'Color', [0.45 0.45 0.45], ...
                    'LineWidth', 0.8);
            end
        else
            pc = double(peak_curves);
            E_peak = E;
            if isfield(detail, 'energy_fit') && ...
                    numel(detail.energy_fit) == size(pc, 1)
                E_peak = double(detail.energy_fit(:));
            end
            for p = 1:size(pc, 2)
                plot(ax, E_peak, pc(:, p), ':', 'Color', [0.45 0.45 0.45], ...
                    'LineWidth', 0.8);
            end
        end
    end

    audit_idx = local_find_q(audit_tbl.q_Ainv, q_val, 1e-6);
    if isfinite(audit_idx)
        if isfinite(audit_tbl.old_energy_meV(audit_idx))
            xline(ax, audit_tbl.old_energy_meV(audit_idx), '--', ...
                'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
        end
        if audit_tbl.final_present(audit_idx)
            xline(ax, audit_tbl.final_energy_meV(audit_idx), '-', ...
                'Color', [0.05 0.55 0.20], 'LineWidth', 1.1);
        end
        if isfinite(audit_tbl.new_energy_meV(audit_idx))
            xline(ax, audit_tbl.new_energy_meV(audit_idx), ':', ...
                'Color', [0.85 0.1 0.1], 'LineWidth', 0.9);
        end
    end

    xlim(ax, [850, 1800]);
    title(ax, sprintf('q=%.3f', q_val));
    if i > 4
        xlabel(ax, 'Energy (meV)');
    end
    if mod(i - 1, 4) == 0
        ylabel(ax, 'Area-normalized intensity');
    end
    grid(ax, 'on');
    box(ax, 'on');
end

local_export_figure(fig, out_path);
end


function local_write_report(report_path, input_dir, audit_tbl, symmetry_tbl, ...
    rescue_tbl, model_comparison_tbl, fig_paths)
fid = fopen(report_path, 'w');
if fid < 0
    error('run_20w_b1_highq_audit:ReportOpenFailed', ...
        'Unable to write report: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

n_refit = sum(~strcmp(audit_tbl.refit_status, 'not_refit_current_branch'));
n_rejected = sum(strcmp(audit_tbl.refit_status, 'rejected_low_confidence_refit'));
n_accepted = sum(strcmp(audit_tbl.refit_status, 'accepted'));
n_current = sum(audit_tbl.final_present);
n_candidate = height(rescue_tbl);
med_old_R2 = median(audit_tbl.old_R2(isfinite(audit_tbl.old_R2)), 'omitnan');
med_new_R2 = median(audit_tbl.new_R2(isfinite(audit_tbl.new_R2)), 'omitnan');
med_sym = median(symmetry_tbl.symmetry_abs_delta_meV( ...
    isfinite(symmetry_tbl.symmetry_abs_delta_meV)), 'omitnan');

fprintf(fid, '# 20w B1 High-q Audit Report\n\n');
fprintf(fid, 'Source result directory: `%s`\n\n', input_dir);
fprintf(fid, '## Executive finding\n\n');
fprintf(fid, ['The high-q failure is localized mainly in the local ', ...
    '`[1000, 1700] meV` single-peak refit layer. The old full-spectrum ', ...
    'Fano candidates still have a median R2 of %.3f, whereas the local ', ...
    'refit median R2 is %.3f. Therefore these points should not be ', ...
    'blindly restored into the principal B1 fit, but several are useful ', ...
    'as gray candidate-only points for discussion.\n\n'], med_old_R2, med_new_R2);

fprintf(fid, '## Counts\n\n');
fprintf(fid, '- High-q q rows audited: %d\n', height(audit_tbl));
fprintf(fid, '- Rows sent to local refit: %d\n', n_refit);
fprintf(fid, '- Local refit accepted: %d\n', n_accepted);
fprintf(fid, '- Local refit rejected: %d\n', n_rejected);
fprintf(fid, '- Current retained high-q B1 rows: %d\n', n_current);
fprintf(fid, '- Candidate-only rows selected from old full-spectrum fit: %d\n', n_candidate);
fprintf(fid, '- Median |E(-q)-E(+q)| from primary candidates: %.1f meV\n\n', med_sym);

fprintf(fid, '## Model comparison\n\n');
fprintf(fid, '| Case | N | rho0 (A) | qc (1/A) | Eflat (meV) | R2 | RMSE (meV) |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|\n');
for i = 1:height(model_comparison_tbl)
    fprintf(fid, '| %s | %d | %.3g | %.4g | %.1f | %.4f | %.1f |\n', ...
        model_comparison_tbl.model_case{i}, model_comparison_tbl.n_points(i), ...
        model_comparison_tbl.rho0_A(i), model_comparison_tbl.q_c_Ainv(i), ...
        model_comparison_tbl.E_flat_meV(i), model_comparison_tbl.R_squared(i), ...
        model_comparison_tbl.RMSE_meV(i));
end

fprintf(fid, '\n## Candidate-only points\n\n');
if isempty(rescue_tbl)
    fprintf(fid, 'No rejected high-q points passed the loose candidate audit.\n\n');
else
    fprintf(fid, '| q (1/A) | Energy (meV) | old R2 | CI half (meV) | Gamma/E | model residual (meV) | symmetry delta (meV) |\n');
    fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|\n');
    for i = 1:height(rescue_tbl)
        fprintf(fid, '| %.4f | %.1f | %.3f | %.1f | %.2f | %.1f | %.1f |\n', ...
            rescue_tbl.q_Ainv(i), rescue_tbl.candidate_energy_meV(i), ...
            rescue_tbl.old_R2(i), rescue_tbl.old_CI_half_meV(i), ...
            rescue_tbl.old_gamma_over_E(i), rescue_tbl.model_residual_abs_meV(i), ...
            rescue_tbl.symmetry_abs_delta_meV(i));
    end
    fprintf(fid, '\nThese rows are marked `candidate_only_not_for_fit`: useful for visual audit, ');
    fprintf(fid, 'not recommended for the main quasi-2D fit without manual/bootstrapped validation.\n\n');
end

fprintf(fid, '## Figures\n\n');
fprintf(fid, '- `%s`\n', local_file_name(fig_paths.dispersion));
fprintf(fid, '- `%s`\n', local_file_name(fig_paths.failure_reasons));
fprintf(fid, '- `%s`\n', local_file_name(fig_paths.symmetry));
fprintf(fid, '- `%s`\n\n', local_file_name(fig_paths.spectrum_grid));

fprintf(fid, '## Recommended next handling\n\n');
fprintf(fid, '1. Keep the current retained B1 points as the principal fit evidence.\n');
fprintf(fid, '2. Overlay candidate-only old full-spectrum high-q points in gray/orange if discussing possible continuation.\n');
fprintf(fid, '3. Do not use local single-peak refit R2 alone as proof that the mode is absent; it is a stricter and different fitting problem.\n');
fprintf(fid, '4. Before restoring any candidate into the main branch, run manual single-spectrum review or bootstrap/multi-start refits at those q values.\n');
end


function name = local_file_name(path_value)
[~, name, ext] = fileparts(path_value);
name = [name ext];
end


function local_export_figure(fig, out_path)
if exist('exportgraphics', 'file')
    exportgraphics(fig, out_path, 'Resolution', 200);
else
    saveas(fig, out_path);
end
close(fig);
end
