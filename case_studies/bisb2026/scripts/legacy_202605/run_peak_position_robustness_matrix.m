function output = run_peak_position_robustness_matrix(options)
%RUN_PEAK_POSITION_ROBUSTNESS_MATRIX Build paired B1 peak-position checks.
%
% This diagnostic keeps the current Fano branch assignment as the reference
% and compares peak positions on the same q channels against:
%   1) Lorentzian refit with the same current q range;
%   2) local maximum in the current processed spectrum;
%   3) local maximum after a simple linear background removal;
%   4) local maximum without area normalization;
%   5) local maximum with two scalar area-normalization windows.

arguments
    options.dateTag {mustBeTextScalar} = "260520"
    options.reuseLorentz (1,1) logical = true
    options.localWindowHalfMeV (1,1) double = 180
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
run(fullfile(project_root, 'startup.m'));

date_tag = char(string(options.dateTag));
out_dir = fullfile(project_root, 'paper_results', ...
    ['peak_position_robustness_' date_tag]);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

sessions = local_sessions();
lorentz_suffix = ['_lorentz_robust' date_tag];
lorentz_outputs = local_load_or_run_lorentz(project_root, sessions, ...
    lorentz_suffix, options.reuseLorentz);

matrix_tbl = table();
for i = 1:numel(sessions)
    fano = local_load_output(project_root, sessions(i).fano_tag);
    lorentz = lorentz_outputs{i};
    matrix_tbl = [matrix_tbl; local_session_matrix( ...
        fano, lorentz, sessions(i), options.localWindowHalfMeV)]; %#ok<AGROW>
end

summary_tbl = local_summary(matrix_tbl);

matrix_path = fullfile(out_dir, 'peak_position_robustness_matrix.csv');
summary_path = fullfile(out_dir, 'peak_position_robustness_summary.csv');
writetable(matrix_tbl, matrix_path);
writetable(summary_tbl, summary_path);

fig_png = fullfile(out_dir, 'peak_position_robustness_summary.png');
fig_pdf = fullfile(out_dir, 'peak_position_robustness_summary.pdf');
local_plot_summary(summary_tbl, fig_png, fig_pdf);

report_path = fullfile(out_dir, 'peak_position_robustness_report.md');
local_write_report(report_path, matrix_tbl, summary_tbl, sessions, ...
    lorentz_suffix, options.localWindowHalfMeV, fig_png, fig_pdf);

output = struct();
output.output_dir = out_dir;
output.matrix = matrix_tbl;
output.summary = summary_tbl;
output.matrix_path = matrix_path;
output.summary_path = summary_path;
output.figure_png = fig_png;
output.figure_pdf = fig_pdf;
output.report_path = report_path;

fprintf('Peak-position robustness matrix complete.\n');
fprintf('  Output: %s\n', out_dir);
fprintf('  Rows: %d\n', height(matrix_tbl));
fprintf('  Summary rows: %d\n', height(summary_tbl));
end


function sessions = local_sessions()
sessions = struct( ...
    'key', {}, 'label', {}, 'fano_tag', {}, 'order', {});

sessions(1).key = '590_PL2_10w';
sessions(1).label = '10w 1film';
sessions(1).fano_tag = '590_gui_history_area_260506';
sessions(1).order = 1;

sessions(2).key = 'n0_PL2_10w_repeat';
sessions(2).label = '10w repeat';
sessions(2).fano_tag = 'n0_PL2_10w_gui_history_area_260506';
sessions(2).order = 2;

sessions(3).key = 'no_PL2_20w_2film';
sessions(3).label = '20w 2film';
sessions(3).fano_tag = 'no_PL2_20w_2film_gui_history_area_260506_highq_refined';
sessions(3).order = 3;
end


function outputs = local_load_or_run_lorentz(project_root, sessions, suffix, reuse_existing)
tags = cell(1, numel(sessions));
for i = 1:numel(sessions)
    tags{i} = [sessions(i).fano_tag suffix];
end
mat_paths = cellfun(@(tag) fullfile(project_root, 'paper_results', ...
    tag, 'analysis_results.mat'), tags, 'UniformOutput', false);
if reuse_existing && all(cellfun(@isfile, mat_paths))
    outputs = cellfun(@(p) local_load_output_by_path(p), mat_paths, ...
        'UniformOutput', false);
    return
end

run_out = run_590_gui_history_area_analysis("all", ...
    qRangeOverride_Ainv=[-0.015 0.015], ...
    outputTagSuffix=string(suffix), peakModelOverride="lorentz");
outputs = run_out.sessions;
end


function output = local_load_output(project_root, tag)
output = local_load_output_by_path(fullfile(project_root, 'paper_results', ...
    tag, 'analysis_results.mat'));
end


function output = local_load_output_by_path(path)
if ~isfile(path)
    error('run_peak_position_robustness_matrix:MissingOutput', ...
        'Missing analysis output: %s', path);
end
loaded = load(path, 'output');
output = loaded.output;
end


function tbl = local_session_matrix(fano, lorentz, session, window_half_mev)
if numel(fano.branches) < 1 || isempty(fano.branches{1})
    tbl = table();
    return
end

fano_b1 = fano.branches{1};
lorentz_b1 = local_branch_or_empty(lorentz.branches, 1);

qe_current = fano.qe_pp;
qe_no_norm = local_variant_qe(fano, "no_norm", NaN, NaN);
qe_norm_wide = local_variant_qe(fano, "area_50_3800", 50, 3800);
qe_norm_mid = local_variant_qe(fano, "area_500_3500", 500, 3500);

rows = local_empty_rows();
for i = 1:size(fano_b1, 1)
    q = fano_b1(i, 1);
    e_fano = fano_b1(i, 2);
    lorentz_row = local_match_q(lorentz_b1, q, max(fano.session.dq_Ainv, 1e-6));

    current = local_peak_at_q(qe_current, q, e_fano, window_half_mev, false);
    linear_bg = local_peak_at_q(qe_current, q, e_fano, window_half_mev, true);
    no_norm = local_peak_at_q(qe_no_norm, q, e_fano, window_half_mev, false);
    norm_wide = local_peak_at_q(qe_norm_wide, q, e_fano, window_half_mev, false);
    norm_mid = local_peak_at_q(qe_norm_mid, q, e_fano, window_half_mev, false);

    row = local_empty_row();
    row.dataset_key = char(session.key);
    row.dataset_label = char(session.label);
    row.dataset_order = session.order;
    row.branch = 1;
    row.q_Ainv = q;
    row.q_abs_Ainv = abs(q);
    row.fano_apex_meV = e_fano;
    row.fano_R2 = fano_b1(i, 4);
    row.fano_gamma_over_E = fano_b1(i, 3) ./ max(e_fano, eps);

    if ~isempty(lorentz_row)
        row.lorentz_energy_meV = lorentz_row(2);
        row.lorentz_delta_meV = lorentz_row(2) - e_fano;
        row.lorentz_R2 = lorentz_row(4);
        if abs(row.lorentz_delta_meV) <= window_half_mev
            row.lorentz_used_delta_meV = row.lorentz_delta_meV;
            row.lorentz_status = "paired";
        else
            row.lorentz_status = "paired_energy_outside_window";
        end
    else
        row.lorentz_status = "missing_q";
    end

    row.current_local_max_meV = current.energy_meV;
    row.current_local_delta_meV = current.energy_meV - e_fano;
    row.linear_bg_local_max_meV = linear_bg.energy_meV;
    row.linear_bg_delta_meV = linear_bg.energy_meV - e_fano;
    row.no_norm_local_max_meV = no_norm.energy_meV;
    row.no_norm_delta_meV = no_norm.energy_meV - e_fano;
    row.area_50_3800_local_max_meV = norm_wide.energy_meV;
    row.area_50_3800_delta_meV = norm_wide.energy_meV - e_fano;
    row.area_500_3500_local_max_meV = norm_mid.energy_meV;
    row.area_500_3500_delta_meV = norm_mid.energy_meV - e_fano;
    row.local_window_half_meV = window_half_mev;
    row.method_note = "local maxima are measured within +/- window around the accepted Fano apex";

    rows(end + 1) = row; %#ok<AGROW>
end
tbl = struct2table(rows);
end


function qe_variant = local_variant_qe(fano, variant_name, norm_min, norm_max)
opts = fano.preprocess_opts;
switch char(variant_name)
    case 'no_norm'
        opts.do_normalize = false;
        opts.do_bg_sub = false;
    otherwise
        opts.do_normalize = true;
        opts.norm_method = 'Area';
        opts.norm_min = norm_min;
        opts.norm_max = norm_max;
        opts.do_bg_sub = false;
end
qe_variant = qe_preprocess(fano.dataset.qe, opts);
end


function br = local_branch_or_empty(branches, idx)
if numel(branches) < idx || isempty(branches{idx})
    br = zeros(0, 12);
else
    br = branches{idx};
end
end


function row = local_match_q(branch, q, dq_tol)
row = [];
if isempty(branch)
    return
end
[delta, idx] = min(abs(branch(:, 1) - q));
if isfinite(delta) && delta <= max(dq_tol * 0.55, 1e-8)
    row = branch(idx, :);
end
end


function peak = local_peak_at_q(qe, q, center_mev, half_width_mev, subtract_linear_bg)
q_axis = double(qe.q_Ainv(:));
[q_delta, q_idx] = min(abs(q_axis - q));
if isempty(q_idx) || ~isfinite(q_delta)
    peak = struct('energy_meV', NaN, 'status', "missing_q");
    return
end

energy = double(qe.energy_meV(:));
spectrum = double(qe.intensity(:, q_idx));
mask = energy >= center_mev - half_width_mev & ...
    energy <= center_mev + half_width_mev & isfinite(spectrum);
if nnz(mask) < 5
    peak = struct('energy_meV', NaN, 'status', "empty_window");
    return
end

e_win = energy(mask);
y_win = spectrum(mask);
if subtract_linear_bg
    y_win = y_win - local_linear_baseline(e_win, y_win);
end

[~, idx] = max(y_win);
peak = struct('energy_meV', e_win(idx), 'status', "ok");
end


function bg = local_linear_baseline(x, y)
n = numel(x);
edge_n = max(2, round(n * 0.12));
idx = unique([1:edge_n, max(1, n - edge_n + 1):n]);
coef = polyfit(x(idx), y(idx), 1);
bg = polyval(coef, x);
end


function rows = local_empty_rows()
rows = repmat(local_empty_row(), 0, 1);
end


function row = local_empty_row()
row = struct();
row.dataset_key = '';
row.dataset_label = '';
row.dataset_order = NaN;
row.branch = NaN;
row.q_Ainv = NaN;
row.q_abs_Ainv = NaN;
row.fano_apex_meV = NaN;
row.fano_R2 = NaN;
row.fano_gamma_over_E = NaN;
row.lorentz_energy_meV = NaN;
row.lorentz_delta_meV = NaN;
row.lorentz_used_delta_meV = NaN;
row.lorentz_R2 = NaN;
row.lorentz_status = "";
row.current_local_max_meV = NaN;
row.current_local_delta_meV = NaN;
row.linear_bg_local_max_meV = NaN;
row.linear_bg_delta_meV = NaN;
row.no_norm_local_max_meV = NaN;
row.no_norm_delta_meV = NaN;
row.area_50_3800_local_max_meV = NaN;
row.area_50_3800_delta_meV = NaN;
row.area_500_3500_local_max_meV = NaN;
row.area_500_3500_delta_meV = NaN;
row.local_window_half_meV = NaN;
row.method_note = "";
end


function summary = local_summary(matrix_tbl)
methods = { ...
    'lorentz', 'lorentz_used_delta_meV'; ...
    'current local max', 'current_local_delta_meV'; ...
    'linear-bg local max', 'linear_bg_delta_meV'; ...
    'no-norm local max', 'no_norm_delta_meV'; ...
    'area 50-3800 local max', 'area_50_3800_delta_meV'; ...
    'area 500-3500 local max', 'area_500_3500_delta_meV'};

summary = table();
datasets = unique(matrix_tbl.dataset_key, 'stable');
for i = 1:numel(datasets)
    sub = matrix_tbl(strcmp(matrix_tbl.dataset_key, datasets{i}), :);
    for j = 1:size(methods, 1)
        method = methods{j, 1};
        col = methods{j, 2};
        vals = abs(double(sub.(col)));
        vals = vals(isfinite(vals));
        if isempty(vals)
            n = 0; med = NaN; p90 = NaN; mx = NaN;
        else
            n = numel(vals);
            med = median(vals, 'omitnan');
            p90 = local_percentile(vals, 90);
            mx = max(vals);
        end
        summary = [summary; table( ...
            {sub.dataset_key{1}}, {sub.dataset_label{1}}, sub.dataset_order(1), ...
            {method}, n, med, p90, mx, ...
            'VariableNames', {'dataset_key', 'dataset_label', ...
            'dataset_order', 'method', 'n_compared', ...
            'median_abs_delta_meV', 'p90_abs_delta_meV', ...
            'max_abs_delta_meV'})]; %#ok<AGROW>
    end
end
summary = sortrows(summary, {'dataset_order', 'method'});
end


function value = local_percentile(vals, pct)
vals = sort(vals(isfinite(vals(:))));
if isempty(vals)
    value = NaN;
    return
end
pos = 1 + (numel(vals) - 1) * pct / 100;
lo = floor(pos);
hi = ceil(pos);
if lo == hi
    value = vals(lo);
else
    value = vals(lo) + (vals(hi) - vals(lo)) * (pos - lo);
end
end


function local_plot_summary(summary_tbl, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 500]);
methods = unique(summary_tbl.method, 'stable');
datasets = unique(summary_tbl.dataset_label, 'stable');
values = NaN(numel(datasets), numel(methods));
for i = 1:numel(datasets)
    for j = 1:numel(methods)
        mask = strcmp(summary_tbl.dataset_label, datasets{i}) & ...
            strcmp(summary_tbl.method, methods{j});
        if any(mask)
            values(i, j) = summary_tbl.median_abs_delta_meV(find(mask, 1));
        end
    end
end

bar(values);
grid on;
box on;
ylabel('Median |Delta E| vs Fano apex (meV)');
set(gca, 'XTickLabel', datasets, 'XTickLabelRotation', 20);
legend(methods, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
    'FontSize', 8);
title('B1 peak-position robustness checks on paired q channels');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_write_report(report_path, matrix_tbl, summary_tbl, sessions, ...
    lorentz_suffix, window_half_mev, fig_png, fig_pdf)
fid = fopen(report_path, 'w');
if fid < 0
    error('run_peak_position_robustness_matrix:ReportOpenFailed', ...
        'Could not open report for writing: %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# B1 peak-position robustness matrix\n\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now, 31));
fprintf(fid, '## Scope\n\n');
fprintf(fid, ['This diagnostic keeps the current area-normalized no-BG Fano apex ' ...
    'B1 branch as the reference and checks paired q channels against Lorentzian refits, ' ...
    'local maxima, a local linear-background variant, and scalar area-normalization variants.\n\n']);
fprintf(fid, '- Lorentz output suffix: `%s`\n', lorentz_suffix);
fprintf(fid, '- Local maximum window: +/- %.1f meV around the accepted Fano apex\n', window_half_mev);
fprintf(fid, '- Lorentz pairs outside this same energy window are retained in the matrix but excluded from the summary statistics\n');
fprintf(fid, '- Branch: B1 only\n\n');

fprintf(fid, '## Inputs\n\n');
for i = 1:numel(sessions)
    fprintf(fid, '- %s: `paper_results/%s`\n', sessions(i).label, sessions(i).fano_tag);
end
fprintf(fid, '\n');

fprintf(fid, '## Summary\n\n');
fprintf(fid, '| Dataset | Method | N | median abs dE (meV) | p90 abs dE (meV) | max abs dE (meV) |\n');
fprintf(fid, '|---|---|---:|---:|---:|---:|\n');
for i = 1:height(summary_tbl)
    fprintf(fid, '| %s | %s | %d | %.1f | %.1f | %.1f |\n', ...
        summary_tbl.dataset_label{i}, summary_tbl.method{i}, ...
        summary_tbl.n_compared(i), summary_tbl.median_abs_delta_meV(i), ...
        summary_tbl.p90_abs_delta_meV(i), summary_tbl.max_abs_delta_meV(i));
end
fprintf(fid, '\n');

fprintf(fid, '## Evidence boundary\n\n');
fprintf(fid, ['- This is a paired robustness diagnostic for B1 peak positions, not a new branch-assignment pipeline.\n' ...
    '- Local maxima are intentionally measured near accepted Fano apex positions, so they test local peak-location stability rather than independently discovering branches.\n' ...
    '- Area normalization is a scalar operation per spectrum; it should not move a local maximum unless preprocessing interactions or numerical edge cases intervene.\n' ...
    '- The result can support the internal consistency of the current Fano-apex extraction only if the observed shifts remain small compared with the B1 energy rise and with the reported confidence intervals.\n']);
fprintf(fid, '\n## Outputs\n\n');
fprintf(fid, '- `peak_position_robustness_matrix.csv` (%d rows)\n', height(matrix_tbl));
fprintf(fid, '- `peak_position_robustness_summary.csv` (%d rows)\n', height(summary_tbl));
fprintf(fid, '- `%s`\n', fig_png);
fprintf(fid, '- `%s`\n', fig_pdf);
end
