function out = run_b1_lorentz_peak_evidence_audit(options)
%RUN_B1_LORENTZ_PEAK_EVIDENCE_AUDIT Audit spectral support for V15 B1 points.
%
% Fixed workflow guards for tests and future readers:
% energyWindowMeV=[300 1800], peakModel='lorentz', runFits=false.
% This is a diagnostic audit only; no physical fit was run.

arguments
    options.runPlots (1,1) logical = true
    options.maxPanelPoints (1,1) double = 24
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

workflow = struct();
workflow.summary_tag = 'b1_lorentz_tracking_peak_evidence_audit_260511';
workflow.v15_tag = '260510_lorentz_tracking_v15_upper_rapidrise_plateau';
workflow.energyWindowMeV = [300 1800];
workflow.peakModel = 'lorentz';
workflow.runFits = false;

summary_dir = fullfile(project_root, 'paper_results', workflow.summary_tag);
if ~isfolder(summary_dir)
    mkdir(summary_dir);
end

sessions = local_session_specs(project_root, workflow.v15_tag);
all_points = table();
all_release = table();
all_robustness = table();
all_competition = table();
figure_paths = table();

for si = 1:numel(sessions)
    fprintf('\nAuditing B1 evidence: %s\n', sessions(si).session_key);
    session_out = local_audit_one_session(sessions(si), workflow, summary_dir);
    all_points = [all_points; session_out.points]; %#ok<AGROW>
    all_release = [all_release; session_out.release]; %#ok<AGROW>
    all_robustness = [all_robustness; session_out.robustness]; %#ok<AGROW>
    all_competition = [all_competition; session_out.competition]; %#ok<AGROW>
    figure_paths = [figure_paths; session_out.figure_paths]; %#ok<AGROW>
end

summary = local_summary_table(all_points);
suspicious = all_points(strcmp(all_points.evidence_class, 'suspicious'), :);

points_csv = fullfile(summary_dir, 'b1_peak_evidence_audit_points.csv');
summary_csv = fullfile(summary_dir, 'b1_peak_evidence_audit_summary.csv');
release_csv = fullfile(summary_dir, 'b1_peak_evidence_constraint_release.csv');
robustness_csv = fullfile(summary_dir, 'b1_peak_evidence_robustness.csv');
competition_csv = fullfile(summary_dir, ...
    'b1_peak_evidence_candidate_competition.csv');
suspicious_csv = fullfile(summary_dir, ...
    'b1_peak_evidence_suspicious_points.csv');
writetable(all_points, points_csv);
writetable(summary, summary_csv);
writetable(all_release, release_csv);
writetable(all_robustness, robustness_csv);
writetable(all_competition, competition_csv);
writetable(suspicious, suspicious_csv);

summary_png = "";
summary_pdf = "";
panel_png = "";
panel_pdf = "";
if options.runPlots
    [summary_png, summary_pdf] = local_plot_three_dataset_summary( ...
        all_points, summary_dir);
    [panel_png, panel_pdf] = local_plot_suspicious_panels( ...
        all_points, summary_dir, options.maxPanelPoints);
end

readme_path = fullfile(summary_dir, ...
    'README_b1_lorentz_peak_evidence_audit.md');
local_write_readme(readme_path, workflow, summary, summary_png, panel_png);
local_update_manifest(project_root, summary_dir);
local_update_results_index(project_root, summary_dir, summary);

out = struct();
out.summary_dir = summary_dir;
out.points_csv = points_csv;
out.summary_csv = summary_csv;
out.release_csv = release_csv;
out.robustness_csv = robustness_csv;
out.competition_csv = competition_csv;
out.suspicious_csv = suspicious_csv;
out.summary_png = summary_png;
out.summary_pdf = summary_pdf;
out.suspicious_panels_png = panel_png;
out.suspicious_panels_pdf = panel_pdf;
out.figure_paths = figure_paths;

fprintf('\nB1 peak evidence audit complete.\n');
fprintf('  Summary: %s\n', summary_dir);
fprintf('  Points: %s\n', points_csv);
fprintf('  no physical fit was run\n');
end


function sessions = local_session_specs(project_root, v15_tag)
base = fullfile(project_root, 'paper_results');
sessions = struct('session_key', {}, 'session_label', {}, 'session_dir', {});
sessions(1).session_key = '590';
sessions(1).session_label = '590 10w defocus 1film';
sessions(1).session_dir = fullfile(base, ...
    ['590_gui_history_area_260506_b1_double_peak_binning_' v15_tag]);
sessions(2).session_key = 'n0';
sessions(2).session_label = 'n0 10w defocus repeat 1film';
sessions(2).session_dir = fullfile(base, ...
    ['n0_PL2_10w_gui_history_area_260506_b1_double_peak_binning_' v15_tag]);
sessions(3).session_key = '20w';
sessions(3).session_label = '20w defocus 2film';
sessions(3).session_dir = fullfile(base, ...
    ['no_PL2_20w_2film_gui_history_area_260506_highq_refined_b1_double_peak_binning_' v15_tag]);
end


function out = local_audit_one_session(spec, workflow, summary_dir)
if ~isfolder(spec.session_dir)
    error('run_b1_lorentz_peak_evidence_audit:MissingSessionDir', ...
        'Missing V15 session directory: %s', spec.session_dir);
end
mat_path = fullfile(spec.session_dir, 'b1_double_peak_binning_results.mat');
if ~isfile(mat_path)
    error('run_b1_lorentz_peak_evidence_audit:MissingResultsMat', ...
        'Missing V15 result MAT: %s', mat_path);
end
loaded = load(mat_path, 'extract', 'extract_opts', 'session', 'input_dir');
extract = loaded.extract;
extract_opts = loaded.extract_opts;
extract_opts.energy_window_meV = workflow.energyWindowMeV;
extract_opts.peak_model = workflow.peakModel;

lower = local_read_table(fullfile(spec.session_dir, ...
    'b1_double_peak_lower_points.csv'));
upper = local_read_table(fullfile(spec.session_dir, ...
    'b1_double_peak_upper_points.csv'));
candidate = local_read_table(fullfile(spec.session_dir, ...
    'b1_double_peak_lorentz_candidate_points.csv'));
path_selection = local_read_table(fullfile(spec.session_dir, ...
    'b1_double_peak_candidate_path_selection.csv'));

if isempty(lower) || isempty(upper)
    error('run_b1_lorentz_peak_evidence_audit:MissingPoints', ...
        'Missing lower/upper points in %s', spec.session_dir);
end

energy_axis = local_energy_axis(loaded, extract);
rows = table();
release_rows = table();
robust_rows = table();
competition_rows = table();

for ui = 1:height(extract.binning_map)
    detail = local_fit_detail(extract, ui);
    if isempty(detail)
        continue
    end
    bin_row = extract.binning_map(ui, :);
    q_value = double(bin_row.q_Ainv(1));
    lower_row = local_match_point(lower, q_value, 'b1_double_peak_lower');
    upper_row = local_match_point(upper, q_value, 'b1_double_peak_upper');
    if isempty(lower_row) || isempty(upper_row)
        continue
    end

    release = local_constraint_release(detail, extract_opts, ...
        lower_row.energy_meV(1), upper_row.energy_meV(1));
    robust = local_robustness_scan(detail, extract_opts, energy_axis, ...
        lower_row.energy_meV(1), upper_row.energy_meV(1));
    comp = local_candidate_competition(candidate, path_selection, q_value);

    [lower_audit, lower_release, lower_robust] = local_branch_audit_row( ...
        spec, bin_row, lower_row, upper_row, detail, release, robust, ...
        comp, 'lower');
    [upper_audit, upper_release, upper_robust] = local_branch_audit_row( ...
        spec, bin_row, lower_row, upper_row, detail, release, robust, ...
        comp, 'upper');
    rows = [rows; lower_audit; upper_audit]; %#ok<AGROW>
    release_rows = [release_rows; lower_release; upper_release]; %#ok<AGROW>
    robust_rows = [robust_rows; lower_robust; upper_robust]; %#ok<AGROW>
    competition_rows = [competition_rows; comp]; %#ok<AGROW>
end

rows = b1_peak_evidence_audit_classify(rows);
rows = local_attach_plot_cache(rows, spec, extract);

session_dir = fullfile(summary_dir, spec.session_key);
if ~isfolder(session_dir)
    mkdir(session_dir);
end
writetable(rows, fullfile(session_dir, ...
    'b1_peak_evidence_audit_points.csv'));
writetable(release_rows, fullfile(session_dir, ...
    'b1_peak_evidence_constraint_release.csv'));
writetable(robust_rows, fullfile(session_dir, ...
    'b1_peak_evidence_robustness.csv'));
writetable(competition_rows, fullfile(session_dir, ...
    'b1_peak_evidence_candidate_competition.csv'));

fig_paths = table();
map_png = "";
map_pdf = "";
if ~isempty(rows)
    [map_png, map_pdf] = local_plot_session_evidence_map( ...
        rows, session_dir, spec);
end
fig_paths = [fig_paths; table({spec.session_key}, {map_png}, {map_pdf}, ...
    'VariableNames', {'session_key', 'evidence_map_png', ...
    'evidence_map_pdf'})];

out = struct();
out.points = rows;
out.release = release_rows;
out.robustness = robust_rows;
out.competition = competition_rows;
out.figure_paths = fig_paths;
end


function tbl = local_read_table(path)
if isfile(path)
    tbl = readtable(path, 'TextType', 'string');
else
    tbl = table();
end
end


function detail = local_fit_detail(extract, idx)
detail = [];
if ~isfield(extract, 'fit_details') || idx > numel(extract.fit_details)
    return
end
detail = extract.fit_details{idx};
end


function energy = local_energy_axis(loaded, extract)
energy = [];
try
    if isfield(loaded, 'input_dir') && ...
            isfile(fullfile(loaded.input_dir, 'analysis_results.mat'))
        raw = load(fullfile(loaded.input_dir, 'analysis_results.mat'), 'output');
        if isfield(raw.output, 'qe_pp') && isfield(raw.output.qe_pp, 'energy_meV')
            energy = double(raw.output.qe_pp.energy_meV(:));
        end
    end
catch
    energy = [];
end
if isempty(energy) && isfield(extract, 'fit_details')
    for i = 1:numel(extract.fit_details)
        d = extract.fit_details{i};
        if ~isempty(d) && isfield(d, 'energy_data')
            energy = double(d.energy_data(:));
            return
        end
    end
end
end


function point = local_match_point(points, q_value, branch_label)
point = table();
if isempty(points)
    return
end
mask = abs(double(points.q_Ainv) - q_value) < 1e-9;
if ismember('branch_label', points.Properties.VariableNames)
    mask = mask & strcmp(string(points.branch_label), branch_label);
end
idx = find(mask, 1, 'first');
if ~isempty(idx)
    point = points(idx, :);
end
end


function release = local_constraint_release(detail, opts, lower_E, upper_E)
release = struct('success', false, 'status', "not_run", ...
    'lower_energy_meV', NaN, 'upper_energy_meV', NaN, ...
    'lower_delta_meV', NaN, 'upper_delta_meV', NaN, 'R2', NaN);
try
    energy = local_detail_energy(detail);
    spectrum = local_detail_spectrum(detail);
    fit = fit_loss_function(energy, spectrum, ...
        'E_min', opts.energy_window_meV(1), ...
        'E_max', opts.energy_window_meV(2), ...
        'max_peaks', 2, ...
        'min_prominence', 0.01, ...
        'smooth_width', 1, ...
        'initial_guesses', [lower_E; upper_E], ...
        'peak_model', 'lorentz', ...
        'pre_subtracted', opts.pre_subtracted, ...
        'bootstrap_ci_samples', 0, ...
        'min_peak_amplitude_fraction', 0);
    if fit.n_peaks ~= 2
        release.status = "not_two_peaks";
        return
    end
    E = sort(double(fit.apex_energy_meV(:)));
    if numel(E) < 2 || any(~isfinite(E(1:2)))
        E = sort(double(fit.omega_p(:)));
    end
    release.success = true;
    release.status = "ok";
    release.lower_energy_meV = E(1);
    release.upper_energy_meV = E(2);
    release.lower_delta_meV = abs(E(1) - lower_E);
    release.upper_delta_meV = abs(E(2) - upper_E);
    release.R2 = fit.R_squared;
catch ME
    release.status = string(ME.identifier);
end
end


function robust = local_robustness_scan(detail, opts, energy_axis, lower_E, upper_E)
robust = struct('lower_max_delta_meV', NaN, 'upper_max_delta_meV', NaN, ...
    'n_success', 0, 'n_attempted', 0, 'status', "not_run");
energy = energy_axis(:);
if isempty(energy) || ~isfield(detail, 'raw_unit_spectrum') || ...
        numel(detail.raw_unit_spectrum) ~= numel(energy)
    energy = local_detail_energy(detail);
    base_spectrum = local_detail_spectrum(detail);
else
    base_spectrum = double(detail.raw_unit_spectrum(:));
end
windows = unique([opts.energy_window_meV(2), 2000, 2100]);
windows = windows(windows >= opts.energy_window_meV(2));
sg_windows = local_robust_sg_windows(detail);
lower_delta = [];
upper_delta = [];
for wi = 1:numel(windows)
    for si = 1:numel(sg_windows)
        robust.n_attempted = robust.n_attempted + 1;
        y = local_sgolay_safe(base_spectrum, sg_windows(si), 3);
        try
            fit = fit_loss_function(energy, y, ...
                'E_min', opts.energy_window_meV(1), ...
                'E_max', windows(wi), ...
                'max_peaks', 2, ...
                'min_prominence', 0.01, ...
                'smooth_width', 1, ...
                'initial_guesses', [lower_E; upper_E], ...
                'peak_model', 'lorentz', ...
                'pre_subtracted', opts.pre_subtracted, ...
                'bootstrap_ci_samples', 0, ...
                'min_peak_amplitude_fraction', 0);
            if fit.n_peaks ~= 2
                continue
            end
            E = sort(double(fit.apex_energy_meV(:)));
            if numel(E) < 2 || any(~isfinite(E(1:2)))
                E = sort(double(fit.omega_p(:)));
            end
            lower_delta(end + 1, 1) = abs(E(1) - lower_E); %#ok<AGROW>
            upper_delta(end + 1, 1) = abs(E(2) - upper_E); %#ok<AGROW>
            robust.n_success = robust.n_success + 1;
        catch
        end
    end
end
if robust.n_success > 0
    robust.lower_max_delta_meV = max(lower_delta, [], 'omitnan');
    robust.upper_max_delta_meV = max(upper_delta, [], 'omitnan');
    robust.status = "ok";
else
    robust.status = "no_successful_refit";
end
end


function windows = local_robust_sg_windows(detail)
current = 91;
if isfield(detail, 'fit_denoise_window') && ...
        isfinite(double(detail.fit_denoise_window))
    current = double(detail.fit_denoise_window);
end
windows = unique(arrayfun(@local_odd_window, max(3, current + [-20 0 20])));
end


function w = local_odd_window(w)
w = round(w);
if mod(w, 2) == 0
    w = w + 1;
end
w = max(3, w);
end


function y = local_sgolay_safe(y, window, order)
y = double(y(:));
window = min(local_odd_window(window), numel(y));
if window <= order || window < 3
    return
end
try
    y = sgolayfilt(y, order, window);
catch
    y = smoothdata(y, 'sgolay', window);
end
end


function comp = local_candidate_competition(candidates, path_selection, q_value)
comp = table();
if isempty(candidates)
    return
end
qmask = abs(double(candidates.q_Ainv) - q_value) < 1e-9;
q_candidates = candidates(qmask, :);
if isempty(q_candidates)
    return
end
[best_R2, best_idx] = max(double(q_candidates.R2), [], 'omitnan');
selected_id = NaN;
selected_source = "";
selected_R2 = NaN;
selected_idx = [];
if ~isempty(path_selection)
    pmask = abs(double(path_selection.q_Ainv) - q_value) < 1e-9;
    pidx = find(pmask, 1, 'first');
    if ~isempty(pidx)
        selected_id = double(path_selection.selected_candidate_id(pidx));
        selected_source = string(path_selection.candidate_source(pidx));
        selected_idx = find(double(q_candidates.candidate_id) == ...
            selected_id, 1, 'first');
    end
end
if isempty(selected_idx)
    selected_idx = best_idx;
    selected_id = double(q_candidates.candidate_id(selected_idx));
    selected_source = string(q_candidates.candidate_source(selected_idx));
end
selected_R2 = double(q_candidates.R2(selected_idx));
comp = table(q_value, selected_id, selected_source, selected_R2, ...
    double(q_candidates.candidate_id(best_idx)), ...
    string(q_candidates.candidate_source(best_idx)), best_R2, ...
    best_R2 - selected_R2, height(q_candidates), ...
    'VariableNames', {'q_Ainv', 'selected_candidate_id', ...
    'selected_candidate_source', 'selected_R2', 'best_R2_candidate_id', ...
    'best_R2_candidate_source', 'best_R2', 'delta_R2_best_minus_selected', ...
    'n_candidates'});
end


function [audit_row, release_row, robust_row] = local_branch_audit_row( ...
    spec, bin_row, lower_row, upper_row, detail, release, robust, comp, branch)
if strcmp(branch, 'lower')
    point = lower_row;
    peak_idx = local_detail_peak_index(detail, point.energy_meV(1), 1);
    release_delta = release.lower_delta_meV;
    robust_delta = robust.lower_max_delta_meV;
    release_energy = release.lower_energy_meV;
else
    point = upper_row;
    peak_idx = local_detail_peak_index(detail, point.energy_meV(1), 2);
    release_delta = release.upper_delta_meV;
    robust_delta = robust.upper_max_delta_meV;
    release_energy = release.upper_energy_meV;
end

support = local_local_support(detail, point.energy_meV(1), point.gamma_meV(1));
ablation = local_component_ablation(detail, peak_idx, point.energy_meV(1), ...
    point.gamma_meV(1));
gamma_over_E = double(point.gamma_meV(1)) ./ ...
    max(abs(double(point.energy_meV(1))), eps);

selected_source = "";
delta_R2 = NaN;
if ~isempty(comp)
    selected_source = string(comp.selected_candidate_source(1));
    delta_R2 = double(comp.delta_R2_best_minus_selected(1));
end

audit_row = table({spec.session_key}, {spec.session_label}, ...
    point.q_Ainv(1), point.q_abs_Ainv(1), {char(point.branch_label(1))}, ...
    point.energy_meV(1), point.gamma_meV(1), gamma_over_E, ...
    point.R2(1), point.amplitude_fit(1), point.source_q_count(1), ...
    {char(point.source_mode(1))}, {char(point.source_q_Ainv(1))}, ...
    {char(point.source_q_index(1))}, point.bin_size_requested(1), ...
    {char(point.peak_model(1))}, {char(point.tracking_mode(1))}, ...
    {char(selected_source)}, delta_R2, support.score, ...
    support.local_snr, support.nearest_peak_delta_meV, ...
    support.has_local_peak_support, support.has_local_shoulder_support, ...
    ablation.sse_increase_fraction, release_delta, robust_delta, ...
    release.R2, release.success, {char(release.status)}, ...
    robust.n_success, robust.n_attempted, {char(robust.status)}, ...
    'VariableNames', {'session_key', 'session_label', 'q_Ainv', ...
    'q_abs_Ainv', 'branch_label', 'energy_meV', 'gamma_meV', ...
    'gamma_over_E', 'R2', 'amplitude_fit', 'source_q_count', ...
    'source_mode', 'source_q_Ainv', 'source_q_index', ...
    'bin_size_requested', 'peak_model', 'tracking_mode', ...
    'selected_candidate_source', 'delta_R2_best_minus_selected', ...
    'local_support_score', 'local_snr', 'nearest_peak_delta_meV', ...
    'has_local_peak_support', 'has_local_shoulder_support', ...
    'component_sse_increase_fraction', 'constraint_release_delta_meV', ...
    'robustness_max_delta_meV', 'constraint_release_R2', ...
    'constraint_release_success', 'constraint_release_status', ...
    'robustness_n_success', 'robustness_n_attempted', ...
    'robustness_status'});

release_row = table({spec.session_key}, point.q_Ainv(1), ...
    {char(point.branch_label(1))}, point.energy_meV(1), release_energy, ...
    release_delta, release.R2, release.success, {char(release.status)}, ...
    'VariableNames', {'session_key', 'q_Ainv', 'branch_label', ...
    'selected_energy_meV', 'release_energy_meV', ...
    'release_delta_meV', 'release_R2', 'release_success', ...
    'release_status'});

robust_row = table({spec.session_key}, point.q_Ainv(1), ...
    {char(point.branch_label(1))}, point.energy_meV(1), robust_delta, ...
    robust.n_success, robust.n_attempted, {char(robust.status)}, ...
    'VariableNames', {'session_key', 'q_Ainv', 'branch_label', ...
    'selected_energy_meV', 'robustness_max_delta_meV', ...
    'robustness_n_success', 'robustness_n_attempted', ...
    'robustness_status'});

% Keep plot-ready cached arrays inside UserData-like text columns after CSV
% writing is done by local_attach_plot_cache.
end


function energy = local_detail_energy(detail)
if isfield(detail, 'energy_data')
    energy = double(detail.energy_data(:));
elseif isfield(detail, 'energy_fit')
    energy = double(detail.energy_fit(:));
else
    energy = (1:numel(local_detail_spectrum(detail))).';
end
end


function spectrum = local_detail_spectrum(detail)
if isfield(detail, 'spectrum_data')
    spectrum = double(detail.spectrum_data(:));
elseif isfield(detail, 'fit_input_spectrum')
    spectrum = double(detail.fit_input_spectrum(:));
else
    spectrum = zeros(0, 1);
end
end


function idx = local_detail_peak_index(detail, energy, fallback_idx)
idx = fallback_idx;
if isfield(detail, 'apex_energy_meV')
    peaks = double(detail.apex_energy_meV(:));
elseif isfield(detail, 'omega_p')
    peaks = double(detail.omega_p(:));
else
    return
end
if isempty(peaks)
    return
end
[~, idx] = min(abs(peaks - energy));
end


function support = local_local_support(detail, peak_E, gamma)
energy = local_detail_energy(detail);
y = local_detail_spectrum(detail);
support = struct('score', 0, 'local_snr', NaN, ...
    'nearest_peak_delta_meV', NaN, 'has_local_peak_support', false, ...
    'has_local_shoulder_support', false);
if isempty(energy) || isempty(y)
    return
end
half_width = max(80, min(180, 0.35 * abs(gamma)));
local_mask = energy >= peak_E - half_width & energy <= peak_E + half_width;
if nnz(local_mask) < 5
    local_mask = abs(energy - peak_E) <= 120;
end
if nnz(local_mask) < 5
    return
end
e = energy(local_mask);
v = y(local_mask);
v_s = smoothdata(v, 'movmean', max(3, min(9, numel(v))));
[~, nearest_idx] = min(abs(e - peak_E));
noise = local_robust_noise(v_s);
local_base = median(v_s, 'omitnan');
peak_value = v_s(nearest_idx);
support.local_snr = (peak_value - local_base) ./ max(noise, eps);
maxima = local_local_maxima(e, v_s);
if ~isempty(maxima)
    [delta, ~] = min(abs(maxima - peak_E));
    support.nearest_peak_delta_meV = delta;
else
    support.nearest_peak_delta_meV = Inf;
end
support.has_local_peak_support = support.nearest_peak_delta_meV <= 80 && ...
    support.local_snr >= 3;
support.has_local_shoulder_support = support.local_snr >= 1.5 && ...
    local_has_shoulder(v_s, nearest_idx);
snr_score = min(max(support.local_snr ./ 5, 0), 1);
peak_score = double(support.has_local_peak_support);
shoulder_score = double(support.has_local_shoulder_support);
support.score = min(1, 0.55 * snr_score + 0.30 * peak_score + ...
    0.15 * shoulder_score);
end


function noise = local_robust_noise(y)
dy = diff(double(y(:)));
dy = dy(isfinite(dy));
if isempty(dy)
    noise = NaN;
    return
end
med = median(dy);
noise = median(abs(dy - med)) ./ 0.954;
if ~isfinite(noise) || noise <= 0
    noise = std(dy, 'omitnan');
end
if ~isfinite(noise) || noise <= 0
    noise = eps;
end
end


function locs = local_local_maxima(e, y)
locs = [];
if numel(y) < 3
    return
end
dy1 = diff(y);
turn = [false; dy1(1:end-1) > 0 & dy1(2:end) <= 0; false];
idx = find(turn);
if isempty(idx)
    [~, idx] = max(y);
end
locs = e(idx);
end


function tf = local_has_shoulder(y, idx)
idx = max(2, min(numel(y) - 1, idx));
left = y(idx) - y(max(1, idx - 2));
right = y(idx) - y(min(numel(y), idx + 2));
tf = (left >= 0 && right >= -0.05 * max(abs(y))) || ...
    (left >= -0.05 * max(abs(y)) && right >= 0);
end


function ablation = local_component_ablation(detail, peak_idx, peak_E, gamma)
ablation = struct('sse_increase_fraction', NaN);
if ~isfield(detail, 'energy_fit') || ~isfield(detail, 'curve_fit') || ...
        ~isfield(detail, 'peak_curves') || peak_idx > numel(detail.peak_curves)
    return
end
energy = local_detail_energy(detail);
y = local_detail_spectrum(detail);
fit_energy = double(detail.energy_fit(:));
total = interp1(fit_energy, double(detail.curve_fit(:)), energy, ...
    'linear', 'extrap');
component = interp1(fit_energy, double(detail.peak_curves{peak_idx}(:)), ...
    energy, 'linear', 'extrap');
half_width = max(120, min(260, 0.45 * abs(gamma)));
mask = energy >= peak_E - half_width & energy <= peak_E + half_width;
if nnz(mask) < 5
    return
end
res_full = y(mask) - total(mask);
res_without = y(mask) - (total(mask) - component(mask));
local_var = sum((y(mask) - mean(y(mask), 'omitnan')).^2, 'omitnan');
sse_full = sum(res_full .^ 2, 'omitnan');
sse_without = sum(res_without .^ 2, 'omitnan');
ablation.sse_increase_fraction = max(0, sse_without - sse_full) ./ ...
    max(local_var, eps);
end


function rows = local_attach_plot_cache(rows, spec, extract)
if isempty(rows)
    return
end
rows.plot_session_dir = repmat({spec.session_dir}, height(rows), 1);
rows.plot_detail_index = NaN(height(rows), 1);
for ui = 1:height(extract.binning_map)
    q = extract.binning_map.q_Ainv(ui);
    mask = abs(rows.q_Ainv - q) < 1e-9;
    rows.plot_detail_index(mask) = ui;
end
end


function [png_path, pdf_path] = local_plot_session_evidence_map(rows, out_dir, spec)
png_path = fullfile(out_dir, 'b1_peak_evidence_map.png');
pdf_path = fullfile(out_dir, 'b1_peak_evidence_map.pdf');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 980 560]);
ax = axes(fig);
hold(ax, 'on');
classes = {'data_supported', 'tracking_assisted', 'suspicious'};
colors = local_class_colors();
markers = containers.Map({'b1_double_peak_lower', 'b1_double_peak_upper'}, ...
    {'o', '^'});
for ci = 1:numel(classes)
    cmask = strcmp(rows.evidence_class, classes{ci});
    for bi = 1:2
        branch = {'b1_double_peak_lower', 'b1_double_peak_upper'};
        bmask = cmask & strcmp(rows.branch_label, branch{bi});
        if any(bmask)
            scatter(ax, rows.q_Ainv(bmask), rows.energy_meV(bmask), 36, ...
                colors(ci, :), markers(branch{bi}), 'filled', ...
                'DisplayName', sprintf('%s %s', classes{ci}, branch{bi}));
        end
    end
end
hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
xlabel(ax, 'q (A^{-1})');
ylabel(ax, 'Energy loss (meV)');
title(ax, sprintf('B1 peak evidence audit: %s', spec.session_label), ...
    'Interpreter', 'none');
legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function [png_path, pdf_path] = local_plot_three_dataset_summary(rows, summary_dir)
png_path = fullfile(summary_dir, 'b1_peak_evidence_three_dataset_summary.png');
pdf_path = fullfile(summary_dir, 'b1_peak_evidence_three_dataset_summary.pdf');
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 500]);
t = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
sessions = unique(rows.session_key, 'stable');
colors = local_class_colors();
classes = {'data_supported', 'tracking_assisted', 'suspicious'};
for si = 1:numel(sessions)
    ax = nexttile(t, si);
    hold(ax, 'on');
    smask = strcmp(rows.session_key, sessions{si});
    for ci = 1:numel(classes)
        cmask = smask & strcmp(rows.evidence_class, classes{ci});
        scatter(ax, rows.q_Ainv(cmask), rows.energy_meV(cmask), 28, ...
            colors(ci, :), 'filled', 'DisplayName', classes{ci});
    end
    hold(ax, 'off');
    box(ax, 'on');
    grid(ax, 'on');
    xlabel(ax, 'q (A^{-1})');
    ylabel(ax, 'Energy loss (meV)');
    title(ax, char(rows.session_label(find(smask, 1, 'first'))), ...
        'Interpreter', 'none');
end
legend(nexttile(t, 3), 'Location', 'eastoutside', 'Interpreter', 'none');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function colors = local_class_colors()
colors = [0.12 0.55 0.20; 0.15 0.35 0.85; 0.85 0.20 0.15];
end


function [png_path, pdf_path] = local_plot_suspicious_panels(rows, summary_dir, max_points)
png_path = fullfile(summary_dir, 'b1_peak_evidence_suspicious_panels.png');
pdf_path = fullfile(summary_dir, 'b1_peak_evidence_suspicious_panels.pdf');
suspicious = rows(strcmp(rows.evidence_class, 'suspicious'), :);
if isempty(suspicious)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 720 360]);
    ax = axes(fig);
    text(ax, 0.5, 0.5, 'No suspicious points', ...
        'HorizontalAlignment', 'center');
    axis(ax, 'off');
    exportgraphics(fig, png_path, 'Resolution', 300);
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    close(fig);
    return
end
suspicious = sortrows(suspicious, {'session_key', 'q_Ainv', 'branch_label'});
n = min(height(suspicious), max_points);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1500 900]);
t = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:min(n, 12)
    ax = nexttile(t, i);
    local_plot_single_panel(ax, suspicious(i, :));
end
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_single_panel(ax, row)
loaded = load(fullfile(row.plot_session_dir{1}, ...
    'b1_double_peak_binning_results.mat'), 'extract');
idx = row.plot_detail_index(1);
detail = loaded.extract.fit_details{idx};
energy = local_detail_energy(detail);
y = local_detail_spectrum(detail);
plot(ax, energy, y, 'k-', 'LineWidth', 0.8, 'DisplayName', 'spectrum');
hold(ax, 'on');
if isfield(detail, 'energy_fit') && isfield(detail, 'curve_fit')
    plot(ax, detail.energy_fit, detail.curve_fit, 'b-', 'LineWidth', 1.0, ...
        'DisplayName', 'full fit');
end
peak_idx = local_detail_peak_index(detail, row.energy_meV(1), 1);
if isfield(detail, 'peak_curves') && peak_idx <= numel(detail.peak_curves)
    total = interp1(detail.energy_fit, detail.curve_fit, energy, ...
        'linear', 'extrap');
    component = interp1(detail.energy_fit, detail.peak_curves{peak_idx}, ...
        energy, 'linear', 'extrap');
    plot(ax, energy, total - component, 'r--', 'LineWidth', 1.0, ...
        'DisplayName', 'without component');
end
xline(ax, row.energy_meV(1), 'Color', [0.85 0.20 0.15], ...
    'LineWidth', 1.1);
hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
xlim(ax, [max(300, row.energy_meV(1) - 300), ...
    min(1800, row.energy_meV(1) + 300)]);
title(ax, sprintf('%s q=%.4g %s', row.session_key{1}, ...
    row.q_Ainv(1), row.branch_label{1}), 'Interpreter', 'none');
xlabel(ax, 'Energy (meV)');
ylabel(ax, 'Intensity');
end


function summary = local_summary_table(points)
if isempty(points)
    summary = table();
    return
end
sessions = unique(points.session_key, 'stable');
classes = {'data_supported', 'tracking_assisted', 'suspicious'};
summary = table();
for si = 1:numel(sessions)
    smask = strcmp(points.session_key, sessions{si});
    n_total = sum(smask);
    counts = zeros(1, numel(classes));
    for ci = 1:numel(classes)
        counts(ci) = sum(smask & strcmp(points.evidence_class, classes{ci}));
    end
    row = table({sessions{si}}, n_total, counts(1), counts(2), counts(3), ...
        counts(3) ./ max(n_total, 1), ...
        'VariableNames', {'session_key', 'n_points', ...
        'n_data_supported', 'n_tracking_assisted', 'n_suspicious', ...
        'suspicious_fraction'});
    summary = [summary; row]; %#ok<AGROW>
end
total = table({'all'}, height(points), ...
    sum(strcmp(points.evidence_class, 'data_supported')), ...
    sum(strcmp(points.evidence_class, 'tracking_assisted')), ...
    sum(strcmp(points.evidence_class, 'suspicious')), ...
    sum(strcmp(points.evidence_class, 'suspicious')) ./ max(height(points), 1), ...
    'VariableNames', summary.Properties.VariableNames);
summary = [summary; total];
end


function local_write_readme(path, workflow, summary, summary_png, panel_png)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# B1 Lorentz Peak Evidence Audit\n\n');
fprintf(fid, '- Source tag: `%s`\n', workflow.v15_tag);
fprintf(fid, '- Energy window: `%g-%g meV`\n', workflow.energyWindowMeV);
fprintf(fid, '- Peak model: `%s`\n', workflow.peakModel);
fprintf(fid, '- Physical fit: no physical fit was run.\n\n');
fprintf(fid, '## Outputs\n\n');
fprintf(fid, '- `b1_peak_evidence_audit_points.csv`\n');
fprintf(fid, '- `b1_peak_evidence_audit_summary.csv`\n');
fprintf(fid, '- `b1_peak_evidence_constraint_release.csv`\n');
fprintf(fid, '- `b1_peak_evidence_robustness.csv`\n');
fprintf(fid, '- `b1_peak_evidence_candidate_competition.csv`\n');
fprintf(fid, '- `b1_peak_evidence_suspicious_points.csv`\n\n');
fprintf(fid, '## Figures\n\n');
fprintf(fid, '- `%s`\n', local_file_name(summary_png));
fprintf(fid, '- `%s`\n\n', local_file_name(panel_png));
fprintf(fid, '## Summary\n\n');
for i = 1:height(summary)
    fprintf(fid, '- `%s`: total %d, supported %d, assisted %d, suspicious %d\n', ...
        summary.session_key{i}, summary.n_points(i), ...
        summary.n_data_supported(i), summary.n_tracking_assisted(i), ...
        summary.n_suspicious(i));
end
end


function name = local_file_name(path)
[~, n, e] = fileparts(char(path));
name = [n e];
end


function local_update_manifest(project_root, summary_dir)
manifest_path = fullfile(project_root, 'paper_results', '00-by_date', ...
    '_manifest.csv');
if isfile(manifest_path)
    manifest = readtable(manifest_path, 'TextType', 'string');
else
    manifest = table('Size', [0 6], ...
        'VariableTypes', {'string', 'string', 'string', 'string', ...
        'string', 'string'}, ...
        'VariableNames', {'date_group', 'date_source', 'result_name', ...
        'last_write_time', 'link_path', 'target_path'});
end
[~, result_name] = fileparts(summary_dir);
row = table("2026-05-11", "name_suffix", string(result_name), ...
    string(datestr(now, 'yyyy-mm-dd HH:MM:SS')), ...
    string(fullfile(project_root, 'paper_results', 'by_date', ...
    '2026-05-11', result_name)), string(summary_dir), ...
    'VariableNames', manifest.Properties.VariableNames);
if any(strcmp(manifest.result_name, result_name))
    manifest(strcmp(manifest.result_name, result_name), :) = [];
end
manifest = [manifest; row];
writetable(manifest, manifest_path);
end


function local_update_results_index(project_root, summary_dir, summary)
index_path = fullfile(project_root, 'case_studies', 'bisb2026', ...
    'RESULTS_INDEX.md');
if isfile(index_path)
    text = fileread(index_path);
else
    text = "# BiSb Results Index" + newline;
end
[~, result_name] = fileparts(summary_dir);
total = summary(strcmp(summary.session_key, 'all'), :);
if isempty(total)
    total_line = '- Total audited points: `0`';
else
    total_line = sprintf(['- Total audited points: `%d`; supported `%d`; ' ...
        'assisted `%d`; suspicious `%d`'], total.n_points(1), ...
        total.n_data_supported(1), total.n_tracking_assisted(1), ...
        total.n_suspicious(1));
end
block = sprintf(['<!-- B1_PEAK_EVIDENCE_AUDIT_260511_START -->\n' ...
    '## 2026-05-11 B1 Lorentz peak evidence audit\n\n' ...
    '- Summary: `%s`\n' ...
    '- Source tag: `260510_lorentz_tracking_v15_upper_rapidrise_plateau`\n' ...
    '- Energy window: `300-1800 meV`\n' ...
    '- Peak model: `lorentz`\n' ...
    '- Physical fit: no physical fit was run in this evidence audit.\n' ...
    '%s\n\n' ...
    '<!-- B1_PEAK_EVIDENCE_AUDIT_260511_END -->\n'], ...
    fullfile('paper_results', result_name), total_line);
text = local_replace_block(text, 'B1_PEAK_EVIDENCE_AUDIT_260511', block);
fid = fopen(index_path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);
end


function text = local_replace_block(text, key, block)
start_tag = ['<!-- ' key '_START -->'];
end_tag = ['<!-- ' key '_END -->'];
start_idx = strfind(text, start_tag);
end_idx = strfind(text, end_tag);
if ~isempty(start_idx) && ~isempty(end_idx) && end_idx(1) > start_idx(1)
    end_pos = end_idx(1) + numel(end_tag) - 1;
    text = [text(1:start_idx(1)-1), block, text(end_pos+1:end)];
else
    if isempty(text) || text(end) ~= newline
        text = [text newline];
    end
    text = [text newline block];
end
text = char(text);
end
