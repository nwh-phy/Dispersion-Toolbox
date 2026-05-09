function output = run_590_gui_history_area_analysis(sessionRequest, options)
%RUN_590_GUI_HISTORY_AREA_ANALYSIS Reproduce GUI-history area analyses.
%   output = run_590_gui_history_area_analysis()
%   output = run_590_gui_history_area_analysis("remaining")
%   output = run_590_gui_history_area_analysis("all")
%
%   Mirrors the current interactive_qe_browser auto-fit path and forces
%   preprocessing normalization to Area over the displayed energy window.

arguments
    sessionRequest {mustBeTextScalar} = "590_PL2_10w"
    options.qRangeOverride_Ainv (1,2) double = [-0.15 0.15]
    options.outputTagSuffix {mustBeTextScalar} = ""
    options.peakModelOverride {mustBeTextScalar} = ""
end

script_path = mfilename('fullpath');
project_root = bisb_find_project_root(fileparts(script_path));
run(fullfile(project_root, 'startup.m'));

sessions = local_requested_sessions(project_root, sessionRequest);
sessions = local_apply_run_options_to_session(sessions, options);
if numel(sessions) > 1
    outputs = cell(1, numel(sessions));
    for i = 1:numel(sessions)
        outputs{i} = local_run_one_session(project_root, sessions(i), options);
    end
    combined_dir = local_write_combined_outputs(project_root, outputs, sessionRequest, options);
    output = struct();
    output.sessions = outputs;
    output.output_dir = combined_dir;
    fprintf('Combined output directory: %s\n', combined_dir);
else
    output = local_run_one_session(project_root, sessions(1), options);
end
end


function output = local_run_one_session(project_root, session, options)
data_dir = session.data_dir;
data_path = fullfile(data_dir, 'eq3D.mat');
history_entry = local_load_history_entry(session);
history_path = history_entry.source_path;
out_dir = fullfile(project_root, 'paper_results', session.output_tag);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

snap = history_entry.snapshot;
snap = local_force_area_snapshot(snap);
if all(isfinite(options.qRangeOverride_Ainv))
    snap.qStart = options.qRangeOverride_Ainv(1);
    snap.qEnd = options.qRangeOverride_Ainv(2);
end
if strlength(string(options.peakModelOverride)) > 0
    snap.peakModel = char(string(options.peakModelOverride));
end
snap.dqOverride = session.dq_Ainv;

dataset = load_qe_dataset(data_path, snap.dqOverride);
qe = dataset.qe;

pp_opts = local_preprocess_opts_from_snapshot(snap);
[qe_pp, bg_diag] = qe_preprocess(qe, pp_opts);
raw_opts = pp_opts;
raw_opts.do_normalize = false;
qe_raw = qe_preprocess(qe, raw_opts);

branch_specs = local_branch_specs_from_snapshot(snap);
[energy_mask, energy_axis] = local_energy_mask(qe_pp, snap);

auto_opts = local_auto_fit_opts_from_snapshot(snap, branch_specs, ...
    energy_mask, energy_axis);
fit_res = qe_auto_fit(qe_pp, qe_raw, auto_opts);

branch_filter_opts = struct();
branch_filter_opts.q_skip_Ainv = 0.005;
branch_filter_opts.min_R2 = 0.3;
branch_filter_opts.max_gamma_ratio = 2.0;
branch_filter_opts.score_column = 12;
assignment = qe_assign_peak_branches_by_windows( ...
    fit_res.all_peaks, branch_specs, branch_filter_opts);
branches = assignment.branches;
[branches, refinement_log] = local_apply_session_refinement( ...
    qe_pp, qe_raw, branches, snap, session.refinement_profile, auto_opts);
assignment = local_update_assignment_summary_from_branches(assignment, branches);

model_results = local_fit_model_suite(branches);
branch_summary = local_branch_summary(branches, model_results, snap.dqOverride);
model_summary = local_model_summary_table(model_results);
single_summary = local_single_spectrum_summary(qe_pp, fit_res, snap);

local_write_branch_tables(branches, out_dir);
writetable(assignment.summary, fullfile(out_dir, 'branch_assignment_summary.csv'));
if ~isempty(assignment.rejected)
    writetable(assignment.rejected, fullfile(out_dir, 'rejected_peaks.csv'));
end
writetable(refinement_log, fullfile(out_dir, 'branch_refinement_log.csv'));
writetable(branch_summary, fullfile(out_dir, 'branch_summary.csv'));
writetable(model_summary, fullfile(out_dir, 'dispersion_model_summary.csv'));
writetable(single_summary, fullfile(out_dir, 'single_spectrum_fit_summary.csv'));

figure_paths = local_export_figures(qe_pp, branches, ...
    model_results, fit_res, snap, out_dir);
report_path = local_write_report(out_dir, dataset, history_entry, snap, ...
    pp_opts, fit_res, branch_summary, model_summary, ...
    single_summary, figure_paths, refinement_log, session.refinement_profile);

output = struct();
output.dataset = dataset;
output.history_path = history_path;
output.session = session;
output.history_entry = history_entry;
output.snap = snap;
output.preprocess_opts = pp_opts;
output.qe_pp = qe_pp;
output.bg_diag = bg_diag;
output.fit_res = fit_res;
output.assignment = assignment;
output.branches = branches;
output.refinement_log = refinement_log;
output.refinement_profile = session.refinement_profile;
output.model_results = model_results;
output.branch_summary = branch_summary;
output.model_summary = model_summary;
output.single_summary = single_summary;
output.figure_paths = figure_paths;
output.report_path = report_path;
output.output_dir = out_dir;

save(fullfile(out_dir, 'analysis_results.mat'), 'output', ...
    'dataset', 'snap', 'pp_opts', 'fit_res', 'assignment', ...
    'branches', 'model_results', 'branch_summary', 'model_summary', ...
    'single_summary', 'figure_paths', 'report_path', ...
    'refinement_log', '-v7.3');

fprintf('Report: %s\n', report_path);
fprintf('Output directory: %s\n', out_dir);
end


function sessions = local_requested_sessions(project_root, sessionRequest)
all_sessions = local_session_configs(project_root);
request = lower(strtrim(char(sessionRequest)));
switch request
    case {'590', '590_pl2_10w', 'primary'}
        keep = strcmp({all_sessions.key}, '590_PL2_10w');
    case {'remaining', 'rest'}
        keep = ~strcmp({all_sessions.key}, '590_PL2_10w');
    case {'all', '*'}
        keep = true(size(all_sessions));
    otherwise
        keys = lower(string({all_sessions.key}));
        keep = keys == string(request);
end
sessions = all_sessions(keep);
if isempty(sessions)
    error('run_590_gui_history_area_analysis:UnknownSession', ...
        'Unknown session request "%s".', char(sessionRequest));
end
end


function sessions = local_apply_run_options_to_session(sessions, options)
if strlength(string(options.outputTagSuffix)) == 0
    return
end
for i = 1:numel(sessions)
    session = sessions(i);
    session.output_tag = [session.output_tag char(options.outputTagSuffix)];
    sessions(i) = session;
end
end


function sessions = local_session_configs(project_root)
base_dir = fullfile(project_root, '20260120 BiSb');
template_history = fullfile(base_dir, '590 PL2 10w 0.004 10sx300', ...
    'op_history_260506.mat');

sessions = repmat(struct('key', '', 'display_name', '', 'data_dir', '', ...
    'dq_Ainv', NaN, 'history_path', '', 'template_history_path', ...
    template_history, 'output_tag', '', ...
    'refinement_profile', local_empty_refinement_profile()), 1, 3);

sessions(1).key = '590_PL2_10w';
sessions(1).display_name = '590 PL2 10w';
sessions(1).data_dir = fullfile(base_dir, '590 PL2 10w 0.004 10sx300');
sessions(1).dq_Ainv = 0.005;
sessions(1).history_path = template_history;
sessions(1).output_tag = '590_gui_history_area_260506';
sessions(1).refinement_profile = local_empty_refinement_profile();

sessions(2).key = 'n0_PL2_10w_repeat';
sessions(2).display_name = 'n0 PL2 10w repeat';
sessions(2).data_dir = fullfile(base_dir, 'n0 pl2 10w 0.004 10s x300');
sessions(2).dq_Ainv = 0.005;
sessions(2).history_path = '';
sessions(2).output_tag = 'n0_PL2_10w_gui_history_area_260506';
sessions(2).refinement_profile = local_empty_refinement_profile();

sessions(3).key = 'no_PL2_20w_2film';
sessions(3).display_name = 'no PL2 20w 2film';
sessions(3).data_dir = fullfile(base_dir, 'no pl2 20w 0.004 10sx300 2film');
sessions(3).dq_Ainv = 0.0025;
% The local 20w history is from an older GUI state and enables background
% subtraction. Use the current 590 saved history as the cross-session
% template, while preserving this dataset's dq calibration.
sessions(3).history_path = '';
sessions(3).output_tag = 'no_PL2_20w_2film_gui_history_area_260506_highq_refined';
sessions(3).refinement_profile = local_20w_highq_refinement_profile();
end


function profile = local_empty_refinement_profile()
profile = struct();
profile.enabled = false;
profile.label = 'none';
profile.branch_index = NaN;
profile.q_min_Ainv = Inf;
profile.ci_half_max_meV = Inf;
profile.gamma_ratio_max = Inf;
profile.refit_window_meV = [NaN NaN];
profile.max_peaks = 1;
profile.prominence = NaN;
profile.smooth_width = NaN;
profile.min_R2 = -Inf;
profile.max_refit_ci_half_meV = Inf;
profile.max_refit_gamma_ratio = Inf;
profile.max_energy_shift_meV = Inf;
profile.window_edge_margin_meV = 0;
end


function profile = local_20w_highq_refinement_profile()
profile = local_empty_refinement_profile();
profile.enabled = true;
profile.label = '20w_B1_highq_refit';
profile.branch_index = 1;
profile.q_min_Ainv = 0.10;
profile.ci_half_max_meV = 200;
profile.gamma_ratio_max = 1.25;
profile.refit_window_meV = [1000 1700];
profile.max_peaks = 1;
profile.prominence = 0.08;
profile.smooth_width = 17;
profile.min_R2 = 0.30;
profile.max_refit_ci_half_meV = 200;
profile.max_refit_gamma_ratio = 1.25;
profile.max_energy_shift_meV = 250;
profile.window_edge_margin_meV = 50;
end


function history_entry = local_load_history_entry(session)
if strlength(string(session.history_path)) > 0 && isfile(session.history_path)
    history_path = session.history_path;
    source_note = 'session_history';
else
    history_path = session.template_history_path;
    source_note = 'template_590_history';
end
loaded = load(history_path, 'opHistory');
history_entry = loaded.opHistory{end};
history_entry.source_path = history_path;
history_entry.source_note = source_note;
if strcmp(source_note, 'template_590_history')
    history_entry.label = sprintf('%s used as template for %s', ...
        char(string(history_entry.label)), session.key);
end
end


function combined_dir = local_write_combined_outputs(project_root, outputs, sessionRequest, options)
request_tag = regexprep(lower(char(string(sessionRequest))), '[^a-z0-9]+', '_');
suffix = char(options.outputTagSuffix);
combined_dir = fullfile(project_root, 'paper_results', ...
    sprintf('gui_history_area_260506_%s%s', request_tag, suffix));
if ~exist(combined_dir, 'dir')
    mkdir(combined_dir);
end

branch_all = table();
model_all = table();
for i = 1:numel(outputs)
    out = outputs{i};
    branch_tbl = out.branch_summary;
    branch_tbl.session = repmat({out.session.key}, height(branch_tbl), 1);
    branch_tbl = movevars(branch_tbl, 'session', 'Before', 1);
    branch_all = [branch_all; branch_tbl]; %#ok<AGROW>

    model_tbl = out.model_summary;
    model_tbl.session = repmat({out.session.key}, height(model_tbl), 1);
    model_tbl = movevars(model_tbl, 'session', 'Before', 1);
    model_all = [model_all; model_tbl]; %#ok<AGROW>
end

writetable(branch_all, fullfile(combined_dir, 'combined_branch_summary.csv'));
writetable(model_all, fullfile(combined_dir, 'combined_model_summary.csv'));

report_path = fullfile(combined_dir, 'combined_report.md');
fid = fopen(report_path, 'w', 'n', 'UTF-8');
if fid < 0
    error('run_590_gui_history_area_analysis:CannotWriteCombinedReport', ...
        'Cannot write %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# Combined area-normalized GUI-history analysis\n\n');
fprintf(fid, '- Request: `%s`\n', char(string(sessionRequest)));
fprintf(fid, '- Normalization: `Area`, using each session display-energy window.\n');
fprintf(fid, '- Heatmap evidence: physical q-E map only, matching the GUI top-left panel.\n\n');

fprintf(fid, '## Session outputs\n\n');
for i = 1:numel(outputs)
    out = outputs{i};
    fprintf(fid, '- `%s`: `%s`\n', out.session.key, out.output_dir);
end
fprintf(fid, '\n');

fprintf(fid, '## Cross-session branch snapshot\n\n');
fprintf(fid, '| Session | Branch | N | |q| coverage | Energy range (meV) | Edge-inner (meV) | Median CI half-width (meV) | Best BIC model |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|---|\n');
for i = 1:height(branch_all)
    fprintf(fid, '| %s | B%d | %d | %.3f-%.3f | %.0f-%.0f | %.1f | %.2f | %s |\n', ...
        branch_all.session{i}, branch_all.branch(i), branch_all.n_points(i), ...
        branch_all.q_abs_min_Ainv(i), branch_all.q_abs_max_Ainv(i), ...
        branch_all.energy_min_meV(i), branch_all.energy_max_meV(i), ...
        branch_all.energy_edge_minus_inner_meV(i), ...
        branch_all.energy_ci_half_median_meV(i), ...
        char(string(branch_all.best_model_by_BIC{i})));
end
fprintf(fid, '\n');

fprintf(fid, '## Files\n\n');
fprintf(fid, '- `combined_branch_summary.csv`\n');
fprintf(fid, '- `combined_model_summary.csv`\n');
end


function snap = local_force_area_snapshot(snap)
snap = local_apply_snapshot_defaults(snap);
snap.normMethod = 'Area';
snap.areaNorm = true;
if ~isfield(snap, 'energyMin') || ~isfinite(snap.energyMin)
    snap.energyMin = 200;
end
if ~isfield(snap, 'energyMax') || ~isfinite(snap.energyMax)
    snap.energyMax = 3876;
end
snap.areaNormMin = snap.energyMin;
snap.areaNormMax = snap.energyMax;
if ~isfield(snap, 'exportMode')
    snap.exportMode = 'PPT PNG';
end
end


function snap = local_apply_snapshot_defaults(snap)
defaults = struct( ...
    'qStart', -0.15, ...
    'qEnd', 0.15, ...
    'qStep', 0.005, ...
    'dqOverride', 0.005, ...
    'energyMin', 200, ...
    'energyMax', 3876, ...
    'refQMin', 0.010, ...
    'refQMax', 0.020, ...
    'smoothMode', 'gaussian', ...
    'smoothWidth', 3, ...
    'sgOrder', 3, ...
    'sgFrameLen', 15, ...
    'areaNorm', true, ...
    'bgSub', false, ...
    'bgMethod', 'Auto', ...
    'bgWin1Lo', 50, ...
    'bgWin1Hi', 300, ...
    'bgDual', false, ...
    'bgWin2Lo', 3000, ...
    'bgWin2Hi', 3500, ...
    'bgIter', false, ...
    'denoise', true, ...
    'denoiseMethod', 'Wiener2D', ...
    'denoiseSigma', 0, ...
    'deconv', false, ...
    'deconvIter', 5, ...
    'prominence', 0.10, ...
    'autoFitSmoothWidth', 25, ...
    'maxPeaks', 3, ...
    'peakModel', 'fano', ...
    'bootstrapCiSamples', 0, ...
    'maxShift', 80, ...
    'branch1Min', 500, ...
    'branch1Max', 2100, ...
    'branch2Min', 1800, ...
    'branch2Max', 2500, ...
    'branch3Min', 2800, ...
    'branch3Max', 3800, ...
    'selectedQIndex', 1, ...
    'selectedQ_Ainv', 0);
fields = fieldnames(defaults);
for i = 1:numel(fields)
    name = fields{i};
    if ~isfield(snap, name) || isempty(snap.(name))
        snap.(name) = defaults.(name);
    end
end
end


function opts = local_preprocess_opts_from_snapshot(snap)
opts = struct();
opts.do_despike = false;
opts.do_normalize = logical(local_snap_value(snap, 'areaNorm', true));
opts.norm_method = 'Area';
opts.norm_min = local_snap_value(snap, 'energyMin', 200);
opts.norm_max = local_snap_value(snap, 'energyMax', 3876);
opts.do_denoise = logical(local_snap_value(snap, 'denoise', true));
opts.denoise_method = char(local_snap_value(snap, 'denoiseMethod', 'Wiener2D'));
opts.denoise_sigma = local_snap_value(snap, 'denoiseSigma', 0);
opts.sg_order = local_snap_value(snap, 'sgOrder', 3);
opts.sg_framelen = local_snap_value(snap, 'sgFrameLen', 15);
opts.do_bg_sub = logical(local_snap_value(snap, 'bgSub', false));
opts.bg_method = char(local_snap_value(snap, 'bgMethod', 'Auto'));
opts.bg_win_lo = [local_snap_value(snap, 'bgWin1Lo', 50), ...
    local_snap_value(snap, 'bgWin1Hi', 300)];
if logical(local_snap_value(snap, 'bgDual', false))
    opts.bg_win_hi = [local_snap_value(snap, 'bgWin2Lo', 3000), ...
        local_snap_value(snap, 'bgWin2Hi', 3500)];
else
    opts.bg_win_hi = [];
end
opts.bg_iterative = logical(local_snap_value(snap, 'bgIter', false));
opts.do_deconv = logical(local_snap_value(snap, 'deconv', false));
opts.deconv_iter = local_snap_value(snap, 'deconvIter', 5);
end


function specs = local_branch_specs_from_snapshot(snap)
windows = [
    local_snap_value(snap, 'branch1Min', 500), local_snap_value(snap, 'branch1Max', 2100)
    local_snap_value(snap, 'branch2Min', 1800), local_snap_value(snap, 'branch2Max', 2500)
    local_snap_value(snap, 'branch3Min', 2800), local_snap_value(snap, 'branch3Max', 3800)
    ];
specs = repmat(struct('name', '', ...
    'energy_window_meV', [0 0], 'enabled', true), 3, 1);
for b = 1:3
    specs(b).name = sprintf('Branch %d', b);
    specs(b).energy_window_meV = sort(double(windows(b, :)));
    specs(b).enabled = true;
end
end


function [mask, energy_axis] = local_energy_mask(qe, snap)
energy_min = min(snap.energyMin, snap.energyMax);
energy_max = max(snap.energyMin, snap.energyMax);
energy_full = double(qe.energy_meV(:));
mask = energy_full >= energy_min & energy_full <= energy_max;
if ~any(mask)
    [~, nearest] = min(abs(energy_full - mean([energy_min, energy_max])));
    mask(nearest) = true;
end
energy_axis = energy_full(mask);
end


function opts = local_auto_fit_opts_from_snapshot(snap, branch_specs, ...
    energy_mask, energy_axis)
opts = struct();
opts.E_min = max(local_snap_value(snap, 'energyMin', 200), 50);
opts.E_max = local_snap_value(snap, 'energyMax', 3876);
q_start = local_snap_value(snap, 'qStart', -0.15);
q_end = local_snap_value(snap, 'qEnd', 0.15);
opts.q_start = min(q_start, q_end);
opts.q_end = max(q_start, q_end);
opts.prominence = local_snap_value(snap, 'prominence', 0.10);
opts.smooth_width = local_snap_value(snap, 'autoFitSmoothWidth', 25);
opts.max_peaks = local_snap_value(snap, 'maxPeaks', 3);
opts.peak_model = char(local_snap_value(snap, 'peakModel', 'fano'));
opts.bootstrap_ci_samples = local_snap_value(snap, 'bootstrapCiSamples', 0);
opts.branch_specs = branch_specs;
opts.pre_subtracted = logical(local_snap_value(snap, 'bgSub', false));
opts.guesses = [];
opts.seed_idx = local_snap_value(snap, 'selectedQIndex', 1);
opts.max_shift = local_snap_value(snap, 'maxShift', 80);
opts.energy_mask = energy_mask;
opts.energy_axis = energy_axis;
opts.R2_threshold = 0.3;
opts.window_seed_branch_indices = 2;
opts.verbose = true;
opts.progress_fn = [];
end


function value = local_snap_value(snap, name, default_value)
if isfield(snap, name) && ~isempty(snap.(name))
    value = snap.(name);
else
    value = default_value;
end
end


function [branches, refinement_log] = local_apply_session_refinement( ...
    qe_pp, qe_raw, branches, snap, profile, auto_opts)
refinement_log = local_empty_refinement_log();
if ~isstruct(profile) || ~isfield(profile, 'enabled') || ~profile.enabled
    return
end

branch_idx = profile.branch_index;
if ~isfinite(branch_idx) || branch_idx < 1 || branch_idx > numel(branches)
    return
end

branch = branches{branch_idx};
if isempty(branch)
    return
end

kept = zeros(0, size(branch, 2));
for i = 1:size(branch, 1)
    old_row = branch(i, :);
    [needs_refine, reason] = local_should_refine_row(old_row, profile);
    if ~needs_refine
        kept = [kept; old_row]; %#ok<AGROW>
        continue
    end

    [new_row, status, detail] = local_refit_branch_row( ...
        qe_pp, qe_raw, old_row, snap, profile, auto_opts);
    if strcmp(status, 'accepted')
        kept = [kept; new_row]; %#ok<AGROW>
        action = 'replaced';
    else
        action = 'rejected';
    end

    refinement_log = [refinement_log; local_refinement_log_row( ...
        profile, old_row, new_row, reason, action, status, detail)]; %#ok<AGROW>
end

branches{branch_idx} = sortrows(kept, 1);
end


function [needs_refine, reason] = local_should_refine_row(row, profile)
needs_refine = false;
reasons = {};
if abs(row(1)) < profile.q_min_Ainv
    reason = '';
    return
end

ci_half = local_row_ci_half(row);
gamma_ratio = local_row_gamma_ratio(row);
if isfinite(ci_half) && ci_half > profile.ci_half_max_meV
    reasons{end+1} = sprintf('CI half %.1f > %.1f meV', ... %#ok<AGROW>
        ci_half, profile.ci_half_max_meV);
end
if isfinite(gamma_ratio) && gamma_ratio > profile.gamma_ratio_max
    reasons{end+1} = sprintf('Gamma/E %.2f > %.2f', ... %#ok<AGROW>
        gamma_ratio, profile.gamma_ratio_max);
end

needs_refine = ~isempty(reasons);
reason = strjoin(reasons, '; ');
end


function [new_row, status, detail] = local_refit_branch_row( ...
    qe_pp, qe_raw, old_row, snap, profile, auto_opts)
new_row = old_row;
detail = '';

q_axis = double(qe_pp.q_Ainv(:));
[q_delta, q_idx] = min(abs(q_axis - old_row(1)));
dq_tol = max(local_snap_value(snap, 'dqOverride', 0.005), 1e-6);
if ~isfinite(q_delta) || q_delta > dq_tol
    status = 'missing_q_channel';
    detail = sprintf('nearest q delta %.4g exceeds tolerance %.4g', q_delta, dq_tol);
    return
end

win = sort(double(profile.refit_window_meV(:)).');
if numel(win) ~= 2 || any(~isfinite(win)) || win(1) >= win(2)
    status = 'invalid_refit_window';
    return
end

energy_axis = auto_opts.energy_axis(:);
energy_mask = auto_opts.energy_mask(:);
spectrum = double(qe_pp.intensity(energy_mask, q_idx));
raw_spectrum = double(qe_raw.intensity(energy_mask, q_idx));
if all(~isfinite(spectrum)) || all(spectrum == 0)
    status = 'empty_spectrum';
    return
end

guess = min(max(old_row(2), win(1)), win(2));
prominence = profile.prominence;
if ~isfinite(prominence)
    prominence = auto_opts.prominence;
end
smooth_width = profile.smooth_width;
if ~isfinite(smooth_width)
    smooth_width = auto_opts.smooth_width;
end

try
    result = fit_loss_function(energy_axis, spectrum, ...
        'E_min', win(1), 'E_max', win(2), ...
        'min_prominence', prominence, ...
        'smooth_width', smooth_width, ...
        'max_peaks', profile.max_peaks, ...
        'initial_guesses', guess, ...
        'peak_model', auto_opts.peak_model, ...
        'pre_subtracted', auto_opts.pre_subtracted, ...
        'bootstrap_ci_samples', auto_opts.bootstrap_ci_samples);
catch ME
    status = 'fit_failed';
    detail = ME.message;
    return
end

if result.n_peaks < 1
    status = 'no_peak';
    return
end

[peak_energy, peak_ci] = local_branch_peak_energy(result);
[~, peak_idx] = min(abs(peak_energy - guess));
if isempty(peak_idx) || ~isfinite(peak_energy(peak_idx))
    status = 'no_finite_peak';
    return
end

new_row = old_row;
new_row(1) = q_axis(q_idx);
new_row(2) = peak_energy(peak_idx);
new_row(3) = result.gamma(peak_idx);
new_row(4) = result.R_squared;
new_row(5) = result.amplitude(peak_idx);
new_row(6:7) = peak_ci(peak_idx, :);
new_row(8:9) = local_result_ci_row(result.gamma_ci, peak_idx, new_row(3));
new_row(10:11) = local_result_ci_row(result.amplitude_ci, peak_idx, new_row(5));
new_row(12) = measure_peak_height( ...
    energy_axis, raw_spectrum, new_row(2), new_row(3));

[passed, detail] = local_refit_quality(new_row, old_row, win, profile);
if passed
    status = 'accepted';
else
    status = 'rejected_low_confidence_refit';
end
end


function [passed, detail] = local_refit_quality(row, old_row, win, profile)
reasons = {};
ci_half = local_row_ci_half(row);
gamma_ratio = local_row_gamma_ratio(row);
if row(4) < profile.min_R2
    reasons{end+1} = sprintf('R2 %.3f < %.3f', row(4), profile.min_R2);
end
if ~isfinite(ci_half) || ci_half > profile.max_refit_ci_half_meV
    reasons{end+1} = sprintf('CI half %.1f > %.1f meV', ... %#ok<AGROW>
        ci_half, profile.max_refit_ci_half_meV);
end
if ~isfinite(gamma_ratio) || gamma_ratio > profile.max_refit_gamma_ratio
    reasons{end+1} = sprintf('Gamma/E %.2f > %.2f', ... %#ok<AGROW>
        gamma_ratio, profile.max_refit_gamma_ratio);
end
energy_shift = abs(row(2) - old_row(2));
if isfinite(profile.max_energy_shift_meV) && ...
        energy_shift > profile.max_energy_shift_meV
    reasons{end+1} = sprintf('energy shift %.1f > %.1f meV', ... %#ok<AGROW>
        energy_shift, profile.max_energy_shift_meV);
end
edge_margin = profile.window_edge_margin_meV;
if isfinite(edge_margin) && edge_margin > 0 && ...
        (row(2) <= win(1) + edge_margin || row(2) >= win(2) - edge_margin)
    reasons{end+1} = sprintf('energy %.1f within %.1f meV of refit window edge', ... %#ok<AGROW>
        row(2), edge_margin);
end
passed = isempty(reasons);
if passed
    detail = 'ok';
else
    detail = strjoin(reasons, '; ');
end
end


function ci = local_result_ci_row(ci_matrix, idx, center)
ci = [NaN NaN];
if idx <= size(ci_matrix, 1) && size(ci_matrix, 2) >= 2
    ci = ci_matrix(idx, 1:2);
end
if ~all(isfinite(ci)) || ci(1) >= center || ci(2) <= center
    half_width = max(1, 0.001 * abs(center));
    ci = [center - half_width, center + half_width];
end
end


function [peak_energy, peak_ci] = local_branch_peak_energy(result)
peak_energy = result.omega_p(:);
peak_ci = local_result_ci(result, 'omega_p_ci', numel(peak_energy));
if isfield(result, 'apex_energy_meV')
    apex = result.apex_energy_meV(:);
    use_apex = isfinite(apex);
    peak_energy(use_apex) = apex(use_apex);
    if isfield(result, 'apex_energy_ci')
        apex_ci = result.apex_energy_ci;
        valid_ci_rows = use_apex & (1:numel(peak_energy))' <= size(apex_ci, 1);
        peak_ci(valid_ci_rows, :) = apex_ci(valid_ci_rows, :);
    end
end
peak_ci = local_fill_energy_ci(peak_energy, peak_ci);
end


function ci = local_result_ci(result, field_name, n_rows)
ci = NaN(n_rows, 2);
if isfield(result, field_name)
    field_ci = result.(field_name);
    n = min(n_rows, size(field_ci, 1));
    if size(field_ci, 2) >= 2
        ci(1:n, :) = field_ci(1:n, 1:2);
    end
end
end


function ci = local_fill_energy_ci(energy, ci)
for i = 1:numel(energy)
    center = energy(i);
    if ~isfinite(center)
        continue
    end
    half_width = max(1, 0.001 * abs(center));
    if ~all(isfinite(ci(i, :))) || ci(i, 1) >= center || ci(i, 2) <= center
        ci(i, :) = [center - half_width, center + half_width];
    else
        ci(i, 1) = min(ci(i, 1), center - half_width);
        ci(i, 2) = max(ci(i, 2), center + half_width);
    end
end
end


function half_width = local_row_ci_half(row)
half_width = NaN;
if numel(row) >= 7 && all(isfinite(row(6:7)))
    half_width = max(abs(row(6:7) - row(2)));
end
if ~isfinite(half_width)
    half_width = max(1, 0.001 * abs(row(2)));
end
end


function gamma_ratio = local_row_gamma_ratio(row)
gamma_ratio = row(3) ./ max(row(2), eps);
end


function log_row = local_refinement_log_row( ...
    profile, old_row, new_row, reason, action, status, detail)
log_row = table({profile.label}, profile.branch_index, old_row(1), ...
    abs(old_row(1)), old_row(2), new_row(2), ...
    local_row_ci_half(old_row), local_row_ci_half(new_row), ...
    local_row_gamma_ratio(old_row), local_row_gamma_ratio(new_row), ...
    old_row(4), new_row(4), {reason}, {action}, {status}, {detail}, ...
    'VariableNames', {'profile', 'branch', 'q_Ainv', 'q_abs_Ainv', ...
    'old_energy_meV', 'new_energy_meV', 'old_CI_half_meV', ...
    'new_CI_half_meV', 'old_gamma_over_E', 'new_gamma_over_E', ...
    'old_R2', 'new_R2', 'trigger_reason', 'action', 'status', 'detail'});
end


function refinement_log = local_empty_refinement_log()
refinement_log = table(cell(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    cell(0,1), cell(0,1), cell(0,1), cell(0,1), ...
    'VariableNames', {'profile', 'branch', 'q_Ainv', 'q_abs_Ainv', ...
    'old_energy_meV', 'new_energy_meV', 'old_CI_half_meV', ...
    'new_CI_half_meV', 'old_gamma_over_E', 'new_gamma_over_E', ...
    'old_R2', 'new_R2', 'trigger_reason', 'action', 'status', 'detail'});
end


function assignment = local_update_assignment_summary_from_branches(assignment, branches)
if ~isfield(assignment, 'summary') || isempty(assignment.summary)
    return
end
n = min(height(assignment.summary), numel(branches));
for i = 1:n
    assignment.summary.assigned_n(i) = size(branches{i}, 1);
end
end


function model_results = local_fit_model_suite(branches)
models = {'quasi2d_plasmon', 'optical_constant', 'optical_quadratic'};
model_results = cell(numel(branches), 1);
for b = 1:numel(branches)
    br = branches{b};
    entries = repmat(struct('model', '', 'success', false, 'fit', struct(), ...
        'error', '', 'aic', NaN, 'bic', NaN), numel(models), 1);
    for mi = 1:numel(models)
        entries(mi).model = models{mi};
        if size(br, 1) < 5
            entries(mi).error = 'insufficient_points';
            continue
        end
        weights = br(:, 4);
        weights(~isfinite(weights) | weights <= 0) = 1;
        try
            fit = fit_dispersion_generic(br(:, 1), br(:, 2), ...
                'model', models{mi}, 'confidence', weights);
            [aic, bic] = local_information_criteria(fit);
            entries(mi).success = true;
            entries(mi).fit = fit;
            entries(mi).aic = aic;
            entries(mi).bic = bic;
        catch ME
            entries(mi).error = ME.message;
        end
    end
    model_results{b} = entries;
end
end


function [aic, bic] = local_information_criteria(fit)
res = double(fit.residuals_meV(:));
res = res(isfinite(res));
n = numel(res);
k = numel(fit.params);
if n == 0
    aic = NaN;
    bic = NaN;
    return
end
rss = max(sum(res .^ 2), eps);
aic = n * log(rss / n) + 2 * k;
bic = n * log(rss / n) + k * log(n);
end


function tbl = local_branch_summary(branches, model_results, dq_Ainv)
rows = cell(numel(branches), 22);
for b = 1:numel(branches)
    br = branches{b};
    stats = local_branch_stats(br, dq_Ainv);
    [best_model, best_bic] = local_best_model(model_results{b});
    rows(b, :) = { ...
        b, size(br, 1), stats.q_min, stats.q_max, ...
        stats.q_abs_min, stats.q_abs_max, ...
        stats.energy_min, stats.energy_max, stats.energy_center_mean, ...
        stats.energy_inner_mean, stats.energy_edge_mean, ...
        stats.energy_edge_minus_center, stats.energy_edge_minus_inner, ...
        stats.gamma_median, stats.gamma_ratio_median, ...
        stats.ci_half_median, stats.ci_half_max, ...
        stats.sym_pair_n, stats.sym_mean_abs_delta, stats.sym_max_abs_delta, ...
        best_model, best_bic};
end
tbl = cell2table(rows, 'VariableNames', { ...
    'branch', 'n_points', 'q_min_Ainv', 'q_max_Ainv', ...
    'q_abs_min_Ainv', 'q_abs_max_Ainv', ...
    'energy_min_meV', 'energy_max_meV', 'energy_center_mean_meV', ...
    'energy_inner_mean_meV', 'energy_edge_mean_meV', ...
    'energy_edge_minus_center_meV', 'energy_edge_minus_inner_meV', ...
    'gamma_median_meV', 'gamma_over_E_median', ...
    'energy_ci_half_median_meV', 'energy_ci_half_max_meV', ...
    'symmetry_pair_n', 'symmetry_mean_abs_delta_meV', ...
    'symmetry_max_abs_delta_meV', 'best_model_by_BIC', 'best_BIC'});
end


function stats = local_branch_stats(br, dq_Ainv)
stats = struct();
if isempty(br)
    fields = {'q_min','q_max','q_abs_min','q_abs_max','energy_min', ...
        'energy_max','energy_center_mean','energy_inner_mean', ...
        'energy_edge_mean','energy_edge_minus_center', ...
        'energy_edge_minus_inner','gamma_median', ...
        'gamma_ratio_median','ci_half_median','ci_half_max', ...
        'sym_pair_n','sym_mean_abs_delta','sym_max_abs_delta'};
    for i = 1:numel(fields)
        stats.(fields{i}) = NaN;
    end
    stats.sym_pair_n = 0;
    return
end

q = br(:, 1);
E = br(:, 2);
G = br(:, 3);
q_abs = abs(q);
stats.q_min = min(q);
stats.q_max = max(q);
stats.q_abs_min = min(q_abs);
stats.q_abs_max = max(q_abs);
stats.energy_min = min(E);
stats.energy_max = max(E);
center_mask = abs(q) <= 0.02;
edge_mask = abs(q) >= 0.12 & abs(q) <= 0.15;
inner_mask = q_abs <= stats.q_abs_min + max(dq_Ainv, 0.005);
stats.energy_center_mean = mean(E(center_mask), 'omitnan');
stats.energy_inner_mean = mean(E(inner_mask), 'omitnan');
stats.energy_edge_mean = mean(E(edge_mask), 'omitnan');
stats.energy_edge_minus_center = stats.energy_edge_mean - stats.energy_center_mean;
stats.energy_edge_minus_inner = stats.energy_edge_mean - stats.energy_inner_mean;
stats.gamma_median = median(G, 'omitnan');
stats.gamma_ratio_median = median(G ./ max(E, eps), 'omitnan');
if size(br, 2) >= 7
    ci_half = 0.5 * (br(:, 7) - br(:, 6));
else
    ci_half = NaN(size(E));
end
stats.ci_half_median = median(ci_half, 'omitnan');
stats.ci_half_max = max(ci_half, [], 'omitnan');
sym = local_symmetry_stats(q, E, dq_Ainv);
stats.sym_pair_n = sym.n_pairs;
stats.sym_mean_abs_delta = sym.mean_abs_delta;
stats.sym_max_abs_delta = sym.max_abs_delta;
end


function sym = local_symmetry_stats(q, E, dq_Ainv)
q_abs = round(abs(q(:)) ./ dq_Ainv) .* dq_Ainv;
levels = unique(q_abs(q_abs > 0));
deltas = [];
tol = max(dq_Ainv * 0.1, 1e-9);
for i = 1:numel(levels)
    level = levels(i);
    neg = q < 0 & abs(abs(q) - level) <= tol;
    pos = q > 0 & abs(abs(q) - level) <= tol;
    if any(neg) && any(pos)
        deltas(end+1, 1) = mean(E(pos), 'omitnan') - ...
            mean(E(neg), 'omitnan'); %#ok<AGROW>
    end
end
sym = struct();
sym.n_pairs = numel(deltas);
if isempty(deltas)
    sym.mean_abs_delta = NaN;
    sym.max_abs_delta = NaN;
else
    sym.mean_abs_delta = mean(abs(deltas), 'omitnan');
    sym.max_abs_delta = max(abs(deltas), [], 'omitnan');
end
end


function [best_model, best_bic] = local_best_model(entries)
best_model = '';
best_bic = NaN;
if isempty(entries)
    return
end
bic = [entries.bic];
ok = [entries.success] & isfinite(bic);
if ~any(ok)
    return
end
idx_ok = find(ok);
[best_bic, rel] = min(bic(ok));
best_model = entries(idx_ok(rel)).model;
end


function tbl = local_model_summary_table(model_results)
rows = {};
for b = 1:numel(model_results)
    entries = model_results{b};
    for i = 1:numel(entries)
        fit = entries(i).fit;
        n = NaN;
        k = NaN;
        r2 = NaN;
        rmse = NaN;
        p1 = NaN;
        p2 = NaN;
        eflat = NaN;
        rho0 = NaN;
        qc = NaN;
        err = entries(i).error;
        if entries(i).success
            n = numel(fit.E_data);
            k = numel(fit.params);
            r2 = fit.R_squared;
            rmse = fit.RMSE_meV;
            p1 = fit.params(1);
            if numel(fit.params) >= 2
                p2 = fit.params(2);
            end
            if isfield(fit, 'E_flat_meV')
                eflat = fit.E_flat_meV;
            end
            if isfield(fit, 'rho0')
                rho0 = fit.rho0;
            end
            if isfield(fit, 'q_c_Ainv')
                qc = fit.q_c_Ainv;
            end
        end
        rows(end+1, :) = {b, entries(i).model, entries(i).success, ...
            n, k, r2, rmse, entries(i).aic, entries(i).bic, ...
            p1, p2, eflat, rho0, qc, err}; %#ok<AGROW>
    end
end
tbl = cell2table(rows, 'VariableNames', { ...
    'branch', 'model', 'success', 'n_points', 'n_params', ...
    'R2', 'RMSE_meV', 'AIC', 'BIC', 'param1', 'param2', ...
    'E_flat_meV', 'rho0_A', 'q_c_Ainv', 'error'});
end


function tbl = local_single_spectrum_summary(qe_pp, fit_res, snap)
q_targets = unique([snap.selectedQ_Ainv, 0.0, 0.05, 0.10, 0.145], 'stable');
rows = {};
for i = 1:numel(q_targets)
    [~, qi] = min(abs(qe_pp.q_Ainv - q_targets(i)));
    detail = fit_res.fit_details{qi};
    if isempty(detail)
        rows(end+1, :) = {q_targets(i), qe_pp.q_Ainv(qi), qi, ...
            0, NaN, NaN, NaN, NaN, ''}; %#ok<AGROW>
        continue
    end
    if isfield(detail, 'apex_energy_meV')
        peak_energy = detail.apex_energy_meV(:);
    else
        peak_energy = detail.omega_p(:);
    end
    fano_q = NaN;
    if isfield(detail, 'fano_q') && ~isempty(detail.fano_q)
        fano_q = median(detail.fano_q, 'omitnan');
    end
    apex_offset = NaN;
    if isfield(detail, 'apex_offset_meV') && ~isempty(detail.apex_offset_meV)
        apex_offset = median(detail.apex_offset_meV, 'omitnan');
    end
    note = '';
    if isfield(detail, 'peak_quality_notes') && ~isempty(detail.peak_quality_notes)
        note = strjoin(cellstr(string(detail.peak_quality_notes)), '; ');
    end
    rows(end+1, :) = {q_targets(i), qe_pp.q_Ainv(qi), qi, ...
        detail.n_peaks, min(peak_energy), max(peak_energy), ...
        detail.R_squared, fano_q, apex_offset, note}; %#ok<AGROW>
end
tbl = cell2table(rows, 'VariableNames', { ...
    'target_q_Ainv', 'nearest_q_Ainv', 'q_index', 'n_peaks', ...
    'min_peak_energy_meV', 'max_peak_energy_meV', 'R2', ...
    'median_fano_q', 'median_apex_offset_meV', 'quality_notes'});
end


function local_write_branch_tables(branches, out_dir)
for b = 1:numel(branches)
    br = branches{b};
    tbl = local_branch_points_table(br, b);
    writetable(tbl, fullfile(out_dir, sprintf('branch%d_points.csv', b)));
end
end


function tbl = local_branch_points_table(br, branch_id)
cols = {'q_Ainv', 'energy_meV', 'gamma_meV', 'R2', 'amplitude_fit', ...
    'E_ci_lo', 'E_ci_hi', 'gamma_ci_lo', 'gamma_ci_hi', ...
    'A_ci_lo', 'A_ci_hi', 'raw_height'};
if isempty(br)
    tbl = cell2table(cell(0, numel(cols) + 2), ...
        'VariableNames', [cols, {'branch', 'E_ci_half_meV'}]);
    return
end
if size(br, 2) < numel(cols)
    br(:, end+1:numel(cols)) = NaN;
end
tbl = array2table(br(:, 1:numel(cols)), 'VariableNames', cols);
tbl.branch = repmat(branch_id, size(br, 1), 1);
tbl.E_ci_half_meV = 0.5 * (tbl.E_ci_hi - tbl.E_ci_lo);
end


function figure_paths = local_export_figures(qe_pp, branches, ...
    model_results, fit_res, snap, out_dir)
figure_paths = struct();
figure_paths.qe_map = fullfile(out_dir, 'qe_area_map_with_branches.png');
figure_paths.dispersion = fullfile(out_dir, 'dispersion_area_fano_branches.png');
figure_paths.single_spectrum_selected = fullfile(out_dir, ...
    'single_spectrum_selected_q_fit.png');
figure_paths.single_spectrum_highq = fullfile(out_dir, ...
    'single_spectrum_high_q_fit.png');

local_plot_qe_map(qe_pp, branches, snap, figure_paths.qe_map, ...
    'Area-normalized physical q-E map (GUI top-left)');
local_plot_dispersion(branches, model_results, figure_paths.dispersion);
local_plot_single_spectrum(qe_pp, fit_res, snap.selectedQ_Ainv, snap, ...
    figure_paths.single_spectrum_selected);
local_plot_single_spectrum(qe_pp, fit_res, 0.145, snap, ...
    figure_paths.single_spectrum_highq);
end


function local_plot_qe_map(qe, branches, snap, out_path, title_text)
q_mask = qe.q_Ainv >= min(snap.qStart, snap.qEnd) & ...
    qe.q_Ainv <= max(snap.qStart, snap.qEnd);
e_mask = qe.energy_meV >= min(snap.energyMin, snap.energyMax) & ...
    qe.energy_meV <= max(snap.energyMin, snap.energyMax);
map = double(qe.intensity(e_mask, q_mask));
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 520]);
ax = axes(fig);
imagesc(ax, qe.q_Ainv(q_mask), qe.energy_meV(e_mask), map);
axis(ax, 'xy');
colormap(ax, turbo);
color_limits = local_color_limits(map);
if all(isfinite(color_limits)) && color_limits(1) < color_limits(2)
    clim(ax, color_limits);
end
colorbar(ax);
hold(ax, 'on');
for b = 1:numel(branches)
    br = branches{b};
    if isempty(br)
        continue
    end
    qe_plot_helpers.plot_branch_scatter(ax, br, ...
        qe_plot_helpers.branch_color(b), sprintf('Branch %d', b));
end
hold(ax, 'off');
xlabel(ax, 'q (1/A)');
ylabel(ax, 'Energy relative to ZLP (meV)');
title(ax, title_text);
legend(ax, 'Location', 'best', 'FontSize', 7);
grid(ax, 'on');
exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function clim = local_color_limits(map)
vals = map(isfinite(map));
if isempty(vals)
    clim = [NaN NaN];
    return
end
vals = sort(vals(:));
lo = local_percentile(vals, 2);
hi = local_percentile(vals, 98);
if lo == hi
    hi = lo + eps;
end
clim = [lo hi];
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


function local_plot_dispersion(branches, model_results, out_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 560]);
ax = axes(fig);
hold(ax, 'on');
for b = 1:numel(branches)
    br = branches{b};
    if isempty(br)
        continue
    end
    col = qe_plot_helpers.branch_color(b);
    qe_plot_helpers.plot_branch_scatter(ax, br, col, sprintf('Branch %d', b));
    best = local_best_model_entry(model_results{b});
    if ~isempty(best)
        local_plot_supported_fit_curve(ax, best.fit, col, b);
    end
end
hold(ax, 'off');
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, 'q (1/A)');
ylabel(ax, 'Energy (meV)');
title(ax, 'Area-normalized Fano apex dispersion');
legend(ax, 'Location', 'best', 'FontSize', 8);
exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function local_plot_supported_fit_curve(ax, fit, col, branch_idx)
q_data = abs(fit.q_data(:));
q_data = q_data(isfinite(q_data) & q_data > 0);
if isempty(q_data)
    return
end
q_min = min(q_data);
q_max = max(q_data);
keep = abs(fit.q_fit) >= q_min & abs(fit.q_fit) <= q_max;
if ~any(keep)
    return
end
label = sprintf('Fit B%d: %s  R^2=%.3f', ...
    branch_idx, char(string(fit.model_name)), fit.R_squared);
plotted_label = false;
for side = [-1, 1]
    if side < 0
        side_keep = keep & fit.q_fit < 0;
    else
        side_keep = keep & fit.q_fit > 0;
    end
    if ~any(side_keep)
        continue
    end
    if plotted_label
        plot(ax, fit.q_fit(side_keep), fit.E_fit(side_keep), '-', ...
            'Color', col, 'LineWidth', 2.0, 'HandleVisibility', 'off');
    else
        plot(ax, fit.q_fit(side_keep), fit.E_fit(side_keep), '-', ...
            'Color', col, 'LineWidth', 2.0, 'DisplayName', label);
        plotted_label = true;
    end
end
end


function entry = local_best_model_entry(entries)
entry = [];
if isempty(entries)
    return
end
bic = [entries.bic];
ok = [entries.success] & isfinite(bic);
if ~any(ok)
    return
end
idx_ok = find(ok);
[~, rel] = min(bic(ok));
entry = entries(idx_ok(rel));
end


function local_plot_single_spectrum(qe_pp, fit_res, target_q, snap, out_path)
[~, qi] = min(abs(qe_pp.q_Ainv - target_q));
[mask, energy_axis] = local_energy_mask(qe_pp, snap);
spectrum = double(qe_pp.intensity(mask, qi));
detail = fit_res.fit_details{qi};
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 520]);
ax = axes(fig);
plot(ax, energy_axis, spectrum, 'Color', [0.1 0.35 0.75], 'LineWidth', 1.1, ...
    'DisplayName', 'area-normalized spectrum');
hold(ax, 'on');
if ~isempty(detail)
    plot(ax, detail.energy_fit, detail.curve_fit, 'r-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%s fit', detail.peak_model_name));
    for p = 1:numel(detail.peak_curves)
        plot(ax, detail.energy_fit, detail.peak_curves{p}, '--', ...
            'LineWidth', 0.9, 'DisplayName', sprintf('peak %d', p));
    end
    if isfield(detail, 'apex_energy_meV')
        yl = ylim(ax);
        for p = 1:numel(detail.apex_energy_meV)
            xline(ax, detail.apex_energy_meV(p), ':', ...
                sprintf('%.0f', detail.apex_energy_meV(p)), ...
                'HandleVisibility', 'off');
        end
        ylim(ax, yl);
    end
end
hold(ax, 'off');
grid(ax, 'on');
xlabel(ax, 'Energy relative to ZLP (meV)');
ylabel(ax, 'Area-normalized intensity');
title(ax, sprintf('Single spectrum fit | q = %.4f 1/A | index %d', ...
    qe_pp.q_Ainv(qi), qi));
legend(ax, 'Location', 'best', 'FontSize', 8);
exportgraphics(fig, out_path, 'Resolution', 300);
close(fig);
end


function report_path = local_write_report(out_dir, dataset, history_entry, snap, ...
    pp_opts, fit_res, branch_summary, model_summary, ...
    single_summary, figure_paths, refinement_log, refinement_profile)
report_path = fullfile(out_dir, 'analysis_report.md');
fid = fopen(report_path, 'w', 'n', 'UTF-8');
if fid < 0
    error('run_590_gui_history_area_analysis:CannotWriteReport', ...
        'Cannot write %s', report_path);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '# %s area-normalized GUI-history analysis\n\n', dataset.label);
fprintf(fid, '## Reproducibility\n\n');
fprintf(fid, '- Dataset: `%s`\n', dataset.label);
fprintf(fid, '- Source: `%s`\n', dataset.source_path);
fprintf(fid, '- History entry: `%s`\n', char(string(history_entry.label)));
fprintf(fid, '- History source: `%s` (`%s`).\n', ...
    history_entry.source_path, history_entry.source_note);
fprintf(fid, '- Normalization forced for this analysis: `%s`, window `[%.0f, %.0f] meV`.\n', ...
    pp_opts.norm_method, pp_opts.norm_min, pp_opts.norm_max);
fprintf(fid, '- Denoise: `%s` via `%s`; background subtraction: `%s`; peak model: `%s`.\n', ...
    local_on_off(pp_opts.do_denoise), pp_opts.denoise_method, ...
    local_on_off(pp_opts.do_bg_sub), char(snap.peakModel));
fprintf(fid, '- Fit range: q `[%.3f, %.3f] 1/A`, E `[%.0f, %.0f] meV`; prominence `%.3f`; max peaks `%d`.\n\n', ...
    snap.qStart, snap.qEnd, snap.energyMin, snap.energyMax, ...
    snap.prominence, snap.maxPeaks);

fprintf(fid, '## Main quantitative results\n\n');
fprintf(fid, '- Auto-fit produced `%d` fitted peaks after R2 filtering; branch assignment kept `%d` points across `%d` GUI branches.\n', ...
    size(fit_res.all_peaks, 1), sum(branch_summary.n_points), height(branch_summary));
fprintf(fid, '- Heatmap evidence uses the GUI top-left physical q-E map. The lower-left comparison map is not used or exported in this report.\n');
fprintf(fid, '- The extraction uses old blind-window logic for B1/B3 and window-seeded propagation for B2, matching the current GUI hybrid logic.\n');
fprintf(fid, '- Energy CIs in the branch tables are apex-energy CIs when available; otherwise the GUI fallback floor is used.\n\n');

if isstruct(refinement_profile) && isfield(refinement_profile, 'enabled') ...
        && refinement_profile.enabled
    n_replaced = sum(strcmp(refinement_log.action, 'replaced'));
    n_rejected = sum(strcmp(refinement_log.action, 'rejected'));
    fprintf(fid, '- Dataset-specific high-q refinement: `%s`; target B%d at |q| >= %.3f 1/A, local refit window `[%.0f, %.0f] meV`.\n', ...
        refinement_profile.label, refinement_profile.branch_index, ...
        refinement_profile.q_min_Ainv, ...
        refinement_profile.refit_window_meV(1), ...
        refinement_profile.refit_window_meV(2));
    fprintf(fid, '- High-q refinement changed `%d` points and rejected `%d` low-confidence points; see `branch_refinement_log.csv`.\n\n', ...
        n_replaced, n_rejected);
end

fprintf(fid, '### Branch summary\n\n');
fprintf(fid, '| Branch | N | |q| coverage | Energy range (meV) | Center mean (meV) | Inner mean (meV) | Edge mean (meV) | Edge-inner (meV) | Median gamma/E | Median CI half-width (meV) | Symmetry mean | Best BIC model |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n');
for i = 1:height(branch_summary)
    fprintf(fid, '| B%d | %d | %.3f-%.3f | %.0f-%.0f | %.1f | %.1f | %.1f | %.1f | %.3f | %.2f | %.1f | %s |\n', ...
        branch_summary.branch(i), branch_summary.n_points(i), ...
        branch_summary.q_abs_min_Ainv(i), branch_summary.q_abs_max_Ainv(i), ...
        branch_summary.energy_min_meV(i), branch_summary.energy_max_meV(i), ...
        branch_summary.energy_center_mean_meV(i), ...
        branch_summary.energy_inner_mean_meV(i), ...
        branch_summary.energy_edge_mean_meV(i), ...
        branch_summary.energy_edge_minus_inner_meV(i), ...
        branch_summary.gamma_over_E_median(i), ...
        branch_summary.energy_ci_half_median_meV(i), ...
        branch_summary.symmetry_mean_abs_delta_meV(i), ...
        char(string(branch_summary.best_model_by_BIC{i})));
end
fprintf(fid, '\n');

fprintf(fid, '### Dispersion model comparison\n\n');
fprintf(fid, '| Branch | Model | R2 | RMSE (meV) | BIC | E_flat (meV) | rho0 (A) |\n');
fprintf(fid, '|---:|---|---:|---:|---:|---:|---:|\n');
for i = 1:height(model_summary)
    if ~model_summary.success(i)
        continue
    end
    fprintf(fid, '| B%d | %s | %.3f | %.1f | %.1f | %.1f | %.2f |\n', ...
        model_summary.branch(i), char(string(model_summary.model{i})), ...
        model_summary.R2(i), model_summary.RMSE_meV(i), ...
        model_summary.BIC(i), model_summary.E_flat_meV(i), ...
        model_summary.rho0_A(i));
end
fprintf(fid, '\n');

fprintf(fid, '### Single-spectrum checkpoints\n\n');
fprintf(fid, '| Target q | Nearest q | q index | Peaks | Energy span (meV) | R2 | Median Fano q | Median apex offset (meV) |\n');
fprintf(fid, '|---:|---:|---:|---:|---:|---:|---:|---:|\n');
for i = 1:height(single_summary)
    fprintf(fid, '| %.4f | %.4f | %d | %d | %.0f-%.0f | %.3f | %.2f | %.1f |\n', ...
        single_summary.target_q_Ainv(i), single_summary.nearest_q_Ainv(i), ...
        single_summary.q_index(i), single_summary.n_peaks(i), ...
        single_summary.min_peak_energy_meV(i), ...
        single_summary.max_peak_energy_meV(i), single_summary.R2(i), ...
        single_summary.median_fano_q(i), ...
        single_summary.median_apex_offset_meV(i));
end
fprintf(fid, '\n');

fprintf(fid, '## Discussion points for advisor meeting\n\n');
local_write_discussion_points(fid, branch_summary, model_summary);

fprintf(fid, '## Generated figures\n\n');
names = fieldnames(figure_paths);
for i = 1:numel(names)
    fprintf(fid, '- `%s`: `%s`\n', names{i}, figure_paths.(names{i}));
end
fprintf(fid, '\n');

fprintf(fid, '## Output tables\n\n');
fprintf(fid, '- `branch1_points.csv`, `branch2_points.csv`, `branch3_points.csv`\n');
fprintf(fid, '- `branch_assignment_summary.csv`\n');
fprintf(fid, '- `rejected_peaks.csv` when any auto-fit points were rejected\n');
fprintf(fid, '- `branch_refinement_log.csv`\n');
fprintf(fid, '- `branch_summary.csv`\n');
fprintf(fid, '- `dispersion_model_summary.csv`\n');
fprintf(fid, '- `single_spectrum_fit_summary.csv`\n');
fprintf(fid, '- `analysis_results.mat`\n');
end


function local_write_discussion_points(fid, branch_summary, model_summary)
for b = 1:height(branch_summary)
    best = branch_summary.best_model_by_BIC{b};
    n = branch_summary.n_points(b);
    edge_shift = branch_summary.energy_edge_minus_inner_meV(b);
    sym_delta = branch_summary.symmetry_mean_abs_delta_meV(b);
    ci_med = branch_summary.energy_ci_half_median_meV(b);
    fprintf(fid, '- B%d: %d points. The measured |q| support is %.3f-%.3f 1/A; edge-inner shift is %.1f meV, median CI half-width is %.2f meV, and +/-q mean mismatch is %.1f meV. Best BIC model in this small suite: `%s`.\n', ...
        b, n, branch_summary.q_abs_min_Ainv(b), ...
        branch_summary.q_abs_max_Ainv(b), edge_shift, ci_med, ...
        sym_delta, char(string(best)));
    if branch_summary.q_abs_min_Ainv(b) > 0.02
        fprintf(fid, '  - B%d has no near-center points inside |q| <= 0.02 1/A, so do not interpret the model curve as a measured q -> 0 extrapolation.\n', b);
    end
end

q2d = model_summary(strcmp(model_summary.model, 'quasi2d_plasmon') & ...
    model_summary.success, :);
for i = 1:height(q2d)
    if isfinite(q2d.E_flat_meV(i))
        fprintf(fid, '- Quasi-2D fit B%d gives E_flat = %.1f meV and rho0 = %.2f A. Treat these as effective parameters of the selected empirical model, not intrinsic constants.\n', ...
            q2d.branch(i), q2d.E_flat_meV(i), q2d.rho0_A(i));
    end
end
fprintf(fid, '- Because Area normalization was used, branch positions and widths are the reliable outputs; branch amplitudes should not be used as physical A(q) evidence in this run.\n');
end


function value = local_on_off(tf)
if tf
    value = 'on';
else
    value = 'off';
end
end
