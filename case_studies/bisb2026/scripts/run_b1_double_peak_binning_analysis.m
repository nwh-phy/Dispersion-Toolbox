function output = run_b1_double_peak_binning_analysis(sessionRequest, options)
%RUN_B1_DOUBLE_PEAK_BINNING_ANALYSIS Mandatory B1 double-peak extraction.
%
% This workflow keeps the old single-peak B1 exports as comparison inputs.
% New outputs are written to b1_double_peak/combined-q/binning directories.

arguments
    sessionRequest {mustBeTextScalar} = "all"
    options.runFits (1,1) logical = true
    options.outputDateTag {mustBeTextScalar} = "260508"
    options.qRangeOverride_Ainv (1,2) double = [-0.15 0.15]
    options.b1EnergyWindowOverrideMeV (1,2) double = [300 2100]
    options.waterfallStartMeV (1,1) double = 0
    options.waterfallEndMeV (1,1) double = NaN
    options.waterfallNormMode {mustBeTextScalar} = "visual"
    options.waterfallAreaNormWindowMeV (1,2) double = [250 3800]
    options.waterfallResidual (1,1) logical = false
    options.waterfallGain (1,1) double = 1
    options.fitDenoiseMethod {mustBeTextScalar} = "none"
    options.fitDenoiseProfile {mustBeTextScalar} = "global"
    options.fitDenoiseWindow (1,1) double = 11
    options.fitDenoiseLowWindow (1,1) double = 11
    options.fitDenoiseHighWindow (1,1) double = 31
    options.fitDenoiseQStartAinv (1,1) double = 0.07
    options.fitDenoiseQEndAinv (1,1) double = 0.15
    options.fitDenoiseOrder (1,1) double = 3
    options.highQForceBinAbsAinv (1,1) double = Inf
    options.binSize (1,1) double = 3
    options.peakModelOverride {mustBeTextScalar} = ""
    options.trackingMode {mustBeTextScalar} = "independent_double_peak"
    options.fallbackSplitCandidatesMeV (1,:) double = NaN
    options.maxTrackingShiftMeV (1,1) double = 180
    options.trackingWindowHalfWidthMeV (1,1) double = 220
    options.trackingWindowHighQHalfWidthMeV (1,1) double = 300
    options.trackingWindowHighQAbsAinv (1,1) double = 0.09
    options.trackingWindowInvalidFallback {mustBeTextScalar} = "fail"
    options.upperTrackingWindowHalfWidthMeV (1,1) double = NaN
    options.upperTrackingWindowHighQHalfWidthMeV (1,1) double = NaN
    options.upperTrackingWindowHighQAbsAinv (1,1) double = NaN
    options.upperQualityRetry (1,1) logical = false
    options.upperMaxGammaOverE (1,1) double = Inf
    options.upperMaxGammaMeV (1,1) double = Inf
    options.upperRetryWindowHalfWidthMeV (1,1) double = 180
    options.upperRetryHighQHalfWidthMeV (1,1) double = 240
    options.upperRetryHighQAbsAinv (1,1) double = 0.09
    options.enableJumpRepair (1,1) logical = false
    options.largeJumpThresholdMeV (1,1) double = 250
    options.referenceLowerPoints table = table()
    options.referenceUpperPoints table = table()
end

script_dir = fileparts(mfilename('fullpath'));
project_root = bisb_find_project_root(script_dir);
addpath(script_dir);
startup_path = fullfile(project_root, 'startup.m');
if isfile(startup_path)
    run(startup_path);
end

sessions = local_requested_sessions(project_root, sessionRequest, options);
session_outputs = cell(1, numel(sessions));
for i = 1:numel(sessions)
    session_outputs{i} = local_run_one_session(project_root, sessions(i), ...
        options);
end

fit_outputs = struct();
if options.runFits
    fit_datasets = local_fit_datasets(session_outputs);
    fit_outputs.lower_physical = run_b1_physical_fit_analysis( ...
        datasets=fit_datasets, ...
        branchFileName='b1_double_peak_lower_points.csv', ...
        outputTag='b1_double_peak_lower_physical_fit_binning_260508', ...
        filePrefix='b1_double_peak_lower_binning');
    fit_outputs.upper_physical = run_b1_physical_fit_analysis( ...
        datasets=fit_datasets, ...
        branchFileName='b1_double_peak_upper_points.csv', ...
        outputTag='b1_double_peak_upper_physical_fit_binning_260508', ...
        filePrefix='b1_double_peak_upper_binning');
    fit_outputs.lower_enhancements = run_b1_physical_fit_enhancements( ...
        datasets=fit_datasets, ...
        branchFileName='b1_double_peak_lower_points.csv', ...
        outputTag='b1_double_peak_lower_physical_fit_enhancements_binning_260508', ...
        filePrefix='b1_double_peak_lower_binning');
    fit_outputs.upper_enhancements = run_b1_physical_fit_enhancements( ...
        datasets=fit_datasets, ...
        branchFileName='b1_double_peak_upper_points.csv', ...
        outputTag='b1_double_peak_upper_physical_fit_enhancements_binning_260508', ...
        filePrefix='b1_double_peak_upper_binning');
end

output = struct();
output.sessions = session_outputs;
output.fit_outputs = fit_outputs;

fprintf('B1 double-peak binning analysis complete.\n');
for i = 1:numel(session_outputs)
    fprintf('  %s: %s\n', session_outputs{i}.session.session_key, ...
        session_outputs{i}.output_dir);
end
end


function sessions = local_requested_sessions(project_root, sessionRequest, options)
all_sessions = local_session_configs(project_root, char(options.outputDateTag));
request = lower(strtrim(char(string(sessionRequest))));
switch request
    case {'all', '*'}
        keep = true(size(all_sessions));
    case {'590', '590_pl2_10w', 'primary'}
        keep = strcmp({all_sessions.session_key}, '590_PL2_10w');
    case {'n0', 'n0_pl2_10w_repeat'}
        keep = strcmp({all_sessions.session_key}, 'n0_PL2_10w_repeat');
    case {'20w', 'no_pl2_20w_2film'}
        keep = strcmp({all_sessions.session_key}, 'no_PL2_20w_2film');
    otherwise
        keys = lower(string({all_sessions.session_key}));
        keep = keys == string(request);
end
sessions = all_sessions(keep);
if isempty(sessions)
    error('run_b1_double_peak_binning_analysis:UnknownSession', ...
        'Unknown session request "%s".', char(string(sessionRequest)));
end

for i = 1:numel(sessions)
    sessions(i).input_dir = local_existing_input_dir(project_root, ...
        sessions(i).preferred_input_tag, sessions(i).fallback_input_tag);
end
end


function sessions = local_session_configs(project_root, date_tag)
results_root = fullfile(project_root, 'paper_results');
sessions = repmat(local_empty_session(), 1, 3);

sessions(1) = local_make_session('590_PL2_10w', ...
    '590 10w defocus 1film', ...
    '590_gui_history_area_260506_wideq030', ...
    '590_gui_history_area_260506', ...
    sprintf('590_gui_history_area_260506_b1_double_peak_binning_%s', date_tag), ...
    '1film', 1, [0.120, 0.470, 0.900], 'o');

sessions(2) = local_make_session('n0_PL2_10w_repeat', ...
    'n0 10w defocus repeat 1film', ...
    'n0_PL2_10w_gui_history_area_260506_wideq030', ...
    'n0_PL2_10w_gui_history_area_260506', ...
    sprintf('n0_PL2_10w_gui_history_area_260506_b1_double_peak_binning_%s', date_tag), ...
    '1film', 1, [0.160, 0.500, 0.220], 'o');

sessions(3) = local_make_session('no_PL2_20w_2film', ...
    '20w defocus 2film', ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined_wideq030', ...
    'no_PL2_20w_2film_gui_history_area_260506_highq_refined', ...
    sprintf('no_PL2_20w_2film_gui_history_area_260506_highq_refined_b1_double_peak_binning_%s', date_tag), ...
    '2film', 2, [0.930, 0.280, 0.300], 'o');

for i = 1:numel(sessions)
    sessions(i).results_root = results_root;
end
end


function session = local_empty_session()
session = struct('session_key', '', 'session_label', '', ...
    'preferred_input_tag', '', 'fallback_input_tag', '', ...
    'input_dir', '', 'output_tag', '', 'results_root', '', ...
    'thickness_class', '', 'thickness_factor', NaN, ...
    'color', [0 0 0], 'marker', 'o');
end


function session = local_make_session(key, label, preferred_tag, ...
    fallback_tag, output_tag, thickness_class, thickness_factor, color, marker)
session = local_empty_session();
session.session_key = key;
session.session_label = label;
session.preferred_input_tag = preferred_tag;
session.fallback_input_tag = fallback_tag;
session.output_tag = output_tag;
session.thickness_class = thickness_class;
session.thickness_factor = thickness_factor;
session.color = color;
session.marker = marker;
end


function input_dir = local_existing_input_dir(project_root, preferred_tag, fallback_tag)
results_root = fullfile(project_root, 'paper_results');
preferred = fullfile(results_root, preferred_tag);
fallback = fullfile(results_root, fallback_tag);
if isfile(fullfile(preferred, 'analysis_results.mat'))
    input_dir = preferred;
elseif isfile(fullfile(fallback, 'analysis_results.mat'))
    input_dir = fallback;
else
    error('run_b1_double_peak_binning_analysis:MissingInput', ...
        'Missing analysis_results.mat in %s and %s.', preferred, fallback);
end
end


function out = local_run_one_session(project_root, session, options)
input_dir = session.input_dir;
result_data = load(fullfile(input_dir, 'analysis_results.mat'), 'output');
old_branch_path = fullfile(input_dir, 'branch1_points.csv');
if ~isfile(old_branch_path)
    error('run_b1_double_peak_binning_analysis:MissingB1', ...
        'Missing old B1 single-peak table: %s', old_branch_path);
end
old_points = readtable(old_branch_path);

qe_pp = result_data.output.qe_pp;
qe_raw = local_reconstruct_raw_qe(result_data.output, qe_pp);
snap = result_data.output.snap;

extract_opts = local_extract_options(snap, old_points, options);
extract = b1_double_peak_binning_extract(qe_pp, qe_raw, extract_opts);
if options.enableJumpRepair
    repaired = b1_double_peak_repair_tracking_points( ...
        extract.lower_points, extract.upper_points, extract.fit_failures, ...
        largeJumpThresholdMeV=options.largeJumpThresholdMeV);
    if isfield(extract, 'repair_log')
        extract.repair_log = [extract.repair_log; repaired.repair_log];
    else
        extract.repair_log = repaired.repair_log;
    end
    extract.exclusion_points = repaired.exclusion_points;
    extract.lower_points = repaired.lower_points;
    extract.upper_points = repaired.upper_points;
    extract.combined_points = repaired.combined_points;
elseif ~isfield(extract, 'exclusion_points')
    extract.exclusion_points = table();
end
if ~isfield(extract, 'repair_log')
    extract.repair_log = table();
end

out_dir = fullfile(project_root, 'paper_results', session.output_tag);
if ~isfolder(out_dir)
    mkdir(out_dir);
end

combined_csv = fullfile(out_dir, 'b1_double_peak_combined_q_points.csv');
lower_csv = fullfile(out_dir, 'b1_double_peak_lower_points.csv');
upper_csv = fullfile(out_dir, 'b1_double_peak_upper_points.csv');
binning_csv = fullfile(out_dir, 'b1_double_peak_binning_map.csv');
noise_csv = fullfile(out_dir, 'b1_double_peak_noise_profile.csv');
failures_csv = fullfile(out_dir, 'b1_double_peak_fit_failures.csv');
repair_log_csv = fullfile(out_dir, 'b1_double_peak_repair_log.csv');
exclusion_csv = fullfile(out_dir, 'b1_double_peak_exclusion_points.csv');
writetable(extract.combined_points, combined_csv);
writetable(extract.lower_points, lower_csv);
writetable(extract.upper_points, upper_csv);
writetable(extract.binning_map, binning_csv);
writetable(extract.noise_profile, noise_csv);
writetable(extract.fit_failures, failures_csv);
writetable(extract.repair_log, repair_log_csv);
writetable(extract.exclusion_points, exclusion_csv);

fig_paths = struct();
fig_paths.binning_spectrum_png = fullfile(out_dir, ...
    'b1_double_peak_binning_spectrum_comparison.png');
fig_paths.binning_spectrum_pdf = fullfile(out_dir, ...
    'b1_double_peak_binning_spectrum_comparison.pdf');
fig_paths.single_vs_double_png = fullfile(out_dir, ...
    'b1_single_peak_vs_double_peak_comparison.png');
fig_paths.single_vs_double_pdf = fullfile(out_dir, ...
    'b1_single_peak_vs_double_peak_comparison.pdf');
fig_paths.fit_denoise_comparison_png = fullfile(out_dir, ...
    'b1_double_peak_fit_spectrum_denoise_comparison.png');
fig_paths.fit_denoise_comparison_pdf = fullfile(out_dir, ...
    'b1_double_peak_fit_spectrum_denoise_comparison.pdf');
fig_paths.double_peak_heatmap_png = fullfile(out_dir, ...
    'b1_double_peak_heatmap_overlay.png');
fig_paths.double_peak_heatmap_pdf = fullfile(out_dir, ...
    'b1_double_peak_heatmap_overlay.pdf');
fig_paths.stacked_signed_q_png = fullfile(out_dir, ...
    'b1_double_peak_stacked_spectra_signed_q.png');
fig_paths.stacked_signed_q_pdf = fullfile(out_dir, ...
    'b1_double_peak_stacked_spectra_signed_q.pdf');
waterfall_name = 'b1_double_peak_waterfall_signed_q';
if strcmp(extract_opts.waterfall_norm_mode, 'area')
    waterfall_name = sprintf('%s_area%s_%snorm', waterfall_name, ...
        local_meV_tag(extract_opts.waterfall_area_norm_window_meV(1)), ...
        local_meV_tag(extract_opts.waterfall_area_norm_window_meV(2)));
end
if extract_opts.waterfall_start_meV > 0
    waterfall_name = sprintf('%s_start%smeV', waterfall_name, ...
        local_meV_tag(extract_opts.waterfall_start_meV));
end
if isfinite(extract_opts.waterfall_end_meV)
    waterfall_name = sprintf('%s_end%smeV', waterfall_name, ...
        local_meV_tag(extract_opts.waterfall_end_meV));
end
if extract_opts.waterfall_residual
    waterfall_name = sprintf('%s_residual', waterfall_name);
end
if abs(extract_opts.waterfall_gain - 1) > 1e-12
    waterfall_name = sprintf('%s_gain%s', waterfall_name, ...
        local_meV_tag(extract_opts.waterfall_gain));
end
fig_paths.waterfall_signed_q_png = fullfile(out_dir, [waterfall_name '.png']);
fig_paths.waterfall_signed_q_pdf = fullfile(out_dir, [waterfall_name '.pdf']);
manual_window_seed_csv = fullfile(out_dir, ...
    'b1_double_peak_manual_window_seed.csv');
plot_q_binning_csv = fullfile(out_dir, ...
    'b1_double_peak_plot_q_binning_map.csv');

local_plot_binning_spectra(qe_pp, extract, extract_opts, ...
    fig_paths.binning_spectrum_png, fig_paths.binning_spectrum_pdf);
local_plot_single_vs_double(old_points, extract, session, ...
    fig_paths.single_vs_double_png, fig_paths.single_vs_double_pdf);
local_plot_fit_denoise_comparison(qe_pp, extract, extract_opts, session, ...
    fig_paths.fit_denoise_comparison_png, ...
    fig_paths.fit_denoise_comparison_pdf);
local_plot_double_peak_heatmap(qe_pp, old_points, extract, extract_opts, ...
    session, fig_paths.double_peak_heatmap_png, ...
    fig_paths.double_peak_heatmap_pdf);
local_plot_stacked_spectra(qe_pp, extract, old_points, extract_opts, session, ...
    'signed_q', fig_paths.stacked_signed_q_png, ...
    fig_paths.stacked_signed_q_pdf);
local_plot_waterfall_spectra(qe_pp, extract, extract_opts, session, ...
    'signed_q', fig_paths.waterfall_signed_q_png, ...
    fig_paths.waterfall_signed_q_pdf);
manual_window_seed = local_manual_window_seed_table(session, extract_opts);
writetable(manual_window_seed, manual_window_seed_csv);
plot_q_binning_map = local_plot_q_binning_map(qe_pp, extract, extract_opts);
writetable(plot_q_binning_map, plot_q_binning_csv);

mat_path = fullfile(out_dir, 'b1_double_peak_binning_results.mat');
save(mat_path, 'extract', 'extract_opts', 'session', 'input_dir', ...
    'old_branch_path', 'fig_paths', 'manual_window_seed_csv', ...
    'plot_q_binning_csv', '-v7.3');

out = struct();
out.session = session;
out.input_dir = input_dir;
out.output_dir = out_dir;
out.extract = extract;
out.combined_csv = combined_csv;
out.lower_csv = lower_csv;
out.upper_csv = upper_csv;
out.binning_csv = binning_csv;
out.noise_csv = noise_csv;
out.failures_csv = failures_csv;
out.repair_log_csv = repair_log_csv;
out.exclusion_csv = exclusion_csv;
out.manual_window_seed_csv = manual_window_seed_csv;
out.plot_q_binning_csv = plot_q_binning_csv;
out.figure_paths = fig_paths;
out.mat_path = mat_path;
end


function qe_raw = local_reconstruct_raw_qe(saved_output, qe_pp)
qe_raw = qe_pp;
if ~isfield(saved_output, 'dataset') || ~isfield(saved_output.dataset, 'qe') || ...
        ~isfield(saved_output, 'preprocess_opts')
    return
end
raw_opts = saved_output.preprocess_opts;
raw_opts.do_normalize = false;
try
    qe_raw = qe_preprocess(saved_output.dataset.qe, raw_opts);
catch
    qe_raw = qe_pp;
end
end


function opts = local_extract_options(snap, old_points, run_options)
opts = struct();
opts.energy_window_meV = sort([local_snap_value(snap, 'branch1Min', 500), ...
    local_snap_value(snap, 'branch1Max', 2100)]);
if all(isfinite(run_options.b1EnergyWindowOverrideMeV))
    opts.energy_window_meV = sort(run_options.b1EnergyWindowOverrideMeV);
end
q_range = sort([local_snap_value(snap, 'qStart', -0.15), ...
    local_snap_value(snap, 'qEnd', 0.15)]);
if all(isfinite(run_options.qRangeOverride_Ainv))
    q_range = sort(run_options.qRangeOverride_Ainv);
end
opts.q_range_Ainv = q_range;
opts.q_skip_Ainv = 0.005;
opts.bin_size = run_options.binSize;
opts.noise_threshold = NaN;
opts.low_q_no_bin_abs_Ainv = 0.05;
opts.peak_model = char(local_snap_value(snap, 'peakModel', 'fano'));
if strlength(string(run_options.peakModelOverride)) > 0
    opts.peak_model = char(string(run_options.peakModelOverride));
end
opts.pre_subtracted = logical(local_snap_value(snap, 'bgSub', false));
opts.min_prominence = local_snap_value(snap, 'prominence', 0.10);
opts.smooth_width = 1;
opts.bootstrap_ci_samples = local_snap_value(snap, 'bootstrapCiSamples', 0);
opts.old_branch_points = old_points;
opts.fallback_split_meV = 180;
if all(isfinite(run_options.fallbackSplitCandidatesMeV))
    opts.fallback_split_candidates_meV = run_options.fallbackSplitCandidatesMeV;
end
opts.min_peak_separation_meV = 10;
opts.tracking_mode = char(string(run_options.trackingMode));
opts.max_tracking_shift_meV = run_options.maxTrackingShiftMeV;
opts.tracking_window_half_width_meV = run_options.trackingWindowHalfWidthMeV;
opts.tracking_window_highq_half_width_meV = run_options.trackingWindowHighQHalfWidthMeV;
opts.tracking_window_highq_abs_Ainv = run_options.trackingWindowHighQAbsAinv;
opts.tracking_window_invalid_fallback = char(string( ...
    run_options.trackingWindowInvalidFallback));
if isfinite(run_options.upperTrackingWindowHalfWidthMeV)
    opts.upper_tracking_window_half_width_meV = ...
        run_options.upperTrackingWindowHalfWidthMeV;
end
if isfinite(run_options.upperTrackingWindowHighQHalfWidthMeV)
    opts.upper_tracking_window_highq_half_width_meV = ...
        run_options.upperTrackingWindowHighQHalfWidthMeV;
end
if isfinite(run_options.upperTrackingWindowHighQAbsAinv)
    opts.upper_tracking_window_highq_abs_Ainv = ...
        run_options.upperTrackingWindowHighQAbsAinv;
end
opts.upper_quality_retry = run_options.upperQualityRetry;
opts.upper_max_gamma_over_E = run_options.upperMaxGammaOverE;
opts.upper_max_gamma_meV = run_options.upperMaxGammaMeV;
opts.upper_retry_window_half_width_meV = ...
    run_options.upperRetryWindowHalfWidthMeV;
opts.upper_retry_highq_half_width_meV = ...
    run_options.upperRetryHighQHalfWidthMeV;
opts.upper_retry_highq_abs_Ainv = run_options.upperRetryHighQAbsAinv;
opts.reference_lower_points = run_options.referenceLowerPoints;
opts.reference_upper_points = run_options.referenceUpperPoints;
opts.waterfall_start_meV = max(0, run_options.waterfallStartMeV);
opts.waterfall_end_meV = run_options.waterfallEndMeV;
opts.waterfall_norm_mode = local_waterfall_norm_mode(run_options.waterfallNormMode);
opts.waterfall_area_norm_window_meV = sort(run_options.waterfallAreaNormWindowMeV);
opts.waterfall_residual = run_options.waterfallResidual;
opts.waterfall_gain = run_options.waterfallGain;
opts.fit_denoise_method = char(string(run_options.fitDenoiseMethod));
opts.fit_denoise_profile = char(string(run_options.fitDenoiseProfile));
opts.fit_denoise_window = run_options.fitDenoiseWindow;
opts.fit_denoise_low_window = run_options.fitDenoiseLowWindow;
opts.fit_denoise_high_window = run_options.fitDenoiseHighWindow;
opts.fit_denoise_q_start_Ainv = run_options.fitDenoiseQStartAinv;
opts.fit_denoise_q_end_Ainv = run_options.fitDenoiseQEndAinv;
opts.fit_denoise_order = run_options.fitDenoiseOrder;
opts.high_q_force_bin_abs_Ainv = run_options.highQForceBinAbsAinv;
end


function mode = local_waterfall_norm_mode(value)
text = lower(strtrim(char(string(value))));
switch text
    case {'visual', 'default'}
        mode = 'visual';
    case {'area', 'area250_3800', 'area250_3800norm', 'area_norm'}
        mode = 'area';
    otherwise
        error('run_b1_double_peak_binning_analysis:UnknownWaterfallNormMode', ...
            'Unknown waterfall normalization mode "%s".', char(string(value)));
end
end


function value = local_snap_value(snap, name, default_value)
if isfield(snap, name) && ~isempty(snap.(name))
    value = snap.(name);
else
    value = default_value;
end
end


function fit_datasets = local_fit_datasets(session_outputs)
fit_datasets = struct('session_key', {}, 'session_label', {}, ...
    'input_dir', {}, 'thickness_class', {}, 'thickness_factor', {}, ...
    'color', {}, 'marker', {});
for i = 1:numel(session_outputs)
    s = session_outputs{i}.session;
    fit_datasets(end + 1) = struct( ... %#ok<AGROW>
        'session_key', s.session_key, ...
        'session_label', s.session_label, ...
        'input_dir', session_outputs{i}.output_dir, ...
        'thickness_class', s.thickness_class, ...
        'thickness_factor', s.thickness_factor, ...
        'color', s.color, ...
        'marker', s.marker);
end
end


function local_plot_binning_spectra(qe, extract, opts, png_path, pdf_path)
energy_axis = double(qe.energy_meV(:));
energy_mask = energy_axis >= opts.energy_window_meV(1) & ...
    energy_axis <= opts.energy_window_meV(2);
combined = extract.binning_map(strcmp(extract.binning_map.source_mode, ...
    'combined_q_binning_3'), :);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1040 720]);
if isempty(combined)
    axes(fig);
    text(0.5, 0.5, 'No q bins exceeded the per-session noise threshold', ...
        'HorizontalAlignment', 'center');
    axis off;
else
    n_plot = min(height(combined), 4);
    t = tiledlayout(fig, n_plot, 1, 'TileSpacing', 'compact', ...
        'Padding', 'compact');
    for i = 1:n_plot
        ax = nexttile(t, i);
        q_idx = local_parse_index_list(combined.source_q_index{i});
        y = double(qe.intensity(energy_mask, q_idx));
        plot(ax, energy_axis(energy_mask), y, '-', 'Color', [0.72 0.72 0.72]);
        hold(ax, 'on');
        plot(ax, energy_axis(energy_mask), mean(y, 2, 'omitnan'), ...
            'k-', 'LineWidth', 1.6);
        hold(ax, 'off');
        box(ax, 'on');
        grid(ax, 'on');
        title(ax, sprintf('combined-q bin %d: q = %s A^{-1}', ...
            i, combined.source_q_Ainv{i}));
        xlabel(ax, 'Energy (meV)');
        ylabel(ax, 'Intensity');
    end
end
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_single_vs_double(old_points, extract, session, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 560]);
ax = axes(fig);
hold(ax, 'on');
scatter(ax, old_points.q_Ainv, old_points.energy_meV, 28, ...
    [0.45 0.45 0.45], 'filled', 'DisplayName', 'old B1 single peak');
scatter(ax, extract.lower_points.q_Ainv, extract.lower_points.energy_meV, ...
    30, [0.10 0.35 0.85], 'filled', 'DisplayName', 'B1 double lower');
scatter(ax, extract.upper_points.q_Ainv, extract.upper_points.energy_meV, ...
    30, [0.85 0.20 0.15], 'filled', 'DisplayName', 'B1 double upper');
hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
xlabel(ax, 'q (A^{-1})');
ylabel(ax, 'Energy (meV)');
title(ax, sprintf('B1 single-peak vs mandatory double-peak: %s', ...
    session.session_label));
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_plot_fit_denoise_comparison(qe, extract, opts, session, ...
    png_path, pdf_path)
details = extract.fit_details(~cellfun(@isempty, extract.fit_details));
details = local_flatten_fit_details(details);
has_denoise = isfield(opts, 'fit_denoise_method') && ...
    ~strcmpi(char(opts.fit_denoise_method), 'none');
if isempty(details) || ~has_denoise
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 760 360]);
    ax = axes(fig);
    text(ax, 0.5, 0.5, 'Fit-spectrum denoise is disabled for this run', ...
        'HorizontalAlignment', 'center');
    axis(ax, 'off');
    title(ax, sprintf('B1 fit-spectrum denoise comparison: %s', ...
        session.session_label), 'Interpreter', 'none');
    exportgraphics(fig, png_path, 'Resolution', 300);
    exportgraphics(fig, pdf_path, 'ContentType', 'vector');
    close(fig);
    return
end

n_show = min(6, numel(details));
pick = unique(round(linspace(1, numel(details), n_show)));
n_show = numel(pick);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1100 680]);
tiledlayout(fig, ceil(n_show / 2), 2, 'Padding', 'compact', ...
    'TileSpacing', 'compact');
energy_axis = double(qe.energy_meV(:));
for i = 1:n_show
    detail = details{pick(i)};
    ax = nexttile;
    raw = local_fit_detail_vector(detail, 'raw_unit_spectrum', energy_axis);
    fit_input = local_fit_detail_vector(detail, 'fit_input_spectrum', ...
        energy_axis);
    plot(ax, energy_axis, raw, '-', 'Color', [0.62 0.62 0.62], ...
        'DisplayName', 'original unit spectrum');
    hold(ax, 'on');
    plot(ax, energy_axis, fit_input, '-', 'Color', [0.04 0.28 0.90], ...
        'LineWidth', 1.1, 'DisplayName', 'denoised fit input');
    hold(ax, 'off');
    xlim(ax, [opts.energy_window_meV(1), opts.energy_window_meV(2)]);
    grid(ax, 'on');
    title(ax, sprintf('unit %d | %s w=%d order=%d', pick(i), ...
        detail.fit_denoise_method, detail.fit_denoise_window, ...
        detail.fit_denoise_order), 'Interpreter', 'none');
    if i == 1
        legend(ax, 'Location', 'best', 'Box', 'off');
    end
    xlabel(ax, 'Energy (meV)');
    ylabel(ax, 'Intensity');
end
sgtitle(fig, sprintf('B1 fit input denoise comparison: %s', ...
    session.session_label), 'Interpreter', 'none');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function details = local_flatten_fit_details(details)
flat = {};
for i = 1:numel(details)
    detail = details{i};
    if isstruct(detail) && isfield(detail, 'raw_unit_spectrum')
        flat{end + 1} = detail; %#ok<AGROW>
    elseif isstruct(detail) && isfield(detail, 'lower_fit')
        if isstruct(detail.lower_fit) && isfield(detail.lower_fit, ...
                'raw_unit_spectrum')
            flat{end + 1} = detail.lower_fit; %#ok<AGROW>
        end
        if isfield(detail, 'upper_fit') && isstruct(detail.upper_fit) && ...
                isfield(detail.upper_fit, 'raw_unit_spectrum')
            flat{end + 1} = detail.upper_fit; %#ok<AGROW>
        end
    end
end
details = flat;
end


function y = local_fit_detail_vector(detail, field_name, energy_axis)
if isfield(detail, field_name) && numel(detail.(field_name)) == numel(energy_axis)
    y = double(detail.(field_name)(:));
else
    y = NaN(size(energy_axis));
end
end


function local_plot_stacked_spectra(qe, extract, old_points, opts, session, mode, ...
    png_path, pdf_path)
[energy_axis, energy_mask] = local_stack_energy_window(qe, opts);
[q_values, traces, source_counts] = local_stack_trace_set(qe, extract, opts, ...
    energy_mask, mode);

switch char(mode)
    case 'positive_q'
        title_text = sprintf('B1 stacked spectra, positive q: %s', ...
            session.session_label);
        y_label = 'q (A^{-1})';
    case 'negative_q'
        title_text = sprintf('B1 stacked spectra, negative q: %s', ...
            session.session_label);
        y_label = 'q (A^{-1})';
    case 'signed_q'
        title_text = sprintf('B1 stacked spectra, signed q centered: %s', ...
            session.session_label);
        y_label = 'signed q (A^{-1}); negative below, positive above';
    case 'absq_combined'
        error('run_b1_double_peak_binning_analysis:DeprecatedAbsQCombined', ...
            ['absq_combined stacked spectra are obsolete because they ', ...
            'average +q and -q. Use signed_q instead.']);
    otherwise
        error('run_b1_double_peak_binning_analysis:UnknownStackMode', ...
            'Unknown stacked-spectrum mode "%s".', char(mode));
end

local_plot_stack_figure(energy_axis, traces, q_values, source_counts, ...
    old_points, char(mode), title_text, y_label, png_path, pdf_path);
end


function local_plot_combined_q_stacked_spectra(qe, extract, opts, session, ...
    png_path, pdf_path)
[energy_axis, energy_mask] = local_stack_energy_window(qe, opts);
combined = extract.binning_map(strcmp(extract.binning_map.source_mode, ...
    'combined_q_binning_3'), :);

traces = zeros(nnz(energy_mask), 0);
q_values = zeros(0, 1);
source_counts = zeros(0, 1);
for i = 1:height(combined)
    q_idx = local_parse_index_list(combined.source_q_index{i});
    q_idx = q_idx(q_idx >= 1 & q_idx <= numel(qe.q_Ainv));
    if isempty(q_idx)
        continue
    end
    y = mean(double(qe.intensity(energy_mask, q_idx)), 2, 'omitnan');
    traces(:, end + 1) = y; %#ok<AGROW>
    q_values(end + 1, 1) = mean(abs(double(qe.q_Ainv(q_idx))), 'omitnan'); %#ok<AGROW>
    source_counts(end + 1, 1) = numel(q_idx); %#ok<AGROW>
end

title_text = sprintf('B1 3-point combined-q stacked spectra: %s', ...
    session.session_label);
local_plot_stack_figure(energy_axis, traces, q_values, source_counts, ...
    table(), 'combined_q_binning_3', title_text, '|q| bin center (A^{-1})', ...
    png_path, pdf_path);
end


function local_plot_waterfall_spectra(qe, extract, opts, session, mode, ...
    png_path, pdf_path)
[energy_axis, energy_mask] = local_waterfall_energy_window(qe, opts);
[q_values, traces, source_counts] = local_waterfall_trace_set(qe, extract, ...
    opts, energy_mask, mode);

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [120 40 760 1120]);
ax = axes(fig);
if isempty(traces) || size(traces, 2) == 0
    text(ax, 0.5, 0.5, 'No spectra available in the requested q range', ...
        'HorizontalAlignment', 'center');
    axis(ax, 'off');
else
    [norm_energy_axis, norm_energy_mask] = local_waterfall_norm_window(qe, opts);
    [~, norm_traces] = local_waterfall_trace_set(qe, extract, opts, ...
        norm_energy_mask, mode);
    normalized = local_visual_normalize_traces(energy_axis, traces, ...
        norm_energy_axis, norm_traces, opts);
    offset = 0.55;
    offsets = (0:(size(normalized, 2) - 1)) .* offset;
    hold(ax, 'on');
    for i = 1:size(normalized, 2)
        plot(ax, energy_axis, normalized(:, i) + offsets(i), '-', ...
            'Color', local_waterfall_color(i, size(normalized, 2)), ...
            'LineWidth', 1.05);
    end
    hold(ax, 'off');
    xlim(ax, [max(0, min(energy_axis)), min(max(energy_axis), ...
        max(opts.energy_window_meV(2), 2000))]);
    ylim(ax, [-0.2, offsets(end) + 1.45]);
    ax.YTick = [];
    ax.YColor = [0 0 0];
    ax.XColor = [0 0 0];
    ax.LineWidth = 1.35;
    ax.FontSize = 18;
    ax.TickDir = 'out';
    ax.XTick = [0 1000 2000];
    xlabel(ax, 'Energy loss (meV)', 'FontSize', 24);
end
title(ax, {'B1 double-peak waterfall', ...
    char(session.session_label), local_waterfall_mode_label(mode, opts)}, ...
    'FontSize', 12, 'Interpreter', 'none');
box(ax, 'off');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function label = local_waterfall_mode_label(mode, opts)
switch char(mode)
    case 'signed_q'
        label = 'signed q';
    case 'absq_combined'
        label = 'obsolete |q| combined';
    case 'combined_q_binning_3'
        label = 'combined-q binning';
    otherwise
        label = char(mode);
end
if isfield(opts, 'waterfall_start_meV') && opts.waterfall_start_meV > 0
    label = sprintf('%s | start %s meV', label, ...
        local_meV_tag(opts.waterfall_start_meV));
end
if isfield(opts, 'waterfall_end_meV') && isfinite(opts.waterfall_end_meV)
    label = sprintf('%s | end %s meV', label, ...
        local_meV_tag(opts.waterfall_end_meV));
end
if isfield(opts, 'waterfall_norm_mode') && strcmp(opts.waterfall_norm_mode, 'area')
    label = sprintf('%s | area %s-%s norm', label, ...
        local_meV_tag(opts.waterfall_area_norm_window_meV(1)), ...
        local_meV_tag(opts.waterfall_area_norm_window_meV(2)));
end
if isfield(opts, 'waterfall_residual') && opts.waterfall_residual
    label = sprintf('%s | residual', label);
end
if isfield(opts, 'waterfall_gain') && isfinite(opts.waterfall_gain) && ...
        abs(opts.waterfall_gain - 1) > 1e-12
    label = sprintf('%s | gain %s', label, local_meV_tag(opts.waterfall_gain));
end
end


function [energy_axis, energy_mask] = local_waterfall_energy_window(qe, opts)
full_energy = double(qe.energy_meV(:));
e_min = max(opts.waterfall_start_meV, max(0, min(full_energy)));
e_max = min(max(full_energy), max(opts.energy_window_meV(2), 2000));
if isfield(opts, 'waterfall_end_meV') && isfinite(opts.waterfall_end_meV)
    e_max = min(max(full_energy), opts.waterfall_end_meV);
end
energy_mask = full_energy >= e_min & full_energy <= e_max;
energy_axis = full_energy(energy_mask);
end


function [energy_axis, energy_mask] = local_waterfall_norm_window(qe, opts)
full_energy = double(qe.energy_meV(:));
if isfield(opts, 'waterfall_area_norm_window_meV')
    window = opts.waterfall_area_norm_window_meV;
else
    window = [250 3800];
end
energy_mask = full_energy >= window(1) & full_energy <= window(2);
if ~any(energy_mask)
    energy_mask = full_energy >= 250;
end
if ~any(energy_mask)
    energy_mask = true(size(full_energy));
end
energy_axis = full_energy(energy_mask);
end


function tag = local_meV_tag(value)
tag = regexprep(sprintf('%.6g', value), '\.', 'p');
end


function [q_values, traces, source_counts] = local_waterfall_trace_set( ...
    qe, extract, opts, energy_mask, mode)
q_axis = double(qe.q_Ainv(:));

switch char(mode)
    case 'signed_q'
        groups = local_binning_map_q_groups(extract, q_axis, 'signed');
        [q_values, traces, source_counts] = local_traces_from_q_groups( ...
            qe, energy_mask, groups, false);
    case 'absq_combined'
        error('run_b1_double_peak_binning_analysis:DeprecatedAbsQCombined', ...
            ['absq_combined waterfall is obsolete because it averages ', ...
            '+q and -q. Use signed_q instead.']);
    case 'combined_q_binning_3'
        combined = extract.binning_map(strcmp(extract.binning_map.source_mode, ...
            'combined_q_binning_3'), :);
        traces = zeros(nnz(energy_mask), 0);
        q_values = zeros(0, 1);
        source_counts = zeros(0, 1);
        for i = 1:height(combined)
            q_idx = local_parse_index_list(combined.source_q_index{i});
            q_idx = q_idx(q_idx >= 1 & q_idx <= numel(q_axis));
            if isempty(q_idx)
                continue
            end
            traces(:, end + 1) = mean(double(qe.intensity(energy_mask, q_idx)), ...
                2, 'omitnan'); %#ok<AGROW>
            q_values(end + 1, 1) = mean(abs(q_axis(q_idx)), 'omitnan'); %#ok<AGROW>
            source_counts(end + 1, 1) = numel(q_idx); %#ok<AGROW>
        end
    otherwise
        error('run_b1_double_peak_binning_analysis:UnknownWaterfallMode', ...
            'Unknown waterfall mode "%s".', char(mode));
end

q_values = double(q_values(:));
source_counts = double(source_counts(:));
end


function normalized = local_visual_normalize_traces(energy_axis, traces, ...
    norm_energy_axis, norm_traces, opts)
if isfield(opts, 'waterfall_norm_mode') && strcmp(opts.waterfall_norm_mode, 'area')
    normalized = local_area_normalize_traces(traces, norm_energy_axis, ...
        norm_traces);
    normalized = local_apply_waterfall_residual_and_gain(energy_axis, ...
        normalized, opts);
    return
end
normalized = zeros(size(traces));
scale_mask = energy_axis >= 250;
if ~any(scale_mask)
    scale_mask = true(size(energy_axis));
end
for i = 1:size(traces, 2)
    y = double(traces(:, i));
    valid_scale = y(scale_mask & isfinite(y));
    if isempty(valid_scale)
        valid_scale = y(isfinite(y));
    end
    if isempty(valid_scale)
        normalized(:, i) = NaN;
        continue
    end
    baseline = prctile(valid_scale, 5);
    z = y - baseline;
    scale_values = abs(z(scale_mask & isfinite(z)));
    if isempty(scale_values)
        scale_values = abs(z(isfinite(z)));
    end
    scale = prctile(scale_values, 95);
    if ~isfinite(scale) || scale <= eps
        scale = max(scale_values, [], 'omitnan');
    end
    if ~isfinite(scale) || scale <= eps
        scale = 1;
    end
    z = z ./ scale;
    z = local_waterfall_soft_limit(z, -0.35, 1.65, 0.55);
    normalized(:, i) = z;
end
normalized = local_apply_waterfall_residual_and_gain(energy_axis, ...
    normalized, opts);
end


function normalized = local_area_normalize_traces(traces, norm_energy_axis, ...
    norm_traces)
normalized = zeros(size(traces));
for i = 1:size(traces, 2)
    y_norm = double(norm_traces(:, i));
    valid = isfinite(norm_energy_axis) & isfinite(y_norm);
    area = NaN;
    if nnz(valid) >= 2
        area = trapz(norm_energy_axis(valid), y_norm(valid));
    end
    if ~isfinite(area) || abs(area) <= eps
        if nnz(valid) >= 2
            area = trapz(norm_energy_axis(valid), abs(y_norm(valid)));
        end
    end
    if ~isfinite(area) || abs(area) <= eps
        area = 1;
    end
    normalized(:, i) = double(traces(:, i)) ./ area;
end
scale_values = abs(normalized(isfinite(normalized)));
if isempty(scale_values)
    return
end
display_scale = prctile(scale_values, 95);
if isfinite(display_scale) && display_scale > eps
    normalized = normalized .* (0.85 ./ display_scale);
end
end


function traces = local_apply_waterfall_residual_and_gain(energy_axis, ...
    traces, opts)
if isfield(opts, 'waterfall_residual') && opts.waterfall_residual
    for i = 1:size(traces, 2)
        y = double(traces(:, i));
        baseline = local_waterfall_asls_baseline(y, 1e6, 0.01, 12);
        traces(:, i) = y - baseline;
    end
    scale_values = abs(traces(isfinite(traces)));
    if ~isempty(scale_values)
        scale = prctile(scale_values, 95);
        if isfinite(scale) && scale > eps
            traces = traces .* (0.85 ./ scale);
        end
    end
end

if isfield(opts, 'waterfall_gain') && isfinite(opts.waterfall_gain) && ...
        opts.waterfall_gain > 0
    traces = traces .* opts.waterfall_gain;
end
end


function baseline = local_waterfall_asls_baseline(y, lambda_value, ...
    asymmetry, iterations)
y = double(y(:));
n = numel(y);
if n < 3
    baseline = y;
    return
end
lambda_value = max(double(lambda_value), 1);
asymmetry = min(max(double(asymmetry), 1e-4), 0.49);
iterations = max(1, round(iterations));
d = diff(speye(n), 2);
penalty = lambda_value * (d' * d);
weights = ones(n, 1);
for iter = 1:iterations %#ok<NASGU>
    w = spdiags(weights, 0, n, n);
    baseline = (w + penalty) \ (weights .* y);
    weights = asymmetry * (y > baseline) + (1 - asymmetry) * (y <= baseline);
end
end


function z = local_waterfall_soft_limit(z, lower_limit, upper_knee, extra_span)
z = max(z, lower_limit);
over = z > upper_knee;
if any(over)
    excess = z(over) - upper_knee;
    max_excess = max(excess, [], 'omitnan');
    if isfinite(max_excess) && max_excess > eps
        z(over) = upper_knee + extra_span .* log1p(excess) ./ log1p(max_excess);
    else
        z(over) = upper_knee;
    end
end
end


function color = local_waterfall_color(index, n_traces)
if n_traces <= 1
    t = 0.5;
else
    t = (index - 1) / (n_traces - 1);
end
red = [0.92 0.02 0.00];
black = [0.02 0.02 0.02];
green = [0.00 0.85 0.08];
if t <= 0.5
    u = t / 0.5;
    color = (1 - u) .* red + u .* black;
else
    u = (t - 0.5) / 0.5;
    color = (1 - u) .* black + u .* green;
end
end


function [energy_axis, energy_mask] = local_stack_energy_window(qe, opts)
energy_axis = double(qe.energy_meV(:));
margin = 150;
energy_mask = energy_axis >= opts.energy_window_meV(1) - margin & ...
    energy_axis <= opts.energy_window_meV(2) + margin;
energy_axis = energy_axis(energy_mask);
end


function [q_values, traces, source_counts] = local_stack_trace_set(qe, extract, opts, ...
    energy_mask, mode)
q_axis = double(qe.q_Ainv(:));

switch char(mode)
    case 'positive_q'
        groups = local_binning_map_q_groups(extract, q_axis, 'positive');
        [q_values, traces, source_counts] = local_traces_from_q_groups( ...
            qe, energy_mask, groups, false);
    case 'negative_q'
        groups = local_binning_map_q_groups(extract, q_axis, 'negative');
        [q_values, traces, source_counts] = local_traces_from_q_groups( ...
            qe, energy_mask, groups, false);
    case 'signed_q'
        groups = local_binning_map_q_groups(extract, q_axis, 'signed');
        [q_values, traces, source_counts] = local_traces_from_q_groups( ...
            qe, energy_mask, groups, false);
    case 'absq_combined'
        error('run_b1_double_peak_binning_analysis:DeprecatedAbsQCombined', ...
            ['absq_combined stacked spectra are obsolete because they ', ...
            'average +q and -q. Use signed_q instead.']);
    otherwise
        q_values = zeros(0, 1);
        traces = zeros(nnz(energy_mask), 0);
        source_counts = zeros(0, 1);
end

q_values = double(q_values(:));
source_counts = double(source_counts(:));
end


function [q_values, traces, source_counts] = local_traces_from_q_groups( ...
    qe, energy_mask, groups, use_abs_q)
traces = zeros(nnz(energy_mask), 0);
q_values = zeros(0, 1);
source_counts = zeros(0, 1);
q_axis = double(qe.q_Ainv(:));

for i = 1:numel(groups)
    if isstruct(groups)
        q_idx = groups(i).indices;
    else
        q_idx = groups{i};
    end
    q_idx = q_idx(q_idx >= 1 & q_idx <= numel(q_axis));
    if isempty(q_idx)
        continue
    end
    traces(:, end + 1) = mean(double(qe.intensity(energy_mask, q_idx)), ...
        2, 'omitnan'); %#ok<AGROW>
    if use_abs_q
        q_values(end + 1, 1) = mean(abs(q_axis(q_idx)), 'omitnan'); %#ok<AGROW>
    else
        q_values(end + 1, 1) = mean(q_axis(q_idx), 'omitnan'); %#ok<AGROW>
    end
    source_counts(end + 1, 1) = numel(q_idx); %#ok<AGROW>
end
end


function groups = local_binning_map_q_groups(extract, q_axis, mode)
groups = local_empty_q_group();
groups = groups([]);
if ~isfield(extract, 'binning_map') || isempty(extract.binning_map)
    return
end
map = extract.binning_map;
if height(map) == 0
    return
end

for i = 1:height(map)
    q_idx = local_parse_index_list(local_table_text(map.source_q_index, i));
    q_idx = q_idx(q_idx >= 1 & q_idx <= numel(q_axis));
    if isempty(q_idx)
        continue
    end
    q_mean = mean(q_axis(q_idx), 'omitnan');
    keep = true;
    switch char(mode)
        case 'positive'
            keep = q_mean > 0;
        case 'negative'
            keep = q_mean < 0;
        case 'combined_only'
            keep = startsWith(local_table_text(map.source_mode, i), ...
                'combined_q_binning');
        case 'signed'
            keep = true;
        otherwise
            keep = true;
    end
    if ~keep
        continue
    end

    item = local_empty_q_group();
    item.indices = q_idx;
    item.source_mode = local_table_text(map.source_mode, i);
    item.q_mean = q_mean;
    item.q_abs_mean = mean(abs(q_axis(q_idx)), 'omitnan');
    groups(end + 1) = item; %#ok<AGROW>
end

if isempty(groups)
    return
end
if strcmp(char(mode), 'negative')
    [~, order] = sort([groups.q_mean], 'ascend');
elseif strcmp(char(mode), 'positive')
    [~, order] = sort([groups.q_mean], 'ascend');
elseif strcmp(char(mode), 'signed')
    [~, order] = sort([groups.q_mean], 'ascend');
else
    [~, order] = sort([groups.q_abs_mean], 'ascend');
end
groups = groups(order);
end


function groups = local_absq_combined_noise_q_groups(extract, q_axis)
pos_groups = local_binning_map_q_groups(extract, q_axis, 'positive');
neg_groups = local_binning_map_q_groups(extract, q_axis, 'negative');
pos_groups = local_sort_struct_groups_by_abs(pos_groups);
neg_groups = local_sort_struct_groups_by_abs(neg_groups);

n_groups = max(numel(pos_groups), numel(neg_groups));
groups = local_empty_q_group();
groups = groups([]);
for i = 1:n_groups
    q_idx = [];
    modes = {};
    if i <= numel(neg_groups)
        q_idx = [q_idx, neg_groups(i).indices]; %#ok<AGROW>
        modes{end + 1} = neg_groups(i).source_mode; %#ok<AGROW>
    end
    if i <= numel(pos_groups)
        q_idx = [q_idx, pos_groups(i).indices]; %#ok<AGROW>
        modes{end + 1} = pos_groups(i).source_mode; %#ok<AGROW>
    end
    if isempty(q_idx)
        continue
    end
    item = local_empty_q_group();
    item.indices = unique(q_idx, 'stable');
    item.source_mode = strjoin(modes, '|');
    item.q_mean = mean(q_axis(item.indices), 'omitnan');
    item.q_abs_mean = mean(abs(q_axis(item.indices)), 'omitnan');
    groups(end + 1) = item; %#ok<AGROW>
end
end


function groups = local_sort_struct_groups_by_abs(groups)
if isempty(groups)
    return
end
[~, order] = sort([groups.q_abs_mean], 'ascend');
groups = groups(order);
end


function item = local_empty_q_group()
item = struct('indices', [], 'source_mode', '', ...
    'q_mean', NaN, 'q_abs_mean', NaN);
end


function text = local_table_text(column, index)
if iscell(column)
    text = char(column{index});
else
    text = char(string(column(index)));
end
end


function local_plot_stack_figure(energy_axis, traces, q_values, source_counts, ...
    old_points, mode, title_text, y_label, png_path, pdf_path)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1120 760]);
ax = axes(fig);

if isempty(traces) || size(traces, 2) == 0
    text(ax, 0.5, 0.5, 'No spectra available in the requested q range', ...
        'HorizontalAlignment', 'center');
    axis(ax, 'off');
else
    offset = qe_auto_stack_offset(traces, 1.25);
    offsets = (0:(size(traces, 2) - 1)) .* offset;
    colors = parula(max(size(traces, 2), 2));
    hold(ax, 'on');
    for i = 1:size(traces, 2)
        plot(ax, energy_axis, traces(:, i) + offsets(i), '-', ...
            'Color', colors(i, :), 'LineWidth', 0.9);
        if ~isempty(old_points)
            old_energy = local_old_point_energy(old_points, q_values(i), mode);
            if isfinite(old_energy)
                marker_y = interp1(energy_axis, traces(:, i), old_energy, ...
                    'linear', NaN) + offsets(i);
                if isfinite(marker_y)
                    plot(ax, old_energy, marker_y, 'o', ...
                        'MarkerSize', 4.5, 'MarkerFaceColor', [0.25 0.25 0.25], ...
                        'MarkerEdgeColor', 'w', 'LineWidth', 0.35);
                end
            end
        end
    end
    hold(ax, 'off');
    tick_idx = unique(round(linspace(1, numel(q_values), ...
        min(numel(q_values), 11))));
    set(ax, 'YTick', offsets(tick_idx), ...
        'YTickLabel', local_q_tick_labels(q_values(tick_idx), ...
        source_counts(tick_idx)));
    ylabel(ax, y_label);
end

box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.14;
xlabel(ax, 'Energy (meV)');
title(ax, title_text);
subtitle(ax, 'Displayed from existing qe_pp spectra; no added smoothing, denoise, or filtering');
exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function labels = local_q_tick_labels(q_values, source_counts)
labels = strings(numel(q_values), 1);
for i = 1:numel(q_values)
    if source_counts(i) > 1
        labels(i) = sprintf('%.4f (%d q)', q_values(i), source_counts(i));
    else
        labels(i) = sprintf('%.4f', q_values(i));
    end
end
labels = cellstr(labels);
end


function energy = local_old_point_energy(old_points, q_value, mode)
energy = NaN;
if isempty(old_points) || ~all(ismember({'q_Ainv', 'energy_meV'}, ...
        old_points.Properties.VariableNames))
    return
end
switch char(mode)
    case 'absq_combined'
        distance = abs(abs(old_points.q_Ainv) - abs(q_value));
    case 'combined_q_binning_3'
        distance = abs(abs(old_points.q_Ainv) - abs(q_value));
    otherwise
        distance = abs(old_points.q_Ainv - q_value);
end
[min_dist, idx] = min(distance);
if ~isempty(idx) && isfinite(min_dist) && min_dist <= 0.006
    energy = old_points.energy_meV(idx);
end
end


function seed = local_manual_window_seed_table(session, opts)
q_bounds = sort(double(opts.q_range_Ainv(:)).');
q_abs_max = min(max(abs(q_bounds)), 0.15);
base_edges = [0.005 0.030; 0.030 0.060; 0.060 0.100; 0.100 0.150];
rows = {};
for i = 1:size(base_edges, 1)
    q_min = base_edges(i, 1);
    q_max = min(base_edges(i, 2), q_abs_max);
    if q_max <= q_min
        continue
    end
    rows(end + 1, :) = {char(session.session_key), char(session.session_label), ...
        q_min, q_max, NaN, 150, NaN, NaN, NaN, 150, NaN, NaN, ...
        'fill from stacked spectra before constrained re-extraction'}; %#ok<AGROW>
end
if isempty(rows)
    rows = {char(session.session_key), char(session.session_label), ...
        NaN, NaN, NaN, 150, NaN, NaN, NaN, 150, NaN, NaN, ...
        'no q bins available under requested range'};
end
seed = cell2table(rows, 'VariableNames', { ...
    'session_key', 'session_label', 'q_abs_min_Ainv', 'q_abs_max_Ainv', ...
    'lower_center_meV', 'lower_half_width_meV', ...
    'lower_window_min_meV', 'lower_window_max_meV', ...
    'upper_center_meV', 'upper_half_width_meV', ...
    'upper_window_min_meV', 'upper_window_max_meV', 'notes'});
end


function plot_map = local_plot_q_binning_map(qe, extract, opts)
q_axis = double(qe.q_Ainv(:));
modes = {'stacked_signed_q', 'waterfall_signed_q'};
rows = {};

for mode_idx = 1:numel(modes)
    mode = modes{mode_idx};
    groups = local_plot_q_groups(q_axis, extract, mode);
    for i = 1:numel(groups)
        if isstruct(groups)
            q_idx = groups(i).indices;
            source_mode = groups(i).source_mode;
        else
            q_idx = groups{i};
            source_mode = '';
        end
        q_idx = q_idx(q_idx >= 1 & q_idx <= numel(q_axis));
        if isempty(q_idx)
            continue
        end
        rows(end + 1, :) = {mode, i, mean(q_axis(q_idx), 'omitnan'), ...
            mean(abs(q_axis(q_idx)), 'omitnan'), numel(q_idx), ...
            source_mode, ...
            local_join_numeric(q_axis(q_idx), '%.6g'), ...
            local_join_numeric(q_idx, '%d')}; %#ok<AGROW>
    end
end

if isempty(rows)
    rows = {'none', NaN, NaN, NaN, 0, '', '', ''};
end
plot_map = cell2table(rows, 'VariableNames', {'plot_mode', 'plot_bin_index', ...
    'q_Ainv_mean', 'q_abs_Ainv_mean', 'source_q_count', ...
    'source_mode', 'source_q_Ainv', 'source_q_index'});
end


function groups = local_plot_q_groups(q_axis, extract, mode)
switch char(mode)
    case 'stacked_positive_q'
        groups = local_binning_map_q_groups(extract, q_axis, 'positive');
    case 'stacked_negative_q'
        groups = local_binning_map_q_groups(extract, q_axis, 'negative');
    case {'stacked_signed_q', 'waterfall_signed_q'}
        groups = local_binning_map_q_groups(extract, q_axis, 'signed');
    case 'waterfall_combined_q'
        groups = local_binning_map_q_groups(extract, q_axis, 'combined_only');
    otherwise
        groups = local_empty_q_group();
        groups = groups([]);
end
end


function text = local_join_numeric(values, format_spec)
if isempty(values)
    text = '';
    return
end
parts = compose(format_spec, values(:));
text = char(strjoin(parts, '|'));
end


function local_plot_double_peak_heatmap(qe, old_points, extract, opts, ...
    session, png_path, pdf_path)
q_axis = double(qe.q_Ainv(:));
energy_axis = double(qe.energy_meV(:));
q_mask = q_axis >= opts.q_range_Ainv(1) & q_axis <= opts.q_range_Ainv(2);
e_margin = 150;
e_mask = energy_axis >= opts.energy_window_meV(1) - e_margin & ...
    energy_axis <= opts.energy_window_meV(2) + e_margin;
map = double(qe.intensity(e_mask, q_mask));

fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1120 680]);
ax = axes(fig);
imagesc(ax, q_axis(q_mask), energy_axis(e_mask), map);
axis(ax, 'xy');
colormap(ax, turbo);
clim_vals = local_color_limits(map);
if all(isfinite(clim_vals)) && clim_vals(1) < clim_vals(2)
    clim(ax, clim_vals);
end
colorbar(ax);
hold(ax, 'on');

scatter(ax, old_points.q_Ainv, old_points.energy_meV, 15, ...
    [0.85 0.85 0.85], 'filled', ...
    'MarkerFaceAlpha', 0.45, 'MarkerEdgeAlpha', 0.45, ...
    'DisplayName', 'old B1 single');

local_overlay_branch_source(ax, extract.lower_points, ...
    'b1_double_peak_lower', [0.05 0.70 1.00]);
local_overlay_branch_source(ax, extract.upper_points, ...
    'b1_double_peak_upper', [1.00 0.20 0.55]);

hold(ax, 'off');
box(ax, 'on');
grid(ax, 'on');
ax.GridAlpha = 0.16;
xlabel(ax, 'q (A^{-1})');
ylabel(ax, 'Energy (meV)');
title(ax, sprintf('B1 mandatory double-peak positions on q-E map: %s', ...
    session.session_label));
legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
    'Box', 'off');

exportgraphics(fig, png_path, 'Resolution', 300);
exportgraphics(fig, pdf_path, 'ContentType', 'vector');
close(fig);
end


function local_overlay_branch_source(ax, points, label, color)
if isempty(points)
    return
end
direct = strcmp(points.source_mode, 'single_q_direct');
combined = startsWith(points.source_mode, 'combined_q_binning');
if any(direct)
    scatter(ax, points.q_Ainv(direct), points.energy_meV(direct), ...
        30, color, 'filled', 'o', ...
        'MarkerEdgeColor', 'w', 'LineWidth', 0.45, ...
        'DisplayName', sprintf('%s direct-q', label));
end
if any(combined)
    scatter(ax, points.q_Ainv(combined), points.energy_meV(combined), ...
        42, color, 's', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.65, ...
        'DisplayName', sprintf('%s combined-q', label));
end
end


function clim_vals = local_color_limits(map)
vals = map(isfinite(map));
if isempty(vals)
    clim_vals = [0 1];
    return
end
lo = prctile(vals, 2);
hi = prctile(vals, 99);
if ~isfinite(lo) || ~isfinite(hi) || lo >= hi
    lo = min(vals);
    hi = max(vals);
end
if lo == hi
    hi = lo + eps;
end
clim_vals = [lo hi];
end


function values = local_parse_index_list(text)
parts = strsplit(char(text), '|');
values = str2double(parts);
values = values(isfinite(values));
values = max(1, round(values(:).'));
end
